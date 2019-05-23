#include <map>
#include <vector>
#include <algorithm>

#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/filesystem.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string-split.hh"
#include "locationdb/locdb.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

struct scan_result_t
{
    acmacs::seqdb::v3::fasta::sequence_t seq;
    std::vector<acmacs::virus::v2::parse_result_t::message_t> messages;
    std::string_view filename;
    size_t line_no;
};

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<bool> all_lab_messages{*this, "all-lab-messages", desc{"otherwise show messages for WHO CCs only"}};

    argument<str_array> filenames{*this, arg_name{"filename"}, mandatory};
};

static inline acmacs::seqdb::fasta::hint_t find_hints(std::string_view filename)
{
    const auto stem = fs::path{filename}.stem().stem().string();
    const auto fields = acmacs::string::split(stem, "-");
    acmacs::seqdb::fasta::hint_t hints;
    hints.lab = ::string::upper(fields[0]);
    if (fields[1] == "h1pdm" || fields[1] == "h1")
        hints.subtype = "H1N1";
    else if (fields[1] == "h3")
        hints.subtype = "H3N2";
    else if (fields[1] == "b" && fields[0] == "niid" && fields.size() == 4) {
        if (fields[3] == "vic")
            hints.lineage = "VICTOIRA";
        else if (fields[3] == "yam")
            hints.lineage = "YAMAGATA";
    }
    return hints;
}

static inline bool whocc_lab(std::string_view lab)
{
    return lab == "CDC" || lab == "Crick" || lab == "NIID" || lab == "VIDRL"; // || lab == "CNIC"
}

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        get_locdb(); // load locbd outside of threading code, it is not thread safe

        std::vector<std::vector<scan_result_t>> sequences_per_file(opt.filenames->size());
#pragma omp parallel for default(shared) schedule(static, 4)
        for (size_t f_no = 0; f_no < opt.filenames->size(); ++f_no) {
            const auto& filename = (*opt.filenames)[f_no];
            const auto hints = find_hints(filename);
            try {
                const std::string file_data_s = acmacs::file::read(filename);
                const std::string_view file_data = file_data_s;
                acmacs::seqdb::v3::fasta::scan_input_t file_input{std::begin(file_data), std::end(file_data)};
                while (!file_input.done()) {
                    acmacs::seqdb::v3::fasta::scan_output_t sequence_ref;
                    std::tie(file_input, sequence_ref) = acmacs::seqdb::v3::fasta::scan(file_input);

                    std::optional<acmacs::seqdb::v3::fasta::sequence_t> seq;
                    for (auto parser : {&acmacs::seqdb::v3::fasta::name_gisaid_spaces, &acmacs::seqdb::v3::fasta::name_gisaid_underscores, &acmacs::seqdb::v3::fasta::name_plain}) {
                        seq = (*parser)(sequence_ref.name, hints, filename, file_input.name_line_no);
                        if (seq.has_value())
                            break;
                    }
                    if (seq.has_value()) {
                        const auto messages = acmacs::seqdb::v3::fasta::normalize_name(*seq);
                        // for (const auto& msg : messages)
                        //     std::cerr << "WARNING: " << filename << ':' << file_input.name_line_no << ": " << msg << '\n';
                        seq->sequence = acmacs::seqdb::v3::fasta::normalize_sequence(sequence_ref.sequence);
                        sequences_per_file[f_no].push_back({*seq, messages, filename, file_input.name_line_no});
                    }
                    else
                        fmt::print(stderr, "WARNING: {}:{}: unable to parse fasta name: {}\n", filename, file_input.name_line_no, sequence_ref.name);
                }
            }
            catch (std::exception& err) {
                throw std::runtime_error(string::concat(filename, ": ", err.what()));
            }
        }

        int errors = 0;
        std::map<std::string, size_t> location_not_found;
        std::map<std::string, size_t> unrecognized_passage;
        std::map<std::string, size_t> labs;
        for (const auto& per_file : sequences_per_file) {
            for (const auto& entry : per_file) {
                if (!entry.seq.lab.empty()) {
                    if (auto [iter, inserted] = labs.emplace(entry.seq.lab, 1UL); !inserted)
                        ++iter->second;
                }
                if (opt.all_lab_messages || whocc_lab(entry.seq.lab)) {
                    for (const auto& msg : entry.messages) {
                        if (msg == acmacs::virus::parse_result_t::message_t::location_not_found) {
                            if (msg.value == "CRIE")
                                fmt::print(stderr, "{}:{}: [CRIE] {}\n", entry.filename, entry.line_no, msg);
                            if (auto [iter, inserted] = location_not_found.emplace(msg.value, 1UL); !inserted)
                                ++iter->second;
                        }
                        else if (msg == acmacs::virus::parse_result_t::message_t::unrecognized_passage) {
                            if (auto [iter, inserted] = unrecognized_passage.emplace(msg.value, 1UL); !inserted)
                                ++iter->second;
                        }
                        else {
                            fmt::print(stderr, "{}:{}: {} --> {}\n", entry.filename, entry.line_no, msg, *entry.seq.name);
                            ++errors;
                        }
                    }
                }
            }
        }

        if (!unrecognized_passage.empty()) {
            fmt::print(stderr, "Unrecognized PASSAGE {}\n", unrecognized_passage.size());
            for (const auto& entry : unrecognized_passage)
                fmt::print(stderr, "  {:3d} {}\n", entry.second, entry.first);
            ++errors;
        }

        if (!location_not_found.empty()) {
            fmt::print(stderr, "LOCATION NOT FOUND {}\n", location_not_found.size());
            for (const auto& entry : location_not_found)
                fmt::print(stderr, "  {:3d} {}\n", entry.second, entry.first);
            ++errors;
        }

        fmt::print(stderr, "\nLABS {}\n", labs.size());
        std::vector<std::pair<std::string, size_t>> labs_entries(std::begin(labs), std::end(labs));
        std::sort(labs_entries.begin(), labs_entries.end(), [](const auto& e1, const auto& e2) { return e1.second > e2.second; });
        for (const auto& entry : labs_entries)
            fmt::print(stderr, "  {:3d} {}\n", entry.second, entry.first);

        // std::vector<acmacs::seqdb::v3::fasta::sequence_t> all_sequences;
        // for (auto& per_file : sequences_per_file)
        //     std::move(per_file.begin(), per_file.end(), std::back_inserter(all_sequences));
        // std::cout << all_sequences.size() << " sequences read\n";

        return errors;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
