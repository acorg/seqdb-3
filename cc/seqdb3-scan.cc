#include <iostream>
#include <set>
#include <map>

#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/enumerate.hh"
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

    argument<str_array> filenames{*this, arg_name{"filename"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        get_locdb();            // load locbd outside of threading code, it is not thread safe
        // acmacs::virus::parse_passage("E1"); // init parse_passage static vars, to avoid conflict in multithreading

        std::vector<std::vector<scan_result_t>> sequences_per_file(opt.filenames->size());
// #pragma omp parallel for default(shared) schedule(static, 4)
        for (size_t f_no = 0; f_no < opt.filenames->size(); ++f_no) {
            const auto& filename = (*opt.filenames)[f_no];
            try {
                const std::string file_data_s = acmacs::file::read(filename);
                const std::string_view file_data = file_data_s;
                acmacs::seqdb::v3::fasta::scan_input_t file_input{std::begin(file_data), std::end(file_data)};
                while (!file_input.done()) {
                    acmacs::seqdb::v3::fasta::scan_output_t sequence_ref;
                    std::tie(file_input, sequence_ref) = acmacs::seqdb::v3::fasta::scan(file_input);

                    std::optional<acmacs::seqdb::v3::fasta::sequence_t> seq;
                    for (auto parser : {&acmacs::seqdb::v3::fasta::name_gisaid_spaces, &acmacs::seqdb::v3::fasta::name_gisaid_underscores, &acmacs::seqdb::v3::fasta::name_plain}) {
                        seq = (*parser)(sequence_ref.name, filename, file_input.name_line_no);
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
                        std::cerr << "WARNING: " << filename << ':' << file_input.name_line_no << ": unable to parse fasta name: " << sequence_ref.name << '\n';
                }
            }
            catch (std::exception& err) {
                throw std::runtime_error(string::concat(filename, err.what()));
            }
        }

        int errors = 0;
        std::set<std::string> location_not_found;
        std::map<std::string, size_t> unrecognized_passage;
        for (const auto& per_file : sequences_per_file) {
            for (const auto& entry : per_file) {
                for (const auto& msg : entry.messages) {
                    if (msg == acmacs::virus::parse_result_t::message_t::location_not_found) {
                        // if (msg.value == "CHU")
                        //     std::cerr << entry.filename << ':' << entry.line_no << ": " << msg << '\n';
                        location_not_found.insert(msg.value);
                    }
                    else if (msg == acmacs::virus::parse_result_t::message_t::unrecognized_passage) {
                        if (auto [iter, inserted] = unrecognized_passage.emplace(msg.value, 1UL); !inserted)
                            ++iter->second;
                    }
                    else {
                        std::cerr << entry.filename << ':' << entry.line_no << ": " << msg << '\n';
                        ++errors;
                    }
                }
            }
        }

        if (!unrecognized_passage.empty()) {
            std::cerr << "Unrecognized PASSAGE " << unrecognized_passage.size() << '\n';
            for (const auto& entry : unrecognized_passage)
                std::cerr << "  " << std::setw(3) << std::right << entry.second << ' ' << entry.first << '\n';
            ++errors;
        }

        if (!location_not_found.empty()) {
            std::cerr << "LOCATION NOT FOUND " << location_not_found.size() << '\n';
            for (const auto& name : location_not_found)
                std::cerr << "  " << name << '\n';
            ++errors;
        }

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
