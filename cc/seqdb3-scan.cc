#include <map>
#include <vector>
#include <algorithm>

#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/filesystem.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string-split.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<bool> all_lab_messages{*this, "all-lab-messages", desc{"otherwise show messages for WHO CCs only"}};
    option<bool> all_subtypes_messages{*this, "all-subtypes-messages", desc{"otherwise show messages for H1, H3, B only"}};

    argument<str_array> filenames{*this, arg_name{"filename"}, mandatory};
};

static int report(const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const Options& opt);

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        auto all_sequences = acmacs::seqdb::fasta::scan(opt.filenames, {});

#pragma omp parallel for default(shared) schedule(static, 256)
        for (size_t e_no = 0; e_no < all_sequences.size(); ++e_no) {
            auto& entry = all_sequences[e_no];
            entry.seq.sequence.translate();
            entry.aligned = entry.seq.sequence.align(entry.seq.type_subtype, entry.seq.fasta_name);
        }

        // for (const auto& seq_e : all_sequences) {
        //     if (seq_e.seq.type_subtype == "A(H3N2)" && !seq_e.seq.sequence.aa().empty()
        //         // && seq_e.seq.sequence.aa().find("ATLCLG") > 50 && seq_e.seq.sequence.aa().find("AMLCLG") > 50
        //         && seq_e.seq.name->find("SINGAPORE/INFIMH-16-0019") != std::string::npos
        //         ) {
        //         // fmt::print(stderr, "{}:{}: {}\n    {}\n  {}\n", seq_e.filename, seq_e.line_no, seq_e.seq.fasta_name, seq_e.seq.sequence.nuc(), seq_e.seq.sequence.aa());
        //         fmt::print(stderr, "{}:{}: {}\n {}\n", seq_e.filename, seq_e.line_no, seq_e.seq.fasta_name, seq_e.seq.sequence.aa());
        //     }
        // }

        const auto errors = report(all_sequences, opt);

        return errors;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------

static inline bool whocc_lab(std::string_view lab)
{
    return lab == "CDC" || lab == "Crick" || lab == "NIID" || lab == "VIDRL"; // || lab == "CNIC"
}

static inline bool our_subtype(std::string_view type_subtype)
{
    return type_subtype == "B" || type_subtype == "A(H1N1)" || type_subtype == "A(H3N2)";
}

template <typename Key> static inline std::vector<std::pair<Key, size_t>> sorted_by_count(const std::map<Key, size_t>& source)
{
    std::vector<std::pair<Key, size_t>> result(std::begin(source), std::end(source));
    std::sort(result.begin(), result.end(), [](const auto& e1, const auto& e2) { return e1.second > e2.second; });
    return result;
}

// ----------------------------------------------------------------------

int report(const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const Options& opt)
{
    std::map<std::string, size_t> location_not_found;
    std::map<std::string, size_t> unrecognized_passage;
    std::map<std::string, size_t> labs;
    std::map<std::string, size_t> subtypes;
    std::map<std::string, std::map<size_t, size_t>> subtypes_sequence_length;
    std::map<std::string, size_t> lineages;
    int errors = 0;

    const auto update = [](auto& counter, const auto& source) {
        if (auto [iter, inserted] = counter.emplace(source, 1UL); !inserted)
            ++iter->second;
    };

    for (const auto& entry : sequences) {
        // if (!entry.seq.lab.empty())
        update(labs, entry.seq.lab);

        update(subtypes, entry.seq.type_subtype);
        if (entry.seq.type_subtype.empty())
            fmt::print(stderr, "{}:{}: No subtype for {}\n", entry.filename, entry.line_no, *entry.seq.name);
        update(subtypes_sequence_length.emplace(entry.seq.type_subtype, std::map<size_t, size_t>{}).first->second, entry.seq.sequence.nuc().size());

        if (!entry.seq.lineage.empty())
            update(lineages, entry.seq.lineage);

        if ((opt.all_lab_messages || whocc_lab(entry.seq.lab)) && (opt.all_subtypes_messages || our_subtype(entry.seq.type_subtype))) {
            for (const auto& msg : entry.messages) {
                if (msg == acmacs::virus::parse_result_t::message_t::location_not_found) {
                    if (msg.value == "CRIE")
                        fmt::print(stderr, "{}:{}: [CRIE] {}\n", entry.filename, entry.line_no, msg);
                    update(location_not_found, msg.value);
                }
                else if (msg == acmacs::virus::parse_result_t::message_t::unrecognized_passage)
                    update(unrecognized_passage, msg.value);
                else {
                    fmt::print(stderr, "{}:{}: {} --> {}\n", entry.filename, entry.line_no, msg, *entry.seq.name);
                    ++errors;
                }
            }
        }
    }

    const auto report_by_count = [](const std::map<std::string, size_t>& source, const char* title) {
        fmt::print(stderr, "{}: {}\n", title, source.size());
        for (const auto& entry : sorted_by_count(source))
            fmt::print(stderr, "{:6d} {}\n", entry.second, entry.first);
        fmt::print(stderr, "\n");
    };

    if (!unrecognized_passage.empty()) {
        report_by_count(unrecognized_passage, "Unrecognized PASSAGE");
        ++errors;
    }

    if (!location_not_found.empty()) {
        report_by_count(location_not_found, "Not found LOCATION");
        ++errors;
    }

    fmt::print(stderr, "======================================================================\n\n");
    fmt::print(stderr, "TOTAL: {}\n\n", sequences.size());
    report_by_count(subtypes, "SUBTYPES");
    report_by_count(lineages, "LINEAGES");
    report_by_count(labs, "LABS");

    fmt::print(stderr, "SUBTYPES and sequence lengths (count:seq-length)\n");
    for (const auto& [subtype, entry] : subtypes_sequence_length) {
        if (subtypes.find(subtype)->second > 1000) {
            fmt::print(stderr, "  {:<10s}", subtype);
            size_t cnt = 0;
            for (const auto& entry : sorted_by_count(entry)) {
                if (entry.second < 100)
                    break;
                if (cnt > 14) {
                    fmt::print(stderr, "\n            ");
                    cnt = 0;
                }
                fmt::print(stderr, " {:6d}:{:4d}", entry.second, entry.first);
                ++cnt;
            }
            fmt::print(stderr, "\n");
        }
    }
    fmt::print(stderr, "\n");

    return errors;

} // report

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
