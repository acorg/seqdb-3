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


        const auto all_sequences = acmacs::seqdb::fasta::scan(opt.filenames, {});
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

static inline bool our_subtype(std::string_view subtype)
{
    return subtype.empty() || subtype == "H1N1" || subtype == "H3N2";
}

// ----------------------------------------------------------------------

int report(const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const Options& opt)
{
    std::map<std::string, size_t> location_not_found;
    std::map<std::string, size_t> unrecognized_passage;
    std::map<std::string, size_t> labs;
    int errors = 0;

    for (const auto& entry : sequences) {
        if (!entry.seq.lab.empty()) {
            if (auto [iter, inserted] = labs.emplace(entry.seq.lab, 1UL); !inserted)
                ++iter->second;
        }
        if ((opt.all_lab_messages || whocc_lab(entry.seq.lab)) && (opt.all_subtypes_messages || our_subtype(entry.seq.a_subtype))) {
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

    fmt::print(stderr, "\nTOTAL {}\n", sequences.size());

    fmt::print(stderr, "\nLABS {}\n", labs.size());
    std::vector<std::pair<std::string, size_t>> labs_entries(std::begin(labs), std::end(labs));
    std::sort(labs_entries.begin(), labs_entries.end(), [](const auto& e1, const auto& e2) { return e1.second > e2.second; });
    for (const auto& entry : labs_entries)
        fmt::print(stderr, "  {:3d} {}\n", entry.second, entry.first);

    return errors;

} // report

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
