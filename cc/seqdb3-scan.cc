#include <cstdio>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/filesystem.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/counter.hh"
#include "seqdb-3/hamming-distance.hh"
#include "seqdb-3/align.hh"
#include "seqdb-3/insertions.hh"
#include "seqdb-3/clades.hh"
#include "seqdb-3/match-hidb.hh"
#include "seqdb-3/create.hh"

// ----------------------------------------------------------------------

static inline bool whocc_lab(std::string_view lab)
{
    return lab == "CDC" || lab == "Crick" || lab == "NIID" || lab == "VIDRL"; // || lab == "CNIC"
}

static inline bool our_subtype(std::string_view type_subtype)
{
    return type_subtype == "B" || type_subtype == "A(H1N1)" || type_subtype == "A(H3N2)";
}

// ----------------------------------------------------------------------

// static inline std::string infer_regex(const std::vector<std::string>& sources)
// {
//     std::vector<std::string> letters(128);
//     for (const auto& src : sources) {
//         for (size_t pos = 0; pos < src.size(); ++pos) {
//             if (ranges::find(letters[pos], src[pos]) == ranges::end(letters[pos]))
//                 letters[pos].append(1, src[pos]);
//         }
//     }
//     std::string res;
//     for (const auto& let : letters) {
//         if (let.empty())
//             break;
//         if (let.size() == 1)
//             res.append(let);
//         else {
//             res.append(1, '[');
//             res.append(let);
//             res.append(1, ']');
//         }
//     }
//     return res;
// }

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<bool> all_lab_messages{*this, "all-lab-messages", desc{"otherwise show messages for WHO CCs only"}};
    option<bool> all_subtypes_messages{*this, "all-subtypes-messages", desc{"otherwise show messages for H1, H3, B only"}};

    option<str>  output_seqdb{*this, 'o', "seqdb", dflt{""}};

    option<str>  print_aa_for{*this, "print-aa-for", dflt{""}};
    option<str>  print_not_aligned_for{*this, "print-not-aligned-for", dflt{""}};
    option<str>  print_counter_for{*this, "print-counter-for", dflt{""}};
    option<str>  print_aligned_for{*this, "print-aligned-for", dflt{""}};
    option<bool> print_aa_sizes{*this, "print-aa-sizes"};
    option<bool> print_stat{*this, "stat"};

    argument<str_array> filenames{*this, arg_name{"filename"}, mandatory};
};

static int report(const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const Options& opt);

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        auto all_sequences = acmacs::seqdb::fasta::scan(opt.filenames, {});
        fmt::print(stderr, "TOTAL sequences upon scanning fasta: {:7d}\n", all_sequences.size());

        // // keep just A(H3
        // all_sequences.erase(std::remove_if(std::begin(all_sequences), std::end(all_sequences), [](const auto& e1) { return e1.fasta.type_subtype.substr(0, 4) != "A(H3"; }),
        // std::end(all_sequences)); fmt::print(stderr, "before aligned (H3 only): {}\n", all_sequences.size());

        acmacs::seqdb::fasta::sort_by_date(all_sequences);
        acmacs::seqdb::translate_align(all_sequences);
        acmacs::seqdb::detect_insertions_deletions(all_sequences);
        acmacs::seqdb::detect_lineages_clades(all_sequences);
        acmacs::seqdb::fasta::sort_by_name(all_sequences);
        acmacs::seqdb::v3::match_hidb(all_sequences);
        if (!opt.output_seqdb->empty())
            acmacs::seqdb::create(opt.output_seqdb, all_sequences, acmacs::seqdb::create_filter::h1_h3_b_aligned);

        fmt::print(stderr, "TOTAL sequences upon translating:    {:7d}  aligned: {}\n", all_sequences.size(), ranges::count_if(all_sequences, acmacs::seqdb::fasta::is_aligned));
        fmt::print(stderr, "\n");

        if (!opt.print_counter_for->empty()) {
            const auto chunk = ::string::upper(*opt.print_counter_for);
            const auto found = [&chunk](size_t limit) { return [&chunk, limit](const auto& sc) { return sc.sequence.aa().find(std::string_view(chunk)) < limit; }; };
            for (auto limit : {50, 100, 150, 200, 1000}) {
                acmacs::Counter<std::string> counter;
                for (const auto& sc : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated) | ranges::view::filter(found(static_cast<size_t>(limit))))
                    counter.count(sc.fasta.type_subtype.h_or_b());
                fmt::print(stderr, "Counter for {} at first {} positions\n{}\n", chunk, limit, counter.report_sorted_max_first());
            }
        }

        if (const auto false_positive = acmacs::seqdb::fasta::report_false_positive(all_sequences, 200); !false_positive.empty())
            fmt::print(stderr, "FALSE POSITIVES {}\n{}\n", ranges::count(false_positive, '\n') / 2, false_positive);

        if (opt.print_counter_for->empty()) {
            acmacs::Counter<std::string> counter_not_aligned, counter_not_aligned_h;
            for (const auto& sc : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated) | ranges::view::filter(acmacs::seqdb::fasta::isnot_aligned)) {
                counter_not_aligned.count(*sc.fasta.type_subtype);
                counter_not_aligned_h.count(sc.fasta.type_subtype.h_or_b());
            }
            fmt::print(stderr, "NOT ALIGNED\n{}\n", counter_not_aligned_h.report_sorted_max_first());
            // fmt::print(stderr, "NOT ALIGNED\n{}\n", counter_not_aligned.report_sorted_max_first());
        }

        if (!opt.print_aa_for->empty()) {
            const auto report = acmacs::seqdb::fasta::report_aa(all_sequences, ::string::upper(*opt.print_aa_for), 99999);
            fmt::print("{} {}\n{}\n", *opt.print_aa_for, ranges::count(report, '\n') / 2, report);
        }

        if (!opt.print_not_aligned_for->empty()) {
            const auto report = acmacs::seqdb::fasta::report_not_aligned(all_sequences, ::string::upper(*opt.print_not_aligned_for), 200);
            fmt::print(stderr, "NOT ALIGNED {} {}\n{}\n", *opt.print_not_aligned_for, ranges::count(report, '\n'), report);
        }

        if (!opt.print_aligned_for->empty()) {
            const auto report = acmacs::seqdb::fasta::report_aa_aligned(all_sequences, ::string::upper(*opt.print_aligned_for));
            fmt::print("ALIGNED {} {}\n{}\n", *opt.print_aligned_for, ranges::count(report, '\n'), report);
        }

        if (opt.print_aa_sizes) {
            std::map<std::string, acmacs::Counter<size_t>> counter; // subtype -> size:count
            for (const auto& sc : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated) | ranges::view::filter(acmacs::seqdb::fasta::is_aligned))
                counter[std::string(sc.sequence.type_subtype().h_or_b())].count(sc.sequence.aa_aligned_length());
            fmt::print("AA sizes\n");
            for (const auto& [subtype, cntr] : counter)
                fmt::print("  {}\n{}\n", subtype, cntr.report_sorted_max_first("    "));
            fmt::print("\n");
        }

        if (opt.print_stat)
            report(all_sequences, opt);

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------

template <typename Key> static inline std::vector<std::pair<Key, size_t>> sorted_by_count(const std::map<Key, size_t>& source)
{
    std::vector<std::pair<Key, size_t>> result(std::begin(source), std::end(source));
    std::sort(result.begin(), result.end(), [](const auto& e1, const auto& e2) { return e1.second > e2.second; });
    return result;
}

// ----------------------------------------------------------------------

int report(const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const Options& opt)
{
    acmacs::Counter<std::string> location_not_found, unrecognized_passage, labs, subtypes, lineages;
    std::map<std::string, std::map<size_t, size_t>> subtypes_sequence_length;
    int errors = 0;

    // const auto update = [](auto& counter, const auto& source) {
    //     if (auto [iter, inserted] = counter.emplace(source, 1UL); !inserted)
    //         ++iter->second;
    // };

    for (const auto& entry : sequences) {
        labs.count(entry.sequence.lab());

        subtypes.count(entry.sequence.type_subtype());
        // if (entry.sequence.type_subtype().empty())
        //     fmt::print(stderr, "{}:{}: No subtype for {}\n", entry.fasta.filename, entry.fasta.line_no, *entry.sequence.name);
        if (auto [iter, inserted] = subtypes_sequence_length.emplace(entry.sequence.type_subtype(), std::map<size_t, size_t>{}).first->second.emplace(entry.sequence.nuc().size(), 1UL); !inserted)
            ++iter->second;

        lineages.count_if(!entry.sequence.lineage().empty(), entry.sequence.lineage());

        if ((opt.all_lab_messages || whocc_lab(entry.sequence.lab())) && (opt.all_subtypes_messages || our_subtype(entry.sequence.type_subtype()))) {
            for (const auto& msg : entry.fasta.messages) {
                if (msg == acmacs::virus::parse_result_t::message_t::location_not_found) {
                    if (msg.value == "CRIE")
                        fmt::print(stderr, "{}:{}: [CRIE] {}\n", entry.fasta.filename, entry.fasta.line_no, msg);
                    location_not_found.count(msg.value);
                }
                else if (msg == acmacs::virus::parse_result_t::message_t::unrecognized_passage)
                    unrecognized_passage.count(msg.value);
                else {
                    fmt::print(stderr, "{}:{}: {} --> {}\n", entry.fasta.filename, entry.fasta.line_no, msg, *entry.sequence.name());
                    ++errors;
                }
            }
        }
    }

    const auto report_by_count = [](const acmacs::Counter<std::string>& source, const char* title) {
        fmt::print(stderr, "{}: {}\n", title, source.size());
        for (const auto& entry : source.sorted_max_first())
            fmt::print(stderr, "{:6d} {}\n", entry->second, entry->first);
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
        if (subtypes[subtype] > 1000) {
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
