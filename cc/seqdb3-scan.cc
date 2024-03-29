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
#include "acmacs-virus/virus-name-normalize.hh"
#include "seqdb-3/hamming-distance.hh"
#include "seqdb-3/scan-align.hh"
#include "seqdb-3/eliminate-identical.hh"
#include "seqdb-3/scan-deletions.hh"
#include "seqdb-3/scan-lineages.hh"
#include "seqdb-3/scan-match-hidb.hh"
#include "seqdb-3/hamming-distance-bins.hh"
#include "seqdb-3/create.hh"

// ----------------------------------------------------------------------

static inline bool whocc_lab(const acmacs::seqdb::scan::sequence_t& seq)
{
    return seq.lab_in({"CDC", "Crick", "NIID", "VIDRL"}); // "CNIC"
}

static inline bool our_subtype(std::string_view type_subtype)
{
    return type_subtype == "B" || type_subtype == "A(H1N1)" || type_subtype == "A(H3N2)";
}

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<bool> all_lab_messages{*this, "all-lab-messages", desc{"otherwise show messages for WHO CCs only"}};
    option<bool> all_subtypes_messages{*this, "all-subtypes-messages", desc{"otherwise show messages for H1, H3, B only"}};

    option<str>  output_seqdb{*this, 'o', "output-dir", dflt{""}};
    option<bool> whocc_only{*this, "whocc-only", desc{"create whocc only db (seqdb.json.xz)"}};
    option<bool> gisaid{*this, "gisaid", desc{"perform gisaid related name fixes and adjustments"}};
    option<str>  ncbi{*this, "ncbi", dflt{""}, desc{"directory with files downloaded from ncbi, see acmacs-whocc/doc/gisaid.org"}};
    option<bool> dont_eliminate_identical{*this, "dont-eliminate-identical", desc{"do not find identical sequences"}};

    option<str>  print_aa_for{*this, "print-aa-for", dflt{""}};
    option<str>  print_not_aligned_for{*this, "print-not-aligned-for", dflt{""}, desc{"ALL or comma separated: H1N,H3,B"}};
    // option<str>  print_counter_for{*this, "print-counter-for", dflt{""}};
    option<str>  print_aligned_for{*this, "print-aligned-for", dflt{""}};
    option<bool> print_aa_sizes{*this, "print-aa-sizes"};
    option<bool> print_stat{*this, "stat"};
    option<bool> print_names{*this, "print-names"};
    option<bool> print_messages{*this, 'm', "print-messages"};

    option<bool> verbose{*this, 'v', "verbose"};

    argument<str_array> filenames{*this, arg_name{"filename"}};
};

static int report(const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& sequences, const Options& opt);
static void report_messages(acmacs::messages::messages_t& messages);
static void report_issues(const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& all_sequences);

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        const acmacs::seqdb::scan::fasta::scan_options_t scan_options(opt.verbose ? acmacs::debug::yes : acmacs::debug::no,
                                                                      opt.gisaid ? acmacs::seqdb::scan::fasta::scan_name_adjustments::gisaid : acmacs::seqdb::scan::fasta::scan_name_adjustments::none,
                                                                      opt.print_names ? acmacs::seqdb::scan::fasta::print_names::yes : acmacs::seqdb::scan::fasta::print_names::no);
        acmacs::seqdb::scan::fasta::scan_results_t scan_results;
        if (!opt.filenames->empty())
            scan_results.merge(acmacs::seqdb::scan::fasta::scan(opt.filenames, scan_options));
        if (opt.ncbi)
            scan_results.merge(acmacs::seqdb::scan::fasta::scan_ncbi(opt.ncbi, scan_options));
        if (opt.print_messages)
            report_messages(scan_results.messages);

        // auto [all_sequences, messages] = acmacs::seqdb::scan::fasta::scan(opt.filenames, scan_options);

        auto& all_sequences = scan_results.results;
        if (all_sequences.empty())
            throw std::runtime_error("no sequences read (no files nor --ncbi in the command line?)");
        AD_INFO("Total sequences upon scanning fasta: {:7d}", all_sequences.size());
        acmacs::seqdb::scan::fasta::remove_without_names(all_sequences);
        acmacs::seqdb::scan::fasta::merge_duplicates(all_sequences);
        acmacs::seqdb::scan::translate_align(all_sequences);
        acmacs::seqdb::scan::detect_insertions_deletions(all_sequences);
        acmacs::seqdb::scan::detect_lineages_clades(all_sequences);
        // acmacs::seqdb::scan::fasta::sort_by_date(all_sequences);
        acmacs::seqdb::scan::match_hidb(all_sequences); // sorts all_sequences by name
        // acmacs::seqdb::scan::hamming_distance_bins_issues(all_sequences); // changes order of all_sequences
        if (!opt.dont_eliminate_identical)            // after hidb matching, because matching may change subtype (e.g. H3 -> H3N2) and it affectes reference to master
            acmacs::seqdb::scan::eliminate_identical(all_sequences); // changes order of all_sequences
        if (!opt.output_seqdb->empty()) {
            acmacs::seqdb::scan::fasta::sort_by_name(all_sequences);
            acmacs::seqdb::create(opt.output_seqdb, all_sequences, opt.whocc_only ? acmacs::seqdb::create_dbs::whocc_only : acmacs::seqdb::create_dbs::all);
        }

        AD_INFO("Total sequences upon translating:    {:7d}  aligned: {}", all_sequences.size(), ranges::count_if(all_sequences, acmacs::seqdb::scan::fasta::is_aligned));
        AD_PRINT("");

        // if (!opt.print_counter_for->empty()) {
        //     const auto chunk = ::string::upper(*opt.print_counter_for);
        //     const auto found = [&chunk](size_t limit) { return [&chunk, limit](const auto& sc) { return sc.sequence.aa().find(std::string_view(chunk)) < limit; }; };
        //     for (auto limit : {50, 100, 150, 200, 1000}) {
        //         acmacs::Counter<std::string> counter;
        //         for (const auto& sc : all_sequences | ranges::views::filter(acmacs::seqdb::scan::fasta::is_translated) | ranges::views::filter(found(static_cast<size_t>(limit))))
        //             counter.count(sc.fasta.type_subtype.h_or_b());
        //         fmt::print(stderr, "Counter for {} at first {} positions\n{}\n", chunk, limit, counter.report_sorted_max_first());
        //     }
        // }

        if (const auto false_positive = acmacs::seqdb::scan::fasta::report_false_positive(all_sequences, 200); !false_positive.empty()) {
            AD_ERROR("FALSE POSITIVES ({})", ranges::count(false_positive, '\n') / 2);
            fmt::print("{}\n", false_positive);
        }

        const auto dates_to_report = acmacs::seqdb::scan::fasta::min_max_dates(all_sequences);
        fmt::print(stderr, "Isolation date range:  {} .. {}\nSubmission date range: {} .. {}\n", dates_to_report.min_isolation_date, dates_to_report.max_isolation_date,
                   dates_to_report.min_submission_date, dates_to_report.max_submission_date);

        report_issues(all_sequences);

        if (!opt.print_aa_for->empty()) {
            const auto report = acmacs::seqdb::scan::fasta::report_aa(all_sequences, ::string::upper(*opt.print_aa_for), 99999);
            fmt::print("{} {}\n{}\n", *opt.print_aa_for, ranges::count(report, '\n') / 2, report);
        }

        if (!opt.print_not_aligned_for->empty()) {
            const auto report = acmacs::seqdb::scan::fasta::report_not_aligned(all_sequences, ::string::upper(*opt.print_not_aligned_for), 200);
            fmt::print(stderr, "NOT ALIGNED {} {} (name at the end)\n{}\n", *opt.print_not_aligned_for, ranges::count(report, '\n'), report);
        }

        if (!opt.print_aligned_for->empty()) {
            const auto report = acmacs::seqdb::scan::fasta::report_aa_aligned(all_sequences, ::string::upper(*opt.print_aligned_for));
            fmt::print("ALIGNED {} {}\n{}\n", *opt.print_aligned_for, ranges::count(report, '\n'), report);
        }

        if (opt.print_aa_sizes) {
            std::map<std::string, acmacs::Counter<size_t>> counter; // subtype -> size:count
            for (const auto& sc : all_sequences | ranges::views::filter(acmacs::seqdb::scan::fasta::is_translated) | ranges::views::filter(acmacs::seqdb::scan::fasta::is_aligned))
                counter[std::string(sc.sequence.type_subtype().h_or_b())].count(sc.sequence.aa_aligned_length());
            fmt::print("AA sizes\n");
            for (const auto& [subtype, cntr] : counter)
                fmt::print("  {}\n{}\n", subtype, cntr.report_sorted_max_first("    "));
            fmt::print("\n");
        }

        if (opt.print_stat)
            report(all_sequences, opt);

        if (opt.filenames->size() == 1 && opt.filenames->front().substr(0, 3) == "/r/") {
            fmt::print("mv -i {} /r/gisaid-{}-{}.fas", opt.filenames->front(), string::replace(dates_to_report.min_submission_date, "-", ""), string::replace(dates_to_report.max_submission_date, "-", ""));
        }

        return 0;
    }
    catch (std::exception& err) {
        AD_ERROR("{}", err);
        return 1;
    }
}

// ----------------------------------------------------------------------

template <typename Key> static inline std::vector<std::pair<Key, size_t>> sorted_by_count(const std::map<Key, size_t>& source)
{
    const auto order_by_value_reverse = [](const auto& e1, const auto& e2) { return e1.second > e2.second; };
    std::vector<std::pair<Key, size_t>> result(std::begin(source), std::end(source));
    std::sort(std::begin(result), std::end(result), order_by_value_reverse);
    return result;
}

// ----------------------------------------------------------------------

void report_messages(acmacs::messages::messages_t& messages)
{
    const auto index = acmacs::messages::make_index(messages);
    AD_INFO("Total messages: {}  keys: {}", messages.size(), index.size());
    for (const auto& [first, last] : index) {
        if (first != last) {
            if (first->key == acmacs::messages::key::location_not_found || first->key == acmacs::messages::key::unrecognized_passage)
                acmacs::messages::report_by_count(first, last);
            else
                acmacs::messages::report(first, last);
        }
    }

} // report_messages

// ----------------------------------------------------------------------

int report(const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& sequences, const Options& opt)
{
    acmacs::Counter<std::string> location_not_found, unrecognized_passage, labs, subtypes, lineages, clades, isolation_dates, submission_dates;
    std::map<std::string, std::map<size_t, size_t>> subtypes_sequence_length;
    int errors = 0;

    // const auto update = [](auto& counter, const auto& source) {
    //     if (auto [iter, inserted] = counter.emplace(source, 1UL); !inserted)
    //         ++iter->second;
    // };

    for (const auto& entry : sequences) {
        ranges::for_each(entry.sequence.lab_ids(), [&labs](const auto& en) { labs.count(en.first); });

        subtypes.count(entry.sequence.type_subtype());
        isolation_dates.count(entry.sequence.date_simulated().substr(0, 7)); // year-month
        if (!entry.sequence.gisaid_last_modified().empty())
            submission_dates.count(entry.sequence.gisaid_last_modified().front().substr(0, 7)); // year-month
        // if (entry.sequence.type_subtype().empty())
        //     fmt::print(stderr, "{}:{}: No subtype for {}\n", entry.fasta.filename, entry.fasta.line_no, *entry.sequence.name);
        if (auto [iter, inserted] = subtypes_sequence_length.emplace(entry.sequence.type_subtype(), std::map<size_t, size_t>{}).first->second.emplace(entry.sequence.nuc().size(), 1UL); !inserted)
            ++iter->second;

        lineages.count_if(!entry.sequence.lineage().empty(), entry.sequence.lineage());
        ranges::for_each(entry.sequence.clades(), [&clades](const auto& clade) { clades.count(clade); });

        if ((opt.all_lab_messages || whocc_lab(entry.sequence)) && (opt.all_subtypes_messages || our_subtype(entry.sequence.type_subtype()))) {
            for (const auto& msg : entry.fasta.messages) {
                if (msg.key == acmacs::messages::key::location_not_found) {
                    if (msg.value == "CRIE")
                        AD_WARNING("CRIE ({}) @@ {}:{}", msg.key, entry.fasta.filename, entry.fasta.line_no);
                    location_not_found.count(msg.value);
                }
                else if (msg.key == acmacs::messages::key::unrecognized_passage)
                    unrecognized_passage.count(msg.value);
                else {
                    AD_WARNING("\"{}\" ({}) -> \"{}\" @@ {}:{}", msg.value, msg.key, *entry.sequence.name(), entry.fasta.filename, entry.fasta.line_no);
                    ++errors;
                }
            }
        }
    }

    const auto report_by_count = [](const acmacs::Counter<std::string>& source, const char* title, bool max_first) {
        fmt::print(stderr, "{}: {}\n", title, source.size());
        if (max_first) {
            for (const auto& entry : source.sorted_max_first())
                fmt::print(stderr, "{:6d} {}\n", entry->second, entry->first);
        }
        else {
            for (const auto& entry : source.counter())
                fmt::print(stderr, "{} {:6d}\n", entry.first, entry.second);
        }
        fmt::print(stderr, "\n");
    };

    if (!unrecognized_passage.empty()) {
        report_by_count(unrecognized_passage, "Unrecognized PASSAGE", true);
        ++errors;
    }

    if (!location_not_found.empty()) {
        report_by_count(location_not_found, "Not found LOCATION", true);
        ++errors;
    }

    fmt::print(stderr, "======================================================================\n\n");
    fmt::print(stderr, "TOTAL: {}\n\n", sequences.size());
    report_by_count(subtypes, "SUBTYPES", true);
    report_by_count(lineages, "LINEAGES", true);
    report_by_count(clades, "CLADES", true);
    report_by_count(labs, "LABS", true);
    report_by_count(isolation_dates, "ISOLATION DATES", false);
    report_by_count(submission_dates, "SUBMISSION DATES", false);

    fmt::print(stderr, "SUBTYPES and sequence lengths (count:seq-length)\n");
    for (const auto& [subtype, entry] : subtypes_sequence_length) {
        if (subtypes[subtype] > 1000) {
            fmt::print(stderr, "  {:<10s}", subtype);
            size_t cnt = 0;
            for (const auto& entry2 : sorted_by_count(entry)) {
                if (entry2.second < 100)
                    break;
                if (cnt > 14) {
                    fmt::print(stderr, "\n            ");
                    cnt = 0;
                }
                fmt::print(stderr, " {:6d}:{:4d}", entry2.second, entry2.first);
                ++cnt;
            }
            fmt::print(stderr, "\n");
        }
    }
    fmt::print(stderr, "\n");

    return errors;

} // report

// ----------------------------------------------------------------------

void report_issues(const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& all_sequences)
{
    using namespace acmacs::seqdb;
    constexpr const auto issue_first = static_cast<size_t>(sequence::issue::not_aligned);

    std::array<acmacs::Counter<std::string>, sequence::number_of_issues> counters;
    for (const auto& sc : all_sequences | ranges::views::filter(scan::fasta::is_translated)) {
        for (auto iss = issue_first; iss < sequence::number_of_issues; ++iss) {
            if (sc.sequence.has_issue(static_cast<sequence::issue>(iss)))
                counters[iss].count(*sc.fasta.type_subtype);
        }
    }
    for (auto iss = issue_first; iss < sequence::number_of_issues; ++iss) {
        if (!counters[iss].empty())
            AD_WARNING("Issue: {}\n{}", sequence::issue_name[iss], counters[iss].report_sorted_max_first());
    }

} // report_issues

// ----------------------------------------------------------------------
