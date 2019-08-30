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
#include "seqdb-3/scan-align.hh"
#include "seqdb-3/scan-deletions.hh"
#include "seqdb-3/scan-clades.hh"
#include "seqdb-3/scan-match-hidb.hh"
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

    option<str>  print_aa_for{*this, "print-aa-for", dflt{""}};
    option<str>  print_not_aligned_for{*this, "print-not-aligned-for", dflt{""}, desc{"ALL or comma separated: H1N,H3,B"}};
    option<str>  print_counter_for{*this, "print-counter-for", dflt{""}};
    option<str>  print_aligned_for{*this, "print-aligned-for", dflt{""}};
    option<bool> print_aa_sizes{*this, "print-aa-sizes"};
    option<bool> print_stat{*this, "stat"};

    option<bool> verbose{*this, 'v', "verbose"};

    argument<str_array> filenames{*this, arg_name{"filename"}, mandatory};
};

static int report(const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& sequences, const Options& opt);
static void report_messages(const acmacs::seqdb::scan::fasta::messages_t& messages);

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        auto [all_sequences, messages] = acmacs::seqdb::scan::fasta::scan(opt.filenames, acmacs::seqdb::scan::fasta::scan_options_t(opt.verbose ? acmacs::debug::yes : acmacs::debug::no));
        fmt::print(stderr, "INFO: Total sequences upon scanning fasta: {:7d}\n", all_sequences.size());

        acmacs::seqdb::scan::fasta::merge_duplicates(all_sequences);
        acmacs::seqdb::scan::fasta::sort_by_date(all_sequences);
        acmacs::seqdb::scan::translate_align(all_sequences);
        acmacs::seqdb::scan::detect_insertions_deletions(all_sequences);
        acmacs::seqdb::scan::detect_lineages_clades(all_sequences);
        acmacs::seqdb::scan::fasta::sort_by_name(all_sequences);
        acmacs::seqdb::scan::match_hidb(all_sequences);
        if (!opt.output_seqdb->empty())
            acmacs::seqdb::create(opt.output_seqdb, all_sequences, opt.whocc_only ? acmacs::seqdb::create_dbs::whocc_only : acmacs::seqdb::create_dbs::all);

        fmt::print(stderr, "INFO: Total sequences upon translating:    {:7d}  aligned: {}\n", all_sequences.size(), ranges::count_if(all_sequences, acmacs::seqdb::scan::fasta::is_aligned));
        fmt::print(stderr, "\n");

        if (!opt.print_counter_for->empty()) {
            const auto chunk = ::string::upper(*opt.print_counter_for);
            const auto found = [&chunk](size_t limit) { return [&chunk, limit](const auto& sc) { return sc.sequence.aa().find(std::string_view(chunk)) < limit; }; };
            for (auto limit : {50, 100, 150, 200, 1000}) {
                acmacs::Counter<std::string> counter;
                for (const auto& sc : all_sequences | ranges::views::filter(acmacs::seqdb::scan::fasta::is_translated) | ranges::views::filter(found(static_cast<size_t>(limit))))
                    counter.count(sc.fasta.type_subtype.h_or_b());
                fmt::print(stderr, "Counter for {} at first {} positions\n{}\n", chunk, limit, counter.report_sorted_max_first());
            }
        }

        if (const auto false_positive = acmacs::seqdb::scan::fasta::report_false_positive(all_sequences, 200); !false_positive.empty())
            fmt::print(stderr, "ERROR: FALSE POSITIVES {}\n{}\n", ranges::count(false_positive, '\n') / 2, false_positive);

        fmt::print(stderr, "{}\n", acmacs::seqdb::scan::fasta::report_dates(all_sequences));
        report_messages(messages);

        if (opt.print_counter_for->empty()) {
            acmacs::Counter<std::string> counter_not_aligned, counter_not_aligned_h;
            for (const auto& sc : all_sequences | ranges::views::filter(acmacs::seqdb::scan::fasta::is_translated) | ranges::views::filter(acmacs::seqdb::scan::fasta::isnot_aligned)) {
                counter_not_aligned.count(*sc.fasta.type_subtype);
                counter_not_aligned_h.count(sc.fasta.type_subtype.h_or_b());
            }
            if (counter_not_aligned_h.total())
                fmt::print(stderr, "WARNING: NOT ALIGNED\n{}\n", counter_not_aligned_h.report_sorted_max_first());
            else
                fmt::print(stderr, "INFO: all aligned\n");
        }

        if (!opt.print_aa_for->empty()) {
            const auto report = acmacs::seqdb::scan::fasta::report_aa(all_sequences, ::string::upper(*opt.print_aa_for), 99999);
            fmt::print("{} {}\n{}\n", *opt.print_aa_for, ranges::count(report, '\n') / 2, report);
        }

        if (!opt.print_not_aligned_for->empty()) {
            const auto report = acmacs::seqdb::scan::fasta::report_not_aligned(all_sequences, ::string::upper(*opt.print_not_aligned_for), 200);
            fmt::print(stderr, "NOT ALIGNED {} {} (name and file reference are at the line end)\n{}\n", *opt.print_not_aligned_for, ranges::count(report, '\n'), report);
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

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------

template <typename Key> static inline acmacs::flat_map_t<Key, size_t> sorted_by_count(const std::map<Key, size_t>& source)
{
    acmacs::flat_map_t<Key, size_t> result(std::begin(source), std::end(source));
    result.sort_by_value_reverse();
    return result;
}

// ----------------------------------------------------------------------

void report_messages(const acmacs::seqdb::scan::fasta::messages_t& messages)
{
    std::map<std::string, acmacs::seqdb::scan::fasta::messages_t, std::less<>> messages_per_key;
    for (const auto& msg : messages)
        messages_per_key.try_emplace(msg.message.key).first->second.push_back(msg);
    for (auto& [key, value] : messages_per_key) {
        // std::sort(std::begin(value), std::end(value));
        if (std::string_view{key} == "location-not-found" || std::string_view{key} == "unrecognized-passage") {
            std::vector<std::string> locations(value.size());
            std::transform(std::begin(value), std::end(value), std::begin(locations), [](const auto& en) { return en.message.value; });
            std::sort(std::begin(locations), std::end(locations));
            const auto end = std::unique(std::begin(locations), std::end(locations));
            fmt::print(stderr, "WARNING: {} ({}):\n", key, end - std::begin(locations));
            fmt::print(stderr, "  \"{}\"\n", ::string::join("\"\n  \"", std::begin(locations), end));
            if (std::string_view{key} == "location-not-found")
                fmt::print(stderr, "locdb \"{}\"\n", ::string::join("\" \"", std::begin(locations), end));
        }
        else {
            fmt::print(stderr, "WARNING: {} ({}):\n", key, value.size());
            for (const auto& val : value)
                fmt::print(stderr, "{}:{}: warning: {} ({})\n", val.filename, val.line_no, val.message.value, key);
        }
        fmt::print(stderr, "\n");
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
        submission_dates.count(entry.sequence.gisaid_last_modified().front().substr(0, 7)); // year-month
        // if (entry.sequence.type_subtype().empty())
        //     fmt::print(stderr, "{}:{}: No subtype for {}\n", entry.fasta.filename, entry.fasta.line_no, *entry.sequence.name);
        if (auto [iter, inserted] = subtypes_sequence_length.emplace(entry.sequence.type_subtype(), std::map<size_t, size_t>{}).first->second.emplace(entry.sequence.nuc().size(), 1UL); !inserted)
            ++iter->second;

        lineages.count_if(!entry.sequence.lineage().empty(), entry.sequence.lineage());
        ranges::for_each(entry.sequence.clades(), [&clades](const auto& clade) { clades.count(clade); });

        if ((opt.all_lab_messages || whocc_lab(entry.sequence)) && (opt.all_subtypes_messages || our_subtype(entry.sequence.type_subtype()))) {
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
