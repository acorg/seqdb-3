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
#include "seqdb-3/fasta.hh"

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

static inline std::string infer_regex(const std::vector<std::string>& sources)
{
    std::vector<std::string> letters(128);
    for (const auto& src : sources) {
        for (size_t pos = 0; pos < src.size(); ++pos) {
            if (ranges::find(letters[pos], src[pos]) == ranges::end(letters[pos]))
                letters[pos].append(1, src[pos]);
        }
    }
    std::string res;
    for (const auto& let : letters) {
        if (let.empty())
            break;
        if (let.size() == 1)
            res.append(let);
        else {
            res.append(1, '[');
            res.append(let);
            res.append(1, ']');
        }
    }
    return res;
}

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<bool> all_lab_messages{*this, "all-lab-messages", desc{"otherwise show messages for WHO CCs only"}};
    option<bool> all_subtypes_messages{*this, "all-subtypes-messages", desc{"otherwise show messages for H1, H3, B only"}};

    option<str>  print_aa_for{*this, "print-aa-for", dflt{""}};
    option<str>  print_not_aligned_for{*this, "print-not-aligned-for", dflt{""}};
    option<str>  print_counter_for{*this, "print-counter-for", dflt{""}};

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

        acmacs::seqdb::fasta::translate_align(all_sequences);
        fmt::print(stderr, "TOTAL sequences upon translating:    {:7d}  aligned: {}\n", all_sequences.size(), ranges::count_if(all_sequences, acmacs::seqdb::fasta::is_aligned));
        fmt::print(stderr, "\n");

        if (!opt.print_counter_for->empty()) {
            acmacs::Counter<std::string> counter;
            const auto chunk = ::string::upper(*opt.print_counter_for);
            for (const auto& sc : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated) | ranges::view::filter([&chunk](const auto& sc) {
                                      const auto pos = sc.sequence.aa().find(std::string_view(chunk));
                                      return pos < 100;
                                  }))
                counter.count(sc.fasta.type_subtype.size() > 4 ? sc.fasta.type_subtype.substr(2, 3) : sc.fasta.type_subtype);
            counter.report_sorted_max_first(fmt::format("Counter for {}\n", chunk), "\n");
        }

        if (const auto false_positive = acmacs::seqdb::fasta::report_false_positive(all_sequences, 200); !false_positive.empty())
            fmt::print(stderr, "FALSE POSITIVES {}\n{}\n", ranges::count(false_positive, '\n') / 2, false_positive);
        // if (const auto not_aligned = acmacs::seqdb::fasta::report_not_aligned(all_sequences, "A(H3N", 200); !not_aligned.empty())
        //     fmt::print(stderr, "H3 NOT ALIGNED {}\n{}\n", ranges::count(not_aligned, '\n') / 2, not_aligned);
        // if (const auto not_aligned = acmacs::seqdb::fasta::report_not_aligned(all_sequences, "A(H4N", 200); !not_aligned.empty())
        //     fmt::print(stderr, "H4 NOT ALIGNED {}\n{}\n", ranges::count(not_aligned, '\n') / 2, not_aligned);

        acmacs::Counter<std::string> counter_not_aligned, counter_not_aligned_h;
        for (const auto& sc : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated) | ranges::view::filter(acmacs::seqdb::fasta::isnot_aligned)) {
            counter_not_aligned.count(sc.fasta.type_subtype); // .substr(2, 3));
            counter_not_aligned_h.count(sc.fasta.type_subtype.size() > 4 ? sc.fasta.type_subtype.substr(2, 3) : sc.fasta.type_subtype);
        }
        counter_not_aligned_h.report_sorted_max_first("NOT ALIGNED\n", "\n");
        counter_not_aligned.report_sorted_max_first("NOT ALIGNED\n", "\n");

        if (!opt.print_aa_for->empty()) {
            const auto report = acmacs::seqdb::fasta::report_aa(all_sequences, ::string::upper(*opt.print_aa_for), 200);
            fmt::print("{} {}\n{}\n", *opt.print_aa_for, ranges::count(report, '\n') / 2, report);
        }

        if (!opt.print_not_aligned_for->empty()) {
            const auto report = acmacs::seqdb::fasta::report_not_aligned(all_sequences, ::string::upper(*opt.print_not_aligned_for), 200);
            fmt::print(stderr, "NOT ALIGNED {} {}\n{}\n", *opt.print_not_aligned_for, ranges::count(report, '\n') / 2, report);
        }

        // std::vector<std::string> qk;
        // for (const auto& sc : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::isnot_aligned)) {
        //     if (const auto pos = sc.sequence.aa().find("QK"); pos < 50 && sc.sequence.aa().substr(pos + 16, 2) == "HH")
        //         qk.emplace_back(sc.sequence.aa().substr(pos, 32));
        // }
        // fmt::print(stderr, "QK: {}\n", infer_regex(qk));

        // ranges::sort(qk);
        // fmt::print(stderr, "QK {}\n{}\n", qk.size(), string::join("\n", qk));

        // acmacs::Counter<std::string> mktii; // type_subtype counter for MKTII
        // acmacs::Counter<std::string> qkip;
        // std::vector<std::pair<std::string,std::string>> mktii_not_h3;
        // for (const auto& seq : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated)) {
        //     const auto aa = seq.sequence.aa();
        //     if (const auto pos = aa.find("MKTII"); pos < 50 && (aa[pos + 16] == 'Q' || aa[pos + 15] == 'A')) { // aa.substr(pos + 15, 2) != "DR") { // DR[ISV]C - start of the B sequence (signal
        //     peptide is 15 aas!)
        //         mktii.count(seq.fasta.type_subtype);
        //         qkip.count(aa.substr(pos + 15, 2));
        //         if (seq.fasta.type_subtype.substr(0, 4) != "A(H3" && seq.fasta.type_subtype.substr(0, 4) != "A(H0")
        //             mktii_not_h3.emplace_back(fmt::format("{} {}", seq.fasta.type_subtype, seq.fasta.entry_name), aa.substr(0, 200));
        //     }
        // }
        // mktii.report_sorted_max_first("MKTII\n", "\n");
        // qkip.report_sorted_max_first("QKIP\n", "\n");

        // for (const auto& e2 : mktii_not_h3)
        //     fmt::print(stderr, "{}\n{}\n", e2.first, e2.second);
        // fmt::print(stderr, "\n");

        // acmacs::Counter<std::string> mkt; // type_subtype counter for MKT
        // fmt::print(stderr, "MKT\n");
        // for (const auto& seq : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated)) {
        //     if (const auto pos = seq.sequence.aa().find("MKT"); pos < 50)
        //         mkt.count(seq.fasta.type_subtype);
        // }
        // mkt.report_sorted_max_first();
        // fmt::print(stderr, "\n");

        // acmacs::Counter<std::string> mkt;
        // for (const auto& seq : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::is_translated)) {
        //     if (const auto pos = seq.sequence.aa().find("MKT"); pos < 100)
        //         mkt.count(seq.sequence.aa().substr(pos, 16));
        //     else if (const auto pos2 = seq.sequence.aa().find("MK"); pos2 < 100)
        //         mkt.count(seq.sequence.aa().substr(pos2, 16));
        // }
        // mkt.report_sorted_max_first();

        // for (const auto& seq : all_sequences | ranges::view::filter(acmacs::seqdb::fasta::isnot_aligned)) {
        //     fmt::print(stderr, "{}\n", seq.sequence.aa().substr(0, 200));
        // }

        const auto errors = 0; // report(all_sequences, opt);

        return errors;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

//        // ----------------------------------------------------------------------
//
//        // std::set<char> all_aa;
//        // for (const auto& seq_e : all_sequences) {
//        //     for (char aa : seq_e.seq.sequence.aa())
//        //         all_aa.insert(aa);
//        // }
//        // std::string all_aa_s(std::begin(all_aa), std::end(all_aa));
//        // ranges::sort(all_aa_s);
//        // fmt::print(stderr, "ALL AA: [{}]\n", all_aa_s);
//
//            // ----------------------------------------------------------------------
//
//        const auto if_aligned = [](const auto& entry) -> bool { return entry.aligned; };
//        const auto if_not_aligned = [](const auto& entry) -> bool { return !entry.aligned; };
//
//        std::vector<std::reference_wrapper<acmacs::seqdb::v3::fasta::scan_result_t>> aligned, not_aligned;
//        ranges::copy(ranges::view::filter(all_sequences, [](const auto& entry) { return entry.aligned; }), ranges::back_inserter(aligned));
//        ranges::copy(ranges::view::filter(all_sequences, [](const auto& entry) { return !entry.aligned; }), ranges::back_inserter(not_aligned));
//        fmt::print(stderr, "ALIGNED: {}  not aligned: {}\n", aligned.size(), not_aligned.size());
//
//        // ----------------------------------------------------------------------
//
//        // const auto update = [](auto& counter, auto&& source) {
//        //     if (auto [iter, inserted] = counter.emplace(source, 1UL); !inserted)
//        //         ++iter->second;
//        // };
//
//        std::vector<acmacs::CounterChar> occurences_per_pos(550);
//        ranges::for_each(aligned, [&occurences_per_pos](const auto& entry) {
//            const auto seq = entry.get().seq.sequence.aa_aligned();
//            for (size_t pos = 0; pos < std::min(seq.size(), occurences_per_pos.size()); ++pos) {
//                if (seq[pos] != 'X' && seq[pos] != '-')
//                    occurences_per_pos[pos].count(seq[pos]);
//            }
//        });
//        const auto occurences_threshold = occurences_per_pos[0].max().second;
//        fmt::print(stderr, "occurences_threshold {} {}\n", occurences_threshold, double(occurences_threshold) / aligned.size());
//
//        std::string common(occurences_per_pos.size(), '.');
//        for (size_t pos = 0; pos < occurences_per_pos.size(); ++pos) {
//            const auto [aa, count] = occurences_per_pos[pos].max();
//            if (count >= occurences_threshold)
//                common[pos] = aa;
//        }
//        fmt::print(stderr, "{}\n", common);
//
//        acmacs::Counter<size_t> hamming_distance_counter;
//        for (auto& entry : not_aligned) {
//            const auto aa = entry.get().seq.sequence.aa();
//            std::vector<std::pair<size_t, size_t>> pos_hamdist;
//            for (auto pos = aa.find(common[0]); pos != std::string::npos; pos = aa.find(common[0], pos + 1)) {
//                pos_hamdist.emplace_back(pos, acmacs::seqdb::hamming_distance_not_considering(aa.substr(pos), common, '.'));
//            }
//            if (const auto me = std::min_element(std::begin(pos_hamdist), std::end(pos_hamdist), [](const auto& e1, const auto& e2) { return e1.second < e2.second; }); me != std::end(pos_hamdist)) {
//                hamming_distance_counter.count(me->second);
//                if (me->second < 30) {
//                    entry.get().aligned = true;
//                    // fmt::print(stderr, "ham dist {} {}   {}\n{}\n", me->second, me->first, entry.get().seq.fasta_name, aa);
//                }
//                else if (me->second < 200) {
//                    fmt::print(stderr, "ham dist {} {}   {}\n{}\n", me->second, me->first, entry.get().seq.fasta_name, aa.substr(0, 150));
//                }
//            }
//        }
//
//        std::vector<std::reference_wrapper<acmacs::seqdb::v3::fasta::scan_result_t>> aligned_2;
//        ranges::copy(ranges::view::filter(all_sequences, [](const auto& entry) { return entry.aligned; }), ranges::back_inserter(aligned_2));
//        fmt::print(stderr, "ALIGNED2: {}\n", aligned_2.size());
//
//        fmt::print(stderr, "hamming_distance_counter {}\n", hamming_distance_counter.counter());
//
//        // std::vector<std::pair<std::string_view, size_t>> chunk_offset;
//        // const auto split_data = acmacs::string::split(common, ".", acmacs::string::Split::RemoveEmpty);
//        // for (const auto& chunk : split_data) {
//        //     if (chunk.size() > 5) {
//        //         chunk_offset.emplace_back(chunk, static_cast<size_t>(chunk.data() - common.data()));
//        //         // fmt::print(stderr, "{:3d} {}\n", chunk.data() - common.data(), chunk);
//        //     }
//        // }
//
//        // for (auto& seq_e : all_sequences) {
//        //     if (! if_aligned(seq_e)) {
//        //         std::vector<std::string_view::size_type> offsets(chunk_offset.size(), std::string_view::npos);
//        //         size_t good_offsets = 0;
//        //         for (size_t no = 0; no < chunk_offset.size(); ++no) {
//        //             offsets[no] = seq_e.seq.sequence.aa().find(chunk_offset[no].first);
//        //             if (offsets[no] != std::string_view::npos && (no == 0 || offsets[no] > offsets[no-1]))
//        //                 ++good_offsets;
//        //         }
//        //         if (offsets[0] != std::string_view::npos && good_offsets > (chunk_offset.size() * 0.8)) {
//        //             seq_e.seq.sequence.set_shift_aa(static_cast<int>(offsets[0]) - static_cast<int>(chunk_offset[0].second));
//        //             seq_e.aligned = true;
//        //         }
//        //     }
//        // }
//
//        fmt::print(stderr, "ALIGNED: {}\n", ranges::count_if(all_sequences, if_aligned));
//        fmt::print(stderr, "H3: {}\n", ranges::count_if(all_sequences, [](const auto& entry) -> bool { return entry.seq.type_subtype == "A(H3N2)"; }));
//
//        fmt::print(stderr, "\nAligned not H3: {}\n", ranges::count_if(all_sequences, [](const auto& entry) -> bool { return entry.aligned && entry.seq.type_subtype.substr(0, 4) != "A(H3" && entry.seq.type_subtype != "A(H0N0)"; }));
//        for (auto& seq_e : all_sequences) {
//            if (if_aligned(seq_e) && seq_e.seq.type_subtype.substr(0, 4) != "A(H3" && seq_e.seq.type_subtype != "A(H0N0)")
//                fmt::print(stderr, "{}\n{}\n", seq_e.seq.fasta_name, seq_e.seq.sequence.aa_aligned());
//        }
//
//        fmt::print(stderr, "\nNot Aligned H3: {}\n", ranges::count_if(all_sequences, [](const auto& entry) -> bool { return !entry.aligned && !entry.seq.sequence.aa().empty() && entry.seq.type_subtype.substr(0, 4) == "A(H3"; }));
//        for (auto& seq_e : all_sequences) {
//            if (!seq_e.aligned && seq_e.seq.type_subtype.substr(0, 4) == "A(H3")
//                fmt::print(stderr, "{}\n{}\n", seq_e.seq.fasta_name, seq_e.seq.sequence.aa().substr(0, 200));
//        }
//
//
//        // ----------------------------------------------------------------------
//
//        // std::vector<std::string_view> aligned;
//        // const auto aligned = ranges::view::all(all_sequences) | ranges::view::filter(if_aligned) | ranges::view::transform([](const auto& entry) { return entry.seq.sequence.aa_aligned(); });
//        // // | ranges::view::for_each([](std::string entry) { return 7.0; });
//        // // aligned | ranges::view::for_each([](std::string entry) { return 7.0; });
//        // for (const auto& seq : aligned) {
//        // }
//        // }
//
//        // FILE* h3_seq_fd = std::fopen("/d/h3.fas", "w");
//        // ranges::for_each(ranges::view::all(all_sequences) | ranges::view::filter(if_aligned), [&h3_seq_fd](const auto& entry) { fmt::print(h3_seq_fd, "{}\n{}\n", *entry.seq.name,
//        // entry.seq.sequence.aa_aligned()); }); std::fclose(h3_seq_fd);
//
//        // FILE* h3_seq_fd = std::fopen("/d/h3.fas", "w");
//
//        // size_t aligned = 0, potential = 0;
//        // for (const auto& seq_e : all_sequences) {
//        //     if (seq_e.seq.type_subtype == "A(H3N2)") {
//        //         if (seq_e.aligned) {
//        //             ++aligned;
//        //             if (const auto yr = acmacs::virus::year(seq_e.seq.name); yr.has_value() && *yr > 2017)
//        //                 fmt::print(h3_seq_fd, "{}\n{}\n", *seq_e.seq.name, seq_e.seq.sequence.aa_aligned());
//        //         }
//        //         else if (seq_e.seq.host->empty() && !seq_e.seq.sequence.aa().empty()) {
//        //             ++potential;
//        //             if (whocc_lab(seq_e.seq.lab))
//        //                 fmt::print(stderr, "!! {} NOT H3? {}\n{}\n", seq_e.seq.lab, seq_e.seq.fasta_name, std::string_view(seq_e.seq.sequence.aa().data(), 200));
//        //             else
//        //                 fmt::print(stderr, "NOT H3? {}\n{}\n", seq_e.seq.fasta_name, std::string_view(seq_e.seq.sequence.aa().data(), 200));
//        //         }
//        //     }
//        // }
//        // std::fclose(h3_seq_fd);
//
//        // fmt::print(stderr, "ALIGNED: {}  Potential: {}\n", aligned, potential);


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
