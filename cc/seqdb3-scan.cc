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

        // ----------------------------------------------------------------------

        // std::set<char> all_aa;
        // for (const auto& seq_e : all_sequences) {
        //     for (char aa : seq_e.seq.sequence.aa())
        //         all_aa.insert(aa);
        // }
        // std::string all_aa_s(std::begin(all_aa), std::end(all_aa));
        // ranges::sort(all_aa_s);
        // fmt::print(stderr, "ALL AA: [{}]\n", all_aa_s);

            // ----------------------------------------------------------------------

        const auto if_aligned = [](const auto& entry) -> bool { return entry.aligned; };

        std::vector<std::reference_wrapper<acmacs::seqdb::v3::fasta::scan_result_t>> aligned, not_aligned;
        ranges::copy(ranges::view::filter(all_sequences, [](const auto& entry) { return entry.aligned; }), ranges::back_inserter(aligned));
        ranges::copy(ranges::view::filter(all_sequences, [](const auto& entry) { return !entry.aligned; }), ranges::back_inserter(not_aligned));
        fmt::print(stderr, "ALIGNED: {}  not aligned: {}\n", aligned.size(), not_aligned.size());

        // ----------------------------------------------------------------------

        // const auto update = [](auto& counter, auto&& source) {
        //     if (auto [iter, inserted] = counter.emplace(source, 1UL); !inserted)
        //         ++iter->second;
        // };

        std::vector<acmacs::CounterChar> occurences_per_pos(550);
        ranges::for_each(aligned, [&occurences_per_pos](const auto& entry) {
            const auto seq = entry.get().seq.sequence.aa_aligned();
            for (size_t pos = 0; pos < std::min(seq.size(), occurences_per_pos.size()); ++pos) {
                if (seq[pos] != 'X' && seq[pos] != '-')
                    occurences_per_pos[pos].count(seq[pos]);
            }
        });
        const auto occurences_threshold = occurences_per_pos[0].max().second;
        fmt::print(stderr, "occurences_threshold {} {}\n", occurences_threshold, double(occurences_threshold) / aligned.size());

        std::string common(occurences_per_pos.size(), '.');
        for (size_t pos = 0; pos < occurences_per_pos.size(); ++pos) {
            const auto [aa, count] = occurences_per_pos[pos].max();
            if (count >= occurences_threshold)
                common[pos] = aa;
        }
        fmt::print(stderr, "{}\n", common);

        acmacs::Counter<size_t> hamming_distance_counter;
        for (auto& entry : not_aligned) {
            const auto aa = entry.get().seq.sequence.aa();
            std::vector<std::pair<size_t, size_t>> pos_hamdist;
            for (auto pos = aa.find(common[0]); pos != std::string::npos; pos = aa.find(common[0], pos + 1)) {
                pos_hamdist.emplace_back(pos, acmacs::seqdb::hamming_distance_not_considering(aa.substr(pos), common, '.'));
            }
            if (const auto me = std::min_element(std::begin(pos_hamdist), std::end(pos_hamdist), [](const auto& e1, const auto& e2) { return e1.second < e2.second; }); me != std::end(pos_hamdist)) {
                hamming_distance_counter.count(me->second);
                if (me->second < 30) {
                    entry.get().aligned = true;
                    // fmt::print(stderr, "ham dist {} {}   {}\n{}\n", me->second, me->first, entry.get().seq.fasta_name, aa);
                }
                else if (me->second < 200) {
                    fmt::print(stderr, "ham dist {} {}   {}\n{}\n", me->second, me->first, entry.get().seq.fasta_name, aa.substr(0, 150));
                }
            }
        }

        std::vector<std::reference_wrapper<acmacs::seqdb::v3::fasta::scan_result_t>> aligned_2;
        ranges::copy(ranges::view::filter(all_sequences, [](const auto& entry) { return entry.aligned; }), ranges::back_inserter(aligned_2));
        fmt::print(stderr, "ALIGNED2: {}\n", aligned_2.size());

        fmt::print(stderr, "hamming_distance_counter {}\n", hamming_distance_counter.counter());

        // std::vector<std::pair<std::string_view, size_t>> chunk_offset;
        // const auto split_data = acmacs::string::split(common, ".", acmacs::string::Split::RemoveEmpty);
        // for (const auto& chunk : split_data) {
        //     if (chunk.size() > 5) {
        //         chunk_offset.emplace_back(chunk, static_cast<size_t>(chunk.data() - common.data()));
        //         // fmt::print(stderr, "{:3d} {}\n", chunk.data() - common.data(), chunk);
        //     }
        // }

        // for (auto& seq_e : all_sequences) {
        //     if (! if_aligned(seq_e)) {
        //         std::vector<std::string_view::size_type> offsets(chunk_offset.size(), std::string_view::npos);
        //         size_t good_offsets = 0;
        //         for (size_t no = 0; no < chunk_offset.size(); ++no) {
        //             offsets[no] = seq_e.seq.sequence.aa().find(chunk_offset[no].first);
        //             if (offsets[no] != std::string_view::npos && (no == 0 || offsets[no] > offsets[no-1]))
        //                 ++good_offsets;
        //         }
        //         if (offsets[0] != std::string_view::npos && good_offsets > (chunk_offset.size() * 0.8)) {
        //             seq_e.seq.sequence.set_shift_aa(static_cast<int>(offsets[0]) - static_cast<int>(chunk_offset[0].second));
        //             seq_e.aligned = true;
        //         }
        //     }
        // }

        fmt::print(stderr, "ALIGNED: {}\n", ranges::count_if(all_sequences, if_aligned));
        fmt::print(stderr, "H3: {}\n", ranges::count_if(all_sequences, [](const auto& entry) -> bool { return entry.seq.type_subtype == "A(H3N2)"; }));

        fmt::print(stderr, "\nAligned not H3: {}\n", ranges::count_if(all_sequences, [](const auto& entry) -> bool { return entry.aligned && entry.seq.type_subtype.substr(0, 4) != "A(H3" && entry.seq.type_subtype != "A(H0N0)"; }));
        for (auto& seq_e : all_sequences) {
            if (if_aligned(seq_e) && seq_e.seq.type_subtype.substr(0, 4) != "A(H3" && seq_e.seq.type_subtype != "A(H0N0)")
                fmt::print(stderr, "{}\n{}\n", seq_e.seq.fasta_name, seq_e.seq.sequence.aa_aligned());
        }

        fmt::print(stderr, "\nNot Aligned H3: {}\n", ranges::count_if(all_sequences, [](const auto& entry) -> bool { return !entry.aligned && !entry.seq.sequence.aa().empty() && entry.seq.type_subtype.substr(0, 4) == "A(H3"; }));
        for (auto& seq_e : all_sequences) {
            if (!seq_e.aligned && !seq_e.seq.sequence.aa().empty() && seq_e.seq.type_subtype.substr(0, 4) == "A(H3")
                fmt::print(stderr, "{}\n{}\n", seq_e.seq.fasta_name, seq_e.seq.sequence.aa().substr(0, 200));
        }


        // ----------------------------------------------------------------------

        // std::vector<std::string_view> aligned;
        // const auto aligned = ranges::view::all(all_sequences) | ranges::view::filter(if_aligned) | ranges::view::transform([](const auto& entry) { return entry.seq.sequence.aa_aligned(); });
        // // | ranges::view::for_each([](std::string entry) { return 7.0; });
        // // aligned | ranges::view::for_each([](std::string entry) { return 7.0; });
        // for (const auto& seq : aligned) {
        // }
        // }

        // FILE* h3_seq_fd = std::fopen("/d/h3.fas", "w");
        // ranges::for_each(ranges::view::all(all_sequences) | ranges::view::filter(if_aligned), [&h3_seq_fd](const auto& entry) { fmt::print(h3_seq_fd, "{}\n{}\n", *entry.seq.name,
        // entry.seq.sequence.aa_aligned()); }); std::fclose(h3_seq_fd);

        // FILE* h3_seq_fd = std::fopen("/d/h3.fas", "w");

        // size_t aligned = 0, potential = 0;
        // for (const auto& seq_e : all_sequences) {
        //     if (seq_e.seq.type_subtype == "A(H3N2)") {
        //         if (seq_e.aligned) {
        //             ++aligned;
        //             if (const auto yr = acmacs::virus::year(seq_e.seq.name); yr.has_value() && *yr > 2017)
        //                 fmt::print(h3_seq_fd, "{}\n{}\n", *seq_e.seq.name, seq_e.seq.sequence.aa_aligned());
        //         }
        //         else if (seq_e.seq.host->empty() && !seq_e.seq.sequence.aa().empty()) {
        //             ++potential;
        //             if (whocc_lab(seq_e.seq.lab))
        //                 fmt::print(stderr, "!! {} NOT H3? {}\n{}\n", seq_e.seq.lab, seq_e.seq.fasta_name, std::string_view(seq_e.seq.sequence.aa().data(), 200));
        //             else
        //                 fmt::print(stderr, "NOT H3? {}\n{}\n", seq_e.seq.fasta_name, std::string_view(seq_e.seq.sequence.aa().data(), 200));
        //         }
        //     }
        // }
        // std::fclose(h3_seq_fd);

        // fmt::print(stderr, "ALIGNED: {}  Potential: {}\n", aligned, potential);

        const auto errors = report(all_sequences, opt);

        return errors;
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
