#include "acmacs-base/counter.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/omp.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/hamming-distance.hh"
#include "seqdb-3/log.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_mean(size_t threshold, size_t size_threshold)
{
    if (threshold > 0 && size_threshold > 0 && !refs_.empty()) {
        struct Entry
        {
            sequence_aligned_ref_t nucs;
            size_t hamming_distance_sum{0};
            size_t ref_index{static_cast<size_t>(-1)};
            std::string_view date;

            void assign(size_t index, const ref& ref)
            {
                nucs = ref.seq().nuc_aligned_master();
                date = ref.entry->date();
                hamming_distance_sum = 0;
                ref_index = index;
            }
        };

        std::vector<Entry> entries(refs_.size());
        for (size_t index{0}; index < refs_.size(); ++index)
            entries[index].assign(index, refs_[index]);
        entries.erase(std::remove_if(std::begin(entries), std::end(entries), [](const auto& en) { return en.nucs.empty(); }), entries.end());
        std::sort(entries.begin(), entries.end(), [](const auto& e1, const auto& e2) { return e1.date > e2.date; }); // most recent first
        if (entries.size() > size_threshold)                                                                         // keep few most recent before comparing
            entries.erase(std::next(entries.begin(), static_cast<ssize_t>(size_threshold)), entries.end());

        size_t total{0};
        for (size_t i1{0}; i1 < entries.size(); ++i1) {
            for (size_t i2{i1 + 1}; i2 < entries.size(); ++i2) {
                const auto hd{hamming_distance(entries[i1].nucs, entries[i2].nucs, hamming_distance_by_shortest::no)};
                entries[i1].hamming_distance_sum += hd;
                entries[i2].hamming_distance_sum += hd;
                ++total;
            }
        }

        const auto base_ref_index = std::min_element(entries.begin(), entries.end(), [](const auto& e1, const auto& e2) { return e1.hamming_distance_sum < e2.hamming_distance_sum; })->ref_index;
        const auto base_seq_id = refs_[base_ref_index].seq_id();
        // AD_INFO("base sequence to exclude by hamming distance {} (--nuc-hamming-distance-mean-threshold {})", base_seq_id, threshold);

        return nuc_hamming_distance_to(threshold, base_seq_id);
    }
    else
        return *this;

} // acmacs::seqdb::v3::subset::nuc_hamming_distance_mean

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_to(size_t threshold, std::string_view seq_id)
{
    if (!seq_id.empty()) {
        const auto& seqdb = acmacs::seqdb::get();
        const auto compare_to = seqdb.select_by_seq_id(seq_id);
        if (compare_to.empty())
            throw std::runtime_error{fmt::format("no sequences with seq-id \"{}\" found (seqdb::v3::subset::nuc_hamming_distance_to)", seq_id)};
        const auto before{refs_.size()};
        refs_.erase(std::remove_if(std::next(std::begin(refs_)), std::end(refs_),
                                   [threshold, &seqdb, comapre_to_seq = compare_to.front().nuc_aligned(seqdb)](auto& en) {
                                       en.hamming_distance = hamming_distance(en.nuc_aligned(seqdb), comapre_to_seq, hamming_distance_by_shortest::no);
                                       return en.hamming_distance >= threshold;
                                   }),
                    std::end(refs_));
        const auto after{refs_.size()};
        if (before - after) {
            AD_INFO("{} sequences removed ({} left) which are too far from {}, threshold: {}", before - after, after, seq_id, threshold);
            if ((before - after) > (before / 4))
                AD_WARNING("too many sequences removed ({} or {:.1f}%) that are too far from {}, hamming distance threshold: {}", before - after,
                           static_cast<double>(before - after) / static_cast<double>(before) * 100.0, seq_id, threshold);
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_hamming_distance_to

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base(size_t threshold, bool do_filter)
{
    if (do_filter) {
        const auto& seqdb = acmacs::seqdb::get();
        const auto before{refs_.size()};
        refs_.erase(std::remove_if(std::next(std::begin(refs_)), std::end(refs_),
                                   [threshold, &seqdb, base_seq = refs_.front().nuc_aligned(seqdb)](auto& en) {
                                       en.hamming_distance = hamming_distance(en.nuc_aligned(seqdb), base_seq, hamming_distance_by_shortest::no);
                                       return en.hamming_distance >= threshold;
                                   }),
                    std::end(refs_));
        const auto after{refs_.size()};
        AD_LOG(acmacs::log::sequences, "{} sequences removed ({} left) which are too far from the base seq, threshold: {}", before - after, after, threshold);
        if ((before - after) > (before / 4))
            AD_WARNING("too many sequences removed ({} or {:.1f}%) that are too far from the base sequence, hamming distance threshold: {}", before - after, static_cast<double>(before - after) / static_cast<double>(before) * 100.0, threshold);
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_hamming_distance(bool do_report)
{
    if (do_report) {
        std::vector<const ref*> refs(refs_.size());
        std::transform(std::begin(refs_), std::end(refs_), std::begin(refs), [](const auto& rr) { return &rr; });
        std::sort(std::begin(refs), std::end(refs), [](const ref* r1, const ref* r2) { return r1->hamming_distance > r2->hamming_distance; });
        for (const auto& en : refs)
            fmt::print("{:4d}  {}\n", en->hamming_distance, en->seq_id());
    }
    return *this;

} // acmacs::seqdb::v3::subset::report_hamming_distance

// ----------------------------------------------------------------------

// Eu's algortihm of subsseting 2019-07-23

// 1. Find first group master sequence. I think good starting sequence
// is the most recent one that matched against hidb. Algorithm also
// prefers matched sequences to make more antigens marked in the sig
// pages.
//
// 2. Compute hamming distance between rest sequences and the master
// sequence, sort rest sequences by hamming distance, smaller first.
//
// 3. Find group end, i.e. first sequence that has hamming distance to
// the group master bigger than dist_threshold. Assign group no to
// this group. Sort group (keep group master first) by number of hi
// names (most number of names first) and by date (most recent first).
//
// 4. Next group master is the first sequence after group end. Repeat
// 2-3-4 until all sequences are processed.
//
// 5. Select masters (first sequences) of every group. If there are
// too many groups, more than output_size, then just used first
// output_size groups. If output_size > number of groups, select the
// second sequence in each group (if group size > 1). Do it until
// output_size sequences selected.

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::group_by_hamming_distance(const Seqdb& seqdb, size_t dist_threshold, size_t output_size)
{
    if (dist_threshold > 0) {
        const auto compute_hamming_distance = [&seqdb](std::string_view master_aa, auto first, auto last) {
            std::for_each(first, last, [&seqdb,master_aa](auto& ref) { ref.hamming_distance = hamming_distance(master_aa, ref.aa_aligned(seqdb)); });
        };

        const auto sort_by_hamming_distance = [](auto first, auto last) { std::sort(first, last, [](const auto& e1, const auto& e2) { return e1.hamming_distance < e2.hamming_distance; }); };

        const auto find_group_end = [dist_threshold](auto first, auto last) { return std::find_if(first, last, [dist_threshold](const auto& en) { return en.hamming_distance >= dist_threshold; }); };

        const auto assign_group_no = [](auto first, auto last, size_t group_no) { std::for_each(first, last, [group_no](auto& en) { en.group_no = group_no; }); };

        const auto sort_by_hi_names = [](auto first, auto last) {
            std::sort(first, last, [](const auto& e1, const auto& e2) {
                return e1.seq().hi_names.size() == e2.seq().hi_names.size() ? e1.entry->date() > e2.entry->date() : e1.seq().hi_names.size() > e2.seq().hi_names.size();
            });
        };

        // ----------------------------------------------------------------------

        std::iter_swap(std::begin(refs_), most_recent_with_hi_name());
        auto group_first = std::begin(refs_);
        acmacs::Counter<ssize_t> counter_group_size;
        for (size_t group_no = 1; group_first != std::end(refs_); ++group_no) {
            const auto group_master_aa_aligned = group_first->aa_aligned(seqdb);
            const auto group_second = std::next(group_first);
            // fmt::print("DEBUG: group {} master: {} {} rest size: {}\n", group_no, group_first->seq_id(), group_first->entry->date(), std::end(refs_) - group_first);
            compute_hamming_distance(group_master_aa_aligned, group_second, std::end(refs_));
            sort_by_hamming_distance(group_second, std::end(refs_));
            const auto group_last = find_group_end(group_second, std::end(refs_));
            assign_group_no(group_first, group_last, group_no);
            sort_by_hi_names(group_no == 1 ? group_second : group_first, group_last);
            counter_group_size.count(group_last - group_first);
            group_first = group_last;
        }
        // fmt::print(stderr, "DEBUG: (num-groups:group-size): {}\n", counter_group_size.report_sorted_max_first(" {second}:{first}"));
        // fmt::print(stderr, "DEBUG: total groups: {}\n", refs_.back().group_no);
        if (refs_.back().group_no > output_size) {
            // too many groups, take one seq from each group starting with group 1, ignore groups with high numbers (furtherst from the recent strain)
            ref_indexes to_remove;
            size_t prev_group = 0;
            for (auto [index, ref] : acmacs::enumerate(refs_)) {
                if (ref.group_no == prev_group)
                    to_remove.push_back(index);
                else {
                    prev_group = ref.group_no;
                    if (prev_group > output_size)
                        to_remove.push_back(index);
                }
            }
            remove(to_remove);
        }
        else {
            // too few groups
            ref_indexes to_keep_indexes;
            size_t to_keep = 0;
            size_t prev_to_keep = output_size;
            while (to_keep < output_size && prev_to_keep != to_keep) {
                prev_to_keep = to_keep;
                size_t group_no = 1;
                for (auto [index, ref] : acmacs::enumerate(refs_)) {
                    if (ref.group_no >= group_no) {
                        to_keep_indexes.push_back(index);
                        ++to_keep;
                        group_no = ref.group_no + 1;
                    }
                    if (to_keep >= output_size)
                        break;
                }
                // fmt::print(stderr, "DEBUG: to_keep {} group_no {}\n", to_keep, group_no);
            }
            keep(to_keep_indexes);
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::group_by_hamming_distance

// ----------------------------------------------------------------------

// davipatti algorithm 2019-07-23 9:58
// > 1. pick a random strain, put in selection
// > 2. pick random strain. if it has a distance < d to anything in in selection then discard it. else, add it to selection.
// > 3. repeat 3 until you have as many strains, n, as you want, or until no more strains to pick
//
// Problems: need to prioritize picking hidb matched sequences.
//
// > parameter d would have to be tuned if d=0, this is just randomly
// > sampling strains if d is very high, only very dissimilar strains will
// > make it into selection, and selection would be small ideally d would
// > be as high as possible such that the number of strains in the
// > selection is close to n
//
// Looks like we need to use a search for d, i.e. we do not stop on
// finding n strains at the step 4 and have to find all to learn how many
// redundant strains there are. And then pick d producing number of
// strains closer to n (I guess having slightly more than n is better
// than having slightly less) and cut it, if necessary.
//
// > i foresee this algorithm being run initially to make a selection when
// > new sequences come in, repeat step 3 above, but just on new strains
// > so, original members stay in selection anything novel enough gets
// > added to the selection selection slowly grows over time
//
// No. The size of selection must be the same (as close to 4k as possible).

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::subset_by_hamming_distance_random(const Seqdb& seqdb, bool do_subset, size_t output_size)
{
    if (do_subset && !refs_.empty()) {
        std::mt19937 generator{std::random_device()()};
        const auto random_from = [&generator](auto first, auto last) {
            std::uniform_int_distribution<ssize_t> distribution(0, last - first - 1);
            return std::next(first, distribution(generator));
        };

        const auto minimal_distance_less_than = [&seqdb](auto first, auto last, std::string_view picked_aa, size_t distance_threshold) -> bool {
            return std::any_of(first, last, [&seqdb, picked_aa, distance_threshold](const auto& en) { return hamming_distance(picked_aa, en.aa_aligned(seqdb)) < distance_threshold; });
        };

        decltype(refs_) best_data;
        for (size_t distance_threshold = 1; distance_threshold < 10; ++distance_threshold) {
            auto data = refs_;
            std::iter_swap(std::begin(data), random_from(std::begin(data), std::end(data)));
            auto selection_start = std::begin(data), selection_end = std::next(selection_start), discarded_start = std::end(data);
            while (discarded_start > selection_end) {
                auto picked = random_from(selection_end, discarded_start);
                if (minimal_distance_less_than(selection_start, selection_end, picked->aa_aligned(seqdb), distance_threshold)) { // discard
                    --discarded_start;
                    std::iter_swap(discarded_start, picked);
                }
                else { // put into selection
                    std::iter_swap(selection_end, picked);
                    ++selection_end;
                }
            }
            fmt::print(stderr, "DEBUG: threshold: {} selection: {}\n", distance_threshold, selection_end - selection_start);
            if (static_cast<size_t>(selection_end - selection_start) < output_size)
                break;          // use previous (best_data)
            best_data.resize(static_cast<size_t>(selection_end - selection_start));
            std::copy(selection_start, selection_end, std::begin(best_data));
        }
        if (best_data.empty())
            throw std::runtime_error(fmt::format("subset_by_hamming_distance_random: best_data is empty"));
        const auto num_seqs = std::min(output_size, best_data.size());
        refs_.resize(num_seqs);
        std::copy(std::begin(best_data), std::next(std::begin(best_data), static_cast<ssize_t>(num_seqs)), std::begin(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::subset_by_hamming_distance_random

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_hamming_bins(const Seqdb& seqdb, size_t bin_size)
{
    if (bin_size > 0) {
        std::vector<std::tuple<std::string, size_t, std::vector<size_t>>> seqids_bins(refs_.size()); // seq_id, max_bin, bins

        auto others = seqdb.all();
        others.subtype(refs_[0].entry->virus_type).host(refs_[0].entry->host()).remove_nuc_duplicates(true, false);

#ifdef _OPENMP
        const int num_threads = omp_get_max_threads();
        // const int slot_size = number_of_antigens() < 1000 ? 4 : 1;
#endif
#pragma omp parallel for default(shared) num_threads(num_threads) firstprivate(others) schedule(static, 1000)
        for (size_t ref_no = 0; ref_no < refs_.size(); ++ref_no) {
            const auto& ref = refs_[ref_no];
            std::get<std::string>(seqids_bins[ref_no]) = ref.seq_id();

            // keep non-zero distances only
            size_t max_distance = 0;
            others.refs_.erase(std::remove_if(std::next(std::begin(others)), std::end(others),
                                              [&seqdb, &max_distance, base_seq = ref.nuc_aligned(seqdb)](auto& en) {
                                                  en.hamming_distance = hamming_distance(en.nuc_aligned(seqdb), base_seq, hamming_distance_by_shortest::yes);
                                                  max_distance = std::max(max_distance, en.hamming_distance);
                                                  return en.hamming_distance == 0;
                                              }),
                               std::end(others));

            const size_t number_of_bins = max_distance / bin_size + 1;
            auto& bins = std::get<2>(seqids_bins[ref_no]);
            bins.resize(number_of_bins, 0ul);
            for (const auto& another : others)
                ++bins[another.hamming_distance / bin_size];
            std::get<1>(seqids_bins[ref_no]) = static_cast<size_t>(std::max_element(std::begin(bins), std::end(bins)) - std::begin(bins));
            if ((ref_no % 1000) == 0)
                AD_PRINT("{}", ref_no);
        }

        seqids_bins.erase(std::remove_if(std::begin(seqids_bins), std::end(seqids_bins), [](const auto& en) { return std::get<1>(en) == 0; }), std::end(seqids_bins));
        std::sort(std::begin(seqids_bins), std::end(seqids_bins), [](const auto& e1, const auto& e2) { return std::get<1>(e1) > std::get<1>(e2); });
        AD_INFO("Total selected: {}  With non-zero max bin: {}", refs_.size(), seqids_bins.size());
        for (const auto& [seq_id, max_bin, bins] : seqids_bins)
            AD_PRINT("  {:2d} {}  {}", max_bin, bins, seq_id);
    }
    return *this;

} // acmacs::seqdb::v3::subset::report_hamming_bins

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
