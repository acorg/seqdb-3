#include "acmacs-base/omp.hh"
#include "acmacs-base/timeit.hh"
#include "seqdb-3/hamming-distance-bins.hh"
#include "seqdb-3/hamming-distance.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3::scan::local
{
    using scan_result_iter = typename std::vector<acmacs::seqdb::v3::scan::fasta::scan_result_t>::iterator;

    constexpr const size_t BIN_SIZE = 200;
    constexpr const size_t MIN_BIN = 1;
    constexpr const size_t MAX_BINS = 2000 / BIN_SIZE + 1;

    // returns {with_issue, total_checked}
    static std::pair<size_t, size_t> set_high_hamming_distance_bin_issue(scan_result_iter first, scan_result_iter last, size_t bin_size, size_t min_bin);
    static std::vector<size_t> hamming_distance_max_bin(scan_result_iter first, scan_result_iter last, size_t bin_size);

} // namespace acmacs::seqdb::inline v3::scan::local

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::hamming_distance_bins_issues(std::vector<fasta::scan_result_t>& sequences)
{
    // 2021-07-21
    // for each sequence that has no issues find hamming distances to all other sequences without issues of the same subtype (for H1 consider only H1pdm)
    // put calculated hamming distances into bins (BIN_SIZE)
    // if a bin number with maximum number of sequences >= MIN_BIN, add high_hamming_distance_bin issue for that sequence

    std::sort(std::begin(sequences), std::end(sequences), [](const auto& e1, const auto& e2) -> bool {
        // sort by subtype, issues and year,
        // having issues first (for a subtype), they will be ignored
        // year is needed to ignore H1 sequences before 2009
        const std::tuple t1{e1.fasta.type_subtype, e1.sequence.good(), e1.sequence.year()}, t2{e2.fasta.type_subtype, e2.sequence.good(), e2.sequence.year()};
        return t1 < t2;
    });

    const auto h1_h3_b = [](auto subtype) {
        using namespace acmacs::virus;
        return subtype == type_subtype_t{"A(H1N1)"} || subtype == type_subtype_t{"A(H3N2)"} || subtype == type_subtype_t{"B"};
    };

    const auto skip_with_issues = [](auto first, auto last) {
        while (!first->sequence.good() && first != last)
            ++first;
        return first;
    };

    const auto check_high_hamming_distance_bin_issue = [h1_h3_b, skip_with_issues](auto first, auto last) {
        if (h1_h3_b(first->fasta.type_subtype)) {
            Timeit ti {fmt::format("Hamming distance bin issues {}", first->fasta.type_subtype)};
            const auto [with_issue, total] = local::set_high_hamming_distance_bin_issue(skip_with_issues(first, last), last, local::BIN_SIZE, local::MIN_BIN);
            ti.message_append(fmt::format(" with-issue:{} total:{}", with_issue, total));
        }
    };

    auto first = std::begin(sequences);
    for (auto cur = std::next(first, 1); cur != std::end(sequences); ++cur) {
        if (cur->fasta.type_subtype != first->fasta.type_subtype) {
            check_high_hamming_distance_bin_issue(first, cur);
            first = cur;
        }
    }

    check_high_hamming_distance_bin_issue(first, std::end(sequences));

} // acmacs::seqdb::v3::scan::hamming_distance_bins_issues

// ----------------------------------------------------------------------

std::pair<size_t, size_t> acmacs::seqdb::v3::scan::local::set_high_hamming_distance_bin_issue(scan_result_iter first, scan_result_iter last, size_t bin_size, size_t min_bin)
{
    if (first->fasta.type_subtype == acmacs::virus::type_subtype_t{"A(H1N1)"}) {
        // skip sequences before 2009
        while (first->sequence.year() < 2009 && first != last)
            ++first;
    }

    const auto max_bin_per_seq = hamming_distance_max_bin(first, last, bin_size);
    size_t with_issue = 0;
    for (auto mb = max_bin_per_seq.begin(); mb != max_bin_per_seq.end(); ++mb) {
        if (*mb >= min_bin) {
            std::next(first, mb - max_bin_per_seq.begin())->sequence.add_issue(sequence::issue::high_hamming_distance_bin);
            // AD_DEBUG("{:60s} {}", std::next(first, mb - max_bin_per_seq.begin())->sequence.name(), std::next(first, mb - max_bin_per_seq.begin())->sequence.issues());
            ++with_issue;
        }
    }

    return {with_issue, static_cast<size_t>(last - first)};

} // acmacs::seqdb::v3::scan::local::set_high_hamming_distance_bin_issue

// ----------------------------------------------------------------------

std::vector<size_t> acmacs::seqdb::v3::scan::local::hamming_distance_max_bin(scan_result_iter first, scan_result_iter last, size_t bin_size)
{
#ifdef _OPENMP
    const int num_threads = omp_get_max_threads();
    // const int slot_size = max_bin.size() / num_threads;
#endif

    const auto num_sequences = static_cast<size_t>(last - first);

    // get all aligned nuc sequences
    std::vector<std::string> nucs(num_sequences);
#pragma omp parallel for default(shared) num_threads(num_threads)
    for (auto nucp = nucs.begin(); nucp != nucs.end(); ++nucp)
        *nucp = std::next(first, nucp - nucs.begin())->sequence.nuc_format();

    // compute pairwise distances between sequences
    using dist_t = uint16_t;
    std::vector<dist_t> distances(num_sequences * num_sequences);
#pragma omp parallel for default(shared) num_threads(num_threads)
    for (size_t s1 = 0; s1 < num_sequences; ++s1) {
        for (size_t s2 = s1 + 1; s2 < num_sequences; ++s2)
            distances[s1 * num_sequences + s2] = hamming_distance<dist_t>(nucs[s1], nucs[s2], hamming_distance_by_shortest::yes);
    }

    std::vector<size_t> max_bin(static_cast<size_t>(num_sequences), 0);

#pragma omp parallel for default(shared) num_threads(num_threads)
    for (size_t s1 = 0; s1 < num_sequences; ++s1) {
        std::array<size_t, MAX_BINS> bins;
        bins.fill(0ul);
        const auto set_bin = [&bins, bin_size](size_t dist) {
            if (dist > 0)
                ++bins[dist / bin_size];
        };
        for (size_t s2 = 0; s2 < s1; ++s2)
            set_bin(distances[s2 * num_sequences + s1]); // s2 < s1
        for (size_t s2 = s1 + 1; s2 < num_sequences; ++s2)
            set_bin(distances[s1 * num_sequences + s2]); // s2 > s1
        // AD_DEBUG("bins {:60s} {}", std::next(first, static_cast<ssize_t>(s1))->sequence.name(), bins);
        max_bin[s1] = static_cast<size_t>(std::max_element(std::begin(bins), std::end(bins)) - std::begin(bins));
        // AD_DEBUG(max_bin[s1] > 0, "bins {:60s} {}", std::next(first, static_cast<ssize_t>(s1))->sequence.name(), bins);
    }

    return max_bin;

} // acmacs::seqdb::v3::scan::local::hamming_distance_max_bin


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
