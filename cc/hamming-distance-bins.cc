#include "acmacs-base/omp.hh"
#include "seqdb-3/hamming-distance-bins.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

namespace local
{
    using scan_result_iter = typename std::vector<acmacs::seqdb::v3::scan::fasta::scan_result_t>::const_iterator;

    constexpr const size_t BIN_SIZE = 200;

    static std::vector<size_t> hamming_distance_max_bin(scan_result_iter first, scan_result_iter last, size_t bin_size);
}

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::hamming_distance_bins_issues(std::vector<fasta::scan_result_t>& sequences)
{
    std::sort(std::begin(sequences), std::end(sequences), [](const auto& e1, const auto& e2) -> bool { return e1.fasta.type_subtype < e2.fasta.type_subtype; });

    auto first = std::begin(sequences);
    for (auto cur = std::next(first, 1); cur != std::end(sequences); ++cur) {
        if (cur->fasta.type_subtype != first->fasta.type_subtype) {
            AD_DEBUG("hamming_distance_bins_issues {} {} {}", cur - first, first->fasta.type_subtype, cur->fasta.type_subtype);
            const auto max_bin_per_seq = local::hamming_distance_bins_issues(first, cur, local::BIN_SIZE);
            first = cur;
        }
    }

    const auto max_bin_per_seq = local::hamming_distance_bins_issues(first, std::end(sequences), local::BIN_SIZE);

} // acmacs::seqdb::v3::scan::hamming_distance_bins_issues

// ----------------------------------------------------------------------

std::vector<size_t> local::hamming_distance_max_bin(scan_result_iter first, scan_result_iter last, size_t bin_size)
{
    std::vector<size_t> max_bin(static_cast<size_t>(last - first), 0);
    return max_bin;

} // local::hamming_distance_max_bin


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
