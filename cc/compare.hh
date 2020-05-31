#pragma once

#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    enum class compare { aa, nuc };

    struct subset_to_compare_t
    {
        std::string name;
        acmacs::seqdb::subset subset;

        subset_to_compare_t(std::string_view a_name) : name{a_name} {}
    };

    using subsets_to_compare_t = std::vector<subset_to_compare_t>;

    struct difference_t
    {
        pos0_t pos;
    };

    using differences_t = std::vector<difference_t>;

    differences_t find_differences(const subsets_to_compare_t& subsets, enum compare cmp_nuc_aa);

    // ----------------------------------------------------------------------

    void update_common(sequence_aligned_t& target, const subset& source, enum compare cmp_nuc_aa);
    sequence_aligned_t find_common(const subsets_to_compare_t& subsets, enum compare cmp_nuc_aa);

    template <typename... Subset> sequence_aligned_t find_common(const Subset&... subsets)
    {
        sequence_aligned_t target;
        for (const auto& ss : {(subsets, ...)}) {
            update_common(target, ss);
        }
        return target;
    }

} // namespace acmacs::seqdb::inline v3

// ----------------------------------------------------------------------
// old API

namespace acmacs::seqdb::inline v3
{
    using subsets_by_title_t = std::map<std::string, subset, std::less<>>;

    std::string compare_report_text(const subset& sequences, enum compare cmp_nuc_aa, size_t split = 0);
    std::string compare_report_text(const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);
    std::string compare_report_text(const subset& set1, const subset& set2, enum compare cmp_nuc_aa);
    std::string compare_report_html(std::string_view title, const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);
    std::string compare_report_html(std::string_view title, const subset& sequences, enum compare cmp_nuc_aa, size_t split = 0);
    std::string compare_report_html(std::string_view title, const subset& set1, const subset& set2, enum compare cmp_nuc_aa);

    // ----------------------------------------------------------------------


} // namespace acmacs::seqdb::inline v3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
