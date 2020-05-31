#pragma once

#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        using subsets_by_title_t = std::map<std::string, subset, std::less<>>;
        enum class compare {aa, nuc};

        std::string compare_report_text(const subset& sequences, enum compare cmp_nuc_aa, size_t split = 0);
        std::string compare_report_text(const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);
        std::string compare_report_text(const subset& set1, const subset& set2, enum compare cmp_nuc_aa);
        std::string compare_report_html(std::string_view title, const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);
        std::string compare_report_html(std::string_view title, const subset& sequences, enum compare cmp_nuc_aa, size_t split = 0);
        std::string compare_report_html(std::string_view title, const subset& set1, const subset& set2, enum compare cmp_nuc_aa);

        // ----------------------------------------------------------------------

        void update_common(sequence_aligned_t& target, const subset& source, enum compare cmp_nuc_aa);
        sequence_aligned_t find_common(const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);

        template <typename ... Subset> sequence_aligned_t find_common(const Subset& ... subsets)
        {
            sequence_aligned_t target;
            for (const auto& ss : {(subsets, ...)}) {
                update_common(target, ss);
            }
            return target;
        }


    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
