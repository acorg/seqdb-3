#pragma once

#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        using subsets_by_title_t = std::map<std::string, subset>;
        enum class compare {aa, nuc};

        std::string compare_report_text(const subset& sequences, enum compare cmp_nuc_aa, size_t split = 0);
        std::string compare_report_text(const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);
        std::string compare_report_text(const subset& set1, const subset& set2, enum compare cmp_nuc_aa);
        std::string compare_report_html(std::string_view title, const subsets_by_title_t& subsets, enum compare cmp_nuc_aa);
        std::string compare_report_html(std::string_view title, const subset& sequences, enum compare cmp_nuc_aa, size_t split = 0);
        std::string compare_report_html(std::string_view title, const subset& set1, const subset& set2, enum compare cmp_nuc_aa);

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
