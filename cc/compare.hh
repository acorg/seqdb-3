#pragma once

#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        using subsets_by_title_t = std::map<std::string, subset>;

        std::string compare_report_text(const subset& sequences, size_t split = 0);
        std::string compare_report_text(const subsets_by_title_t& subsets);
        std::string compare_report_text(const subset& set1, const subset& set2);
        std::string compare_report_html(std::string_view title, const subsets_by_title_t& subsets);
        std::string compare_report_html(std::string_view title, const subset& sequences, size_t split = 0);
        std::string compare_report_html(std::string_view title, const subset& set1, const subset& set2);

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
