#pragma once

#include <vector>
#include <string_view>

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace scan::fasta { struct scan_result_t; }

        void create(std::string_view prefix, std::vector<scan::fasta::scan_result_t>& sequences);

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
