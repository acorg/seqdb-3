#pragma once

#include <vector>
#include <string_view>

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace fasta { struct scan_result_t; }

        enum class create_filter { all_aligned, h1_h3_b_aligned, whocc_aligned };

        void create(std::string_view filename, std::vector<fasta::scan_result_t>& sequences, create_filter filter);

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
