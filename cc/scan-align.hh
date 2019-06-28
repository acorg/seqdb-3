#pragma once

#include <vector>

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace scan
        {
            namespace fasta
            {
                struct scan_result_t;
            }

            // removes not translated
            void translate_align(std::vector<fasta::scan_result_t>& sequences);
        } // namespace scan

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
