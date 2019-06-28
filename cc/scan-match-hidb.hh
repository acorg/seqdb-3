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

            // sequences msut be sorted by name!
            void match_hidb(std::vector<fasta::scan_result_t>& sequences);

        } // namespace scan
    }     // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
