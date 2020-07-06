#pragma once

#include <vector>

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
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
} // namespace acmacs::seqdb::inline v3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
