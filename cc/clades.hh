#pragma once

#include <vector>

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace fasta { struct scan_result_t; }

        void detect_lineages_clades(std::vector<fasta::scan_result_t>& sequences);
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
