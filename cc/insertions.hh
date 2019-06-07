#pragma once

#include <vector>

#include "seqdb-3/fasta.hh"

namespace acmacs::seqdb
{
    inline namespace v3
    {
        std::map<std::string, const sequence_t*> masters_per_subtype(const std::vector<fasta::scan_result_t>& sequences);
        void insertions_deletions(std::vector<std::reference_wrapper<seqdb::sequence_t>>& sequences);
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
