#pragma once

#include <vector>

#include "seqdb-3/fasta.hh"

namespace acmacs::seqdb
{
    inline namespace v3
    {
        using subtype_master_t = std::map<std::string, const sequence_t*, std::less<>>;

        subtype_master_t masters_per_subtype(const std::vector<fasta::scan_result_t>& sequences);

        void insertions_deletions(sequence_t& to_align, const sequence_t& master);
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
