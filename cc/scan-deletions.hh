#pragma once

#include <map>

#include "acmacs-base/debug.hh"
#include "seqdb-3/scan-fasta.hh"

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace scan
        {
            void detect_insertions_deletions(std::vector<fasta::scan_result_t>& sequence_data);

            // ----------------------------------------------------------------------

            void deletions_insertions(const sequence_t& master, sequence_t& to_align);
            deletions_insertions_t deletions_insertions(std::string_view master, std::string_view to_align, acmacs::debug dbg = acmacs::debug::no);

        } // namespace scan
    }     // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
