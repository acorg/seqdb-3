#pragma once

#include <map>

#include "acmacs-base/debug.hh"
#include "seqdb-3/fasta.hh"

namespace acmacs::seqdb
{
    inline namespace v3
    {

        void detect_insertions_deletions(std::vector<fasta::scan_result_t>& sequence_data);

        // ----------------------------------------------------------------------

        struct deletions_insertions_t
        {
            struct pos_num_t
            {
                size_t pos;
                size_t num;
            };

            std::vector<pos_num_t> deletions, insertions;
        };

        void deletions_insertions(const sequence_t& master, sequence_t& to_align);
        deletions_insertions_t deletions_insertions(std::string_view master, std::string_view to_align, acmacs::debug dbg = acmacs::debug::no);

        std::string format(const std::vector<deletions_insertions_t::pos_num_t>& pos_num, std::string_view sequence, char deletion_symbol = '-');
        std::string format(const deletions_insertions_t& deletions);

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
