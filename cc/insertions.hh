#pragma once

#include <map>

#include "seqdb-3/fasta.hh"

namespace acmacs::seqdb
{
    inline namespace v3
    {
        using subtype_master_t = std::map<std::string, const sequence_t*, std::less<>>;

        subtype_master_t masters_per_subtype(const std::vector<fasta::scan_result_t>& sequences);

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
        deletions_insertions_t deletions_insertions(std::string_view master, std::string_view to_align);

        std::string format(const std::vector<deletions_insertions_t::pos_num_t>& pos_num, std::string_view sequence, char deletion_symbol = '-');

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
