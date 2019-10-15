#pragma once

#include <tuple>
#include <vector>
#include <algorithm>

#include "seqdb-3/sequence.hh"

// ----------------------------------------------------------------------

namespace rjson::inline v2 { class value; }

namespace acmacs::seqdb::inline v3
{
    using amino_acid_at_pos1_t = std::tuple<pos1_t, char, bool>; // pos (1-based), aa, equal/not-equal
    using amino_acid_at_pos1_list_t = std::vector<amino_acid_at_pos1_t>;

    amino_acid_at_pos1_list_t extract_aa_at_pos1(const rjson::value& source);

    constexpr inline bool matches(sequence_aligned_ref_t seq, const amino_acid_at_pos1_list_t& aa_at_pos1)
    {
        return std::all_of(std::begin(aa_at_pos1), std::end(aa_at_pos1),
                           [seq](const amino_acid_at_pos1_t& pos1_aa) { return (at_pos(seq, std::get<pos1_t>(pos1_aa)) == std::get<char>(pos1_aa)) == std::get<bool>(pos1_aa); });
    }

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
