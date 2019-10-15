#pragma once

#include <tuple>
#include <vector>

// ----------------------------------------------------------------------

namespace rjson::inline v2 { class value; }

namespace acmacs::seqdb::inline v3
{
    using amino_acid_at_pos_t = std::tuple<size_t, char, bool>; // pos, aa, equal/not-equal
    using amino_acid_at_pos_list_t = std::vector<amino_acid_at_pos_t>;

    amino_acid_at_pos_list_t extract_aa_at_pos(const rjson::value& source);
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
