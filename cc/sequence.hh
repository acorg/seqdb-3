#pragma once

#include <tuple>

#include "acmacs-base/named-type.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    using pos0_t = named_size_t<struct seqdb_pos0_tag_t>;
    using pos1_t = named_size_t<struct seqdb_pos1_tag_t>;

    using sequence_aligned_t = named_string_t<struct seqdb_sequence_aligned_tag_t>;
    using sequence_aligned_ref_t = named_string_view_t<struct seqdb_sequence_aligned_ref_tag_t>;

    using alignment_t = named_int_from_string_t<struct seqdb_alignment_tag_t>;
    using sequence_with_alignment_ref_t = std::tuple<std::string_view, alignment_t>;

    constexpr sequence_aligned_ref_t aligned(sequence_with_alignment_ref_t source, size_t length = std::string_view::npos)
    {
        // shift is negative in seqdb for historical reasons
        return sequence_aligned_ref_t{std::get<std::string_view>(source).substr(static_cast<size_t>(std::abs(std::get<alignment_t>(source).as_number())), length)};
    }

    constexpr char at_pos(sequence_aligned_ref_t seq, pos0_t pos0) { return seq->at(*pos0); }
    constexpr char at_pos(sequence_aligned_ref_t seq, pos1_t pos1) { return seq->at(*pos1 - 1); }
    constexpr char at_pos(sequence_with_alignment_ref_t seq, pos0_t pos0) { return aligned(seq)->at(*pos0); }
    constexpr char at_pos(sequence_with_alignment_ref_t seq, pos1_t pos1) { return aligned(seq)->at(*pos1 - 1); }

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
