#pragma once

#include <tuple>

#include "acmacs-base/named-type.hh"
#include "seqdb-3/error.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    struct pos0_t;

    struct pos1_t : public named_size_t<struct seqdb_pos1_tag_t>
    {
        using named_size_t<struct seqdb_pos1_tag_t>::named_size_t;
        constexpr pos1_t(pos0_t pos0);
    };

    struct pos0_t : public named_size_t<struct seqdb_pos0_tag_t>
    {
        using named_size_t<struct seqdb_pos0_tag_t>::named_size_t;
        constexpr pos0_t(pos1_t pos1) : named_size_t<struct seqdb_pos0_tag_t>{*pos1 - 1} {}
        constexpr pos0_t& operator=(const pos0_t&) = default;
        constexpr pos0_t& operator=(pos1_t pos1) { return operator=(pos0_t{*pos1 - 1}); }
        constexpr pos0_t nuc_to_aa() const { return pos0_t{get() / 3}; }
        constexpr pos0_t aa_to_nuc() const { return pos0_t{get() * 3}; }
        constexpr size_t nuc_offset() const { return get() % 3; }
    };

    constexpr inline pos1_t::pos1_t(pos0_t pos0) : named_size_t<struct seqdb_pos1_tag_t>{*pos0 + 1} {}

    template <typename P1, typename P2> using enable_if_different_pos_t = std::enable_if_t<(std::is_same_v<P1, pos0_t> || std::is_same_v<P1, pos1_t>) && (std::is_same_v<P2, pos0_t> || std::is_same_v<P2, pos1_t>) && !std::is_same_v<P1, P2>, char>;

    template <typename P1, typename P2, typename = enable_if_different_pos_t<P1, P2>> constexpr inline bool operator==(P1 p1, P2 p2) { return p1 == P1{p2}; }
    template <typename P1, typename P2, typename = enable_if_different_pos_t<P1, P2>> constexpr inline bool operator!=(P1 p1, P2 p2) { return !operator==(p1, p2); }
    template <typename P1, typename P2, typename = enable_if_different_pos_t<P1, P2>> constexpr inline bool operator>(P1 p1, P2 p2) { return p1 > P1{p2}; }
    template <typename P1, typename P2, typename = enable_if_different_pos_t<P1, P2>> constexpr inline bool operator>=(P1 p1, P2 p2) { return p1 >= P1{p2}; }
    template <typename P1, typename P2, typename = enable_if_different_pos_t<P1, P2>> constexpr inline bool operator<(P1 p1, P2 p2) { return p1 < P1{p2}; }
    template <typename P1, typename P2, typename = enable_if_different_pos_t<P1, P2>> constexpr inline bool operator<=(P1 p1, P2 p2) { return p1 <= P1{p2}; }

    // --------------------------------------------------

    struct sequence_aligned_t : public named_string_t<struct seqdb_sequence_aligned_ref_tag_t>
    {
        using base = named_string_t<struct seqdb_sequence_aligned_ref_tag_t>;
        using base::named_string_t;
        /*constexpr*/ char at(pos0_t pos0) const noexcept { return pos0 < size() ? operator[](*pos0) : ' '; }
        /*constexpr*/ pos0_t size() const noexcept { return pos0_t{base::size()}; }
        /*constexpr*/ void set(seqdb::pos0_t pos0, char aa) noexcept { if (pos0 < size()) get()[*pos0] = aa; }
        /*constexpr*/ void resize(seqdb::pos0_t new_size) { get().resize(*new_size); }
    };

    // not owning reference to a sequence
    struct sequence_aligned_ref_t : public named_string_view_t<struct seqdb_sequence_aligned_ref_tag_t>
    {
        using base = named_string_view_t<struct seqdb_sequence_aligned_ref_tag_t>;
        using base::named_string_view_t;
        /*constexpr*/ char at(pos0_t pos0) const noexcept { return pos0 < size() ? operator[](*pos0) : ' '; }
        /*constexpr*/ pos0_t size() const noexcept { return pos0_t{base::size()}; }
    };

    using alignment_t = named_int_from_string_t<struct seqdb_alignment_tag_t>;

    struct sequence_with_alignment_ref_t : public std::tuple<std::string_view, alignment_t>
    {
        using std::tuple<std::string_view, alignment_t>::tuple;
        constexpr sequence_with_alignment_ref_t() : std::tuple<std::string_view, alignment_t>{std::string_view{}, alignment_t{0}} {}
        bool empty() const noexcept { return std::get<std::string_view>(*this).empty(); }
        constexpr sequence_aligned_ref_t aligned(size_t length = std::string_view::npos) const;
    };

    inline size_t aligned_length(sequence_with_alignment_ref_t source) // std::abs() is not constexpr
    {
        // shift is negative in seqdb for historical reasons
        return std::get<std::string_view>(source).size() - static_cast<size_t>(std::abs(std::get<alignment_t>(source).as_number()));
    }

    constexpr inline sequence_aligned_ref_t aligned(sequence_with_alignment_ref_t source, size_t length = std::string_view::npos) noexcept
    {
        // shift is negative in seqdb for historical reasons
        return sequence_aligned_ref_t{std::get<std::string_view>(source).substr(static_cast<size_t>(std::abs(std::get<alignment_t>(source).as_number())), length)};
    }

    constexpr inline sequence_aligned_ref_t sequence_with_alignment_ref_t::aligned(size_t length) const { return acmacs::seqdb::aligned(*this, length); }

    inline char at_pos(sequence_aligned_ref_t seq, pos0_t pos0) noexcept { return seq.at(pos0); }
    inline char at_pos(sequence_aligned_ref_t seq, pos1_t pos1) noexcept { return seq.at(pos1); }
    inline char at_pos(sequence_with_alignment_ref_t seq, pos0_t pos0) noexcept { return aligned(seq).at(pos0); }
    inline char at_pos(sequence_with_alignment_ref_t seq, pos1_t pos1) noexcept { return aligned(seq).at(pos1); }

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------

template <> struct fmt::formatter<acmacs::seqdb::pos1_t> : fmt::formatter<size_t> {
    template <typename FormatCtx> auto format(const acmacs::seqdb::pos1_t& pos1, FormatCtx& ctx) { return fmt::formatter<size_t>::format(pos1.get(), ctx); }
};

template <> struct fmt::formatter<acmacs::seqdb::pos0_t> : fmt::formatter<size_t> {
    template <typename FormatCtx> auto format(const acmacs::seqdb::pos0_t& pos0, FormatCtx& ctx) { return fmt::formatter<size_t>::format(pos0.get() + 1, ctx); }
};

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
