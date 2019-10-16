#pragma once

#include <tuple>
#include <vector>
#include <algorithm>

#include "seqdb-3/sequence.hh"

// ----------------------------------------------------------------------

namespace rjson::inline v2 { class value; }

namespace acmacs::seqdb::inline v3
{
    using amino_acid_at_pos0_t = std::tuple<pos0_t, char>; // pos (0-based), aa
    using amino_acid_at_pos0_list_t = std::vector<amino_acid_at_pos0_t>;

    struct amino_acid_at_pos0_eq_t : public std::tuple<pos0_t, char, bool> // pos (0-based), aa, equal/not-equal
    {
        using std::tuple<pos0_t, char, bool>::tuple;
        constexpr amino_acid_at_pos0_eq_t() : std::tuple<pos0_t, char, bool>{pos0_t{0}, ' ', false} {}
    };

    using amino_acid_at_pos0_eq_list_t = std::vector<amino_acid_at_pos0_eq_t>;

    // --------------------------------------------------

    using amino_acid_at_pos1_t = std::tuple<pos1_t, char>; // pos (1-based), aa
    using amino_acid_at_pos1_list_t = std::vector<amino_acid_at_pos1_t>;

    struct amino_acid_at_pos1_eq_t : public std::tuple<pos1_t, char, bool> // pos (1-based), aa, equal/not-equal
    {
        using std::tuple<pos1_t, char, bool>::tuple;
        constexpr amino_acid_at_pos1_eq_t() : std::tuple<pos1_t, char, bool>{pos1_t{0}, ' ', false} {}
    };

    using amino_acid_at_pos1_eq_list_t = std::vector<amino_acid_at_pos1_eq_t>;

    // --------------------------------------------------

    amino_acid_at_pos1_eq_list_t extract_aa_at_pos1_eq_list(const rjson::value& source);
    amino_acid_at_pos1_eq_t extract_aa_at_pos1_eq(std::string_view source); // "!183P"
    amino_acid_at_pos1_eq_list_t extract_aa_at_pos1_eq_list(std::string_view source); // space or comma separated list, e.g. "183P 141E !123K"
    inline amino_acid_at_pos1_eq_list_t extract_aa_at_pos1_eq_list(const std::string& source) { return extract_aa_at_pos1_eq_list(std::string_view{source}); }

    inline bool matches(sequence_aligned_ref_t seq, const amino_acid_at_pos1_list_t& aa_at_pos1) // all_of is not constexpr until c++20 (g++9.1)
    {
        return std::all_of(std::begin(aa_at_pos1), std::end(aa_at_pos1),
                           [seq](const amino_acid_at_pos1_t& pos1_aa) { return at_pos(seq, std::get<pos1_t>(pos1_aa)) == std::get<char>(pos1_aa); });
    }

    inline bool matches(sequence_aligned_ref_t seq, const amino_acid_at_pos1_eq_list_t& aa_at_pos1) // all_of is not constexpr until c++20 (g++9.1)
    {
        return std::all_of(std::begin(aa_at_pos1), std::end(aa_at_pos1),
                           [seq](const amino_acid_at_pos1_eq_t& pos1_aa) { return (at_pos(seq, std::get<pos1_t>(pos1_aa)) == std::get<char>(pos1_aa)) == std::get<bool>(pos1_aa); });
    }

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------

template <> struct fmt::formatter<acmacs::seqdb::amino_acid_at_pos1_eq_t>
{
    template <typename ParseContext> constexpr auto parse(ParseContext& ctx) { return ctx.begin(); }
    template <typename FormatContext> auto format(const acmacs::seqdb::amino_acid_at_pos1_eq_t& pos1_aa_eq, FormatContext& ctx)
    {
        if (std::get<bool>(pos1_aa_eq))
            return format_to(ctx.out(), "{}{}", std::get<acmacs::seqdb::pos1_t>(pos1_aa_eq), std::get<char>(pos1_aa_eq));
        else
            return format_to(ctx.out(), "!{}{}", std::get<acmacs::seqdb::pos1_t>(pos1_aa_eq), std::get<char>(pos1_aa_eq));
    }
};

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
