#pragma once

#include <tuple>
#include <vector>
#include <algorithm>

#include "seqdb-3/sequence.hh"

// ----------------------------------------------------------------------

namespace rjson::v3 { class value; }

namespace acmacs::seqdb::inline v3
{
    class extract_at_pos_error : public error { public: using error::error; };

    // ======================================================================

    struct nucleotide_at_pos0_t : public std::tuple<pos0_t, char> // pos (0-based), aa
    {
        using std::tuple<pos0_t, char>::tuple;
    };
    using nucleotide_at_pos0_list_t = std::vector<nucleotide_at_pos0_t>;

    struct nucleotide_at_pos0_eq_t : public std::tuple<pos0_t, char, bool> // pos (0-based), aa, equal/not-equal
    {
        using std::tuple<pos0_t, char, bool>::tuple;
        constexpr nucleotide_at_pos0_eq_t() : std::tuple<pos0_t, char, bool>{pos0_t{0}, ' ', false} {}
    };

    using nucleotide_at_pos0_eq_list_t = std::vector<nucleotide_at_pos0_eq_t>;

    // --------------------------------------------------

    struct nucleotide_at_pos1_t : public std::tuple<pos1_t, char> // pos (1-based), aa
    {
        using std::tuple<pos1_t, char>::tuple;
    };
    using nucleotide_at_pos1_list_t = std::vector<nucleotide_at_pos1_t>;

    struct nucleotide_at_pos1_eq_t : public std::tuple<pos1_t, char, bool> // pos (1-based), aa, equal/not-equal
    {
        using std::tuple<pos1_t, char, bool>::tuple;
        constexpr nucleotide_at_pos1_eq_t() : std::tuple<pos1_t, char, bool>{pos1_t{0}, ' ', false} {}
    };

    using nucleotide_at_pos1_eq_list_t = std::vector<nucleotide_at_pos1_eq_t>;

    nucleotide_at_pos1_eq_list_t extract_nuc_at_pos1_eq_list(const rjson::v3::value& source);
    nucleotide_at_pos1_eq_t extract_nuc_at_pos1_eq(std::string_view source);
    nucleotide_at_pos1_eq_list_t extract_nuc_at_pos1_eq_list(std::string_view source); // space or comma separated list, e.g. "1703A 384C 618C !1010G"

    // ======================================================================

    struct amino_acid_at_pos0_t : public std::tuple<pos0_t, char> // pos (0-based), aa
    {
        using std::tuple<pos0_t, char>::tuple;
    };
    using amino_acid_at_pos0_list_t = std::vector<amino_acid_at_pos0_t>;

    struct amino_acid_at_pos0_eq_t : public std::tuple<pos0_t, char, bool> // pos (0-based), aa, equal/not-equal
    {
        using std::tuple<pos0_t, char, bool>::tuple;
        constexpr amino_acid_at_pos0_eq_t() : std::tuple<pos0_t, char, bool>{pos0_t{0}, ' ', false} {}
    };

    using amino_acid_at_pos0_eq_list_t = std::vector<amino_acid_at_pos0_eq_t>;

    // --------------------------------------------------

    struct amino_acid_at_pos1_t : public std::tuple<pos1_t, char> // pos (1-based), aa
    {
        using std::tuple<pos1_t, char>::tuple;
    };
    using amino_acid_at_pos1_list_t = std::vector<amino_acid_at_pos1_t>;

    struct amino_acid_at_pos1_eq_t : public std::tuple<pos1_t, char, bool> // pos (1-based), aa, equal/not-equal
    {
        using std::tuple<pos1_t, char, bool>::tuple;
        constexpr amino_acid_at_pos1_eq_t() : std::tuple<pos1_t, char, bool>{pos1_t{0}, ' ', false} {}
    };

    using amino_acid_at_pos1_eq_list_t = std::vector<amino_acid_at_pos1_eq_t>;

    // --------------------------------------------------

    using pos1_list_t = std::vector<pos1_t>;
    using pos0_list_t = std::vector<pos0_t>;
    pos1_list_t extract_pos1_list(std::string_view source);

    // --------------------------------------------------

    amino_acid_at_pos1_eq_list_t extract_aa_at_pos1_eq_list(const rjson::v3::value& source);
    amino_acid_at_pos1_eq_t extract_aa_at_pos1_eq(std::string_view source); // "!183P"
    amino_acid_at_pos1_eq_list_t extract_aa_at_pos1_eq_list(std::string_view source); // space or comma separated list, e.g. "183P 141E !123K"

    // ======================================================================

    inline bool matches(sequence_aligned_ref_t seq, const nucleotide_at_pos1_list_t& nuc_at_pos1) // all_of is not constexpr until c++20 (g++9.1)
    {
        return std::all_of(std::begin(nuc_at_pos1), std::end(nuc_at_pos1),
                           [seq](const nucleotide_at_pos1_t& pos1_nuc) { return at_pos(seq, std::get<pos1_t>(pos1_nuc)) == std::get<char>(pos1_nuc); });
    }

    inline bool matches(sequence_aligned_ref_t seq, const nucleotide_at_pos1_eq_list_t& nuc_at_pos1) // all_of is not constexpr until c++20 (g++9.1)
    {
        return std::all_of(std::begin(nuc_at_pos1), std::end(nuc_at_pos1),
                           [seq](const nucleotide_at_pos1_eq_t& pos1_nuc) { return (at_pos(seq, std::get<pos1_t>(pos1_nuc)) == std::get<char>(pos1_nuc)) == std::get<bool>(pos1_nuc); });
    }

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

    // ======================================================================

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    struct amino_acid_nucleotide_at_pos1_eq_t_tag {};
}

template <> struct fmt::formatter<acmacs::seqdb::amino_acid_nucleotide_at_pos1_eq_t_tag> : public fmt::formatter<acmacs::fmt_default_formatter>
{
    template <typename AA_NUC, typename FormatContext> auto format(const AA_NUC& pos1_aa_eq, FormatContext& ctx)
    {
        if (std::get<bool>(pos1_aa_eq))
            return format_to(ctx.out(), "{}{}", std::get<acmacs::seqdb::pos1_t>(pos1_aa_eq), std::get<char>(pos1_aa_eq));
        else
            return format_to(ctx.out(), "!{}{}", std::get<acmacs::seqdb::pos1_t>(pos1_aa_eq), std::get<char>(pos1_aa_eq));
    }
};

template <> struct fmt::formatter<acmacs::seqdb::amino_acid_at_pos1_eq_t> : public fmt::formatter<acmacs::seqdb::amino_acid_nucleotide_at_pos1_eq_t_tag>
{
    // template <typename FormatContext> auto format(const acmacs::seqdb::amino_acid_at_pos1_eq_t& pos1_aa_eq, FormatContext& ctx)
    // {
    //     return fmt::formatter<acmacs::seqdb::amino_acid_nucleotide_at_pos1_eq_t_tag>::format(pos1_aa_eq, ctx);
    // }
};

template <> struct fmt::formatter<acmacs::seqdb::nucleotide_at_pos1_eq_t> : public fmt::formatter<acmacs::seqdb::amino_acid_nucleotide_at_pos1_eq_t_tag>
{
    // template <typename FormatContext> auto format(const acmacs::seqdb::nucleotide_at_pos1_eq_t& pos1_aa_eq, FormatContext& ctx)
    // {
    //     return fmt::formatter<acmacs::seqdb::amino_acid_nucleotide_at_pos1_eq_t_tag>::format(pos1_aa_eq, ctx);
    // }
};

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
