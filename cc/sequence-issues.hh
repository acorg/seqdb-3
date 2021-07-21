#pragma once

#include <bitset>
#include "acmacs-base/fmt.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3::sequence
{
    enum class issue {
        not_aligned,
        has_insertions, // not detected, unused
        too_short,
        garbage_at_the_beginning,
        garbage_at_the_end,
        high_hamming_distance_bin,

        _last
    };

    constexpr const auto number_of_issues{static_cast<size_t>(issue::_last)};
    constexpr const std::array<const char*, number_of_issues> issue_name{"Not aligned", "Has insertions", "Too short", "garbage_at_the_beginning", "garbage_at_the_end", "high_hamming_distance_bin"};
    constexpr const std::array<char, number_of_issues> issue_name_char{'A', 'i', 's', 'b', 'e', 'h'};

    using issues_t = std::bitset<number_of_issues>;

    inline void set(issues_t& issues, issue iss) { issues.set(static_cast<size_t>(iss)); }
    inline void reset(issues_t& issues, issue iss) { issues.reset(static_cast<size_t>(iss)); }
    inline bool has(const issues_t& issues, issue iss) { return issues[static_cast<size_t>(iss)]; }

} // namespace acmacs::seqdb::inline v3::sequence

template <> struct fmt::formatter<acmacs::seqdb::v3::sequence::issues_t> : fmt::formatter<acmacs::fmt_helper::default_formatter> {
    template <typename FormatCtx> auto format(const acmacs::seqdb::v3::sequence::issues_t& issues, FormatCtx& ctx)
    {
        using namespace acmacs::seqdb::v3::sequence;
        bool first{true};
        for (auto iss = static_cast<size_t>(issue::not_aligned); iss < number_of_issues; ++iss) {
            if (has(issues, static_cast<issue>(iss))) {
                if (!first)
                    format_to(ctx.out(), " ");
                else
                    first = false;
                format_to(ctx.out(), "{}", issue_name[iss]);
            }
        }
        return ctx.out();
    }
};

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
