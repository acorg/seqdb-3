#pragma once

#include <bitset>

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

    using issues_t = std::bitset<static_cast<size_t>(issue::_last)>;

    inline void set(issues_t& issues, issue iss) { issues.set(static_cast<size_t>(iss)); }
    inline void reset(issues_t& issues, issue iss) { issues.reset(static_cast<size_t>(iss)); }
    inline bool has(const issues_t& issues, issue iss) { return issues[static_cast<size_t>(iss)]; }

} // namespace acmacs::seqdb::inline v3::sequence

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
