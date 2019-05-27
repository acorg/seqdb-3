#pragma once

namespace acmacs::seqdb
{
    inline namespace v3
    {
        template <typename S1, typename S2> inline size_t hamming_distance(S1&& s1, S2&& s2)
        {
            auto f1 = std::begin(s1);
            auto f2 = std::begin(s2);
            size_t dist = 0;
            for (; f1 != std::end(s1) && f2 != std::end(s2); ++f1, ++f2) {
                if (*f1 != *f2)
                    ++dist;
            }
            dist += static_cast<size_t>(std::end(s1) - f1) + static_cast<size_t>(std::end(s2) - f2);
            return dist;
        }

        template <typename S1, typename S2, typename C> inline size_t hamming_distance_not_considering(S1&& s1, S2&& s2, C not_consider)
        {
            const auto f1 = std::begin(s1);
            const auto f2 = std::begin(s2);
            size_t dist = 0;
            for (; f1 != std::end(s1) && f2 != std::end(s2); ++f1, ++f2) {
                if (*f1 != *f2 && *f1 != not_consider && *f2 != not_consider)
                    ++dist;
            }
            dist += static_cast<size_t>(std::end(s1) - f1) + static_cast<size_t>(std::end(s2) - f2);
            return dist;
        }

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
