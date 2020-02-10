#pragma once

#include "seqdb-3/sequence.hh"

namespace acmacs::seqdb::inline v3
{
    enum class hamming_distance_by_shortest { no, yes };

    inline size_t hamming_distance(std::string_view s1, std::string_view s2, hamming_distance_by_shortest shortest = hamming_distance_by_shortest::no)
    {
        auto f1 = std::begin(s1);
        auto f2 = std::begin(s2);
        size_t dist = 0;
        for (; f1 != std::end(s1) && f2 != std::end(s2); ++f1, ++f2) {
            if (*f1 != *f2)
                ++dist;
        }
        if (shortest == hamming_distance_by_shortest::no)
            dist += static_cast<size_t>(std::end(s1) - f1) + static_cast<size_t>(std::end(s2) - f2);
        return dist;
    }

    inline size_t hamming_distance(sequence_aligned_ref_t s1, sequence_aligned_ref_t s2)
    {
        return hamming_distance(*s1, *s2);
    }

    inline size_t hamming_distance(sequence_with_alignment_ref_t s1, sequence_with_alignment_ref_t s2)
    {
        return hamming_distance(s1.aligned(), s2.aligned());
    }

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
