#pragma once

#include "acmacs-base/string-matcher.hh"
#include "acmacs-base/flat-set.hh"
#include "acmacs-virus/passage.hh"
#include "acmacs-virus/reassortant.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    struct score_size_t
    {
        string_match::score_t score;
        size_t len;

        constexpr bool operator<(const score_size_t& a) const { return score < a.score; }
    };

    struct score_seq_found_t : public score_size_t
    {
        size_t seq_no;
        size_t found_no;

        score_seq_found_t(const score_size_t& ss, size_t sn, size_t fn) : score_size_t{ss.score, ss.len}, seq_no{sn}, found_no{fn} {}
        constexpr bool operator<(const score_seq_found_t& a) const { return score > a.score; }
    };

    // ----------------------------------------------------------------------

    std::optional<score_size_t> match_antigen(const acmacs::virus::Reassortant& seq_reassortant, flat_set_t<acmacs::virus::Passage>& seq_passages, const acmacs::virus::Reassortant& hi_reassortant, acmacs::virus::Passage& hi_passage);

} // namespace acmacs::seqdb::inline v3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
