#pragma once

#include "acmacs-chart-2/log.hh"

// ----------------------------------------------------------------------

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#endif

namespace acmacs::log::inline v1
{
    const log_key_t sequences{"seq"};
    const log_key_t fasta{"fasta"};
    const log_key_t hi_name_matching{"seqdb-matching"};

} // namespace acmacs::log::inline v1

#pragma GCC diagnostic pop

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
