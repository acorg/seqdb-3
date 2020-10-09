#pragma once

#include "acmacs-chart-2/log.hh"

// ----------------------------------------------------------------------

namespace acmacs::log
{
    enum {
        sequences = 12,
        fasta,
        hi_name_matching
    };

    inline void register_enabler_seqdb3()
    {
        using namespace std::string_view_literals;
        register_enabler_acmacs_chart();
        register_enabler("seq"sv, sequences);
        register_enabler("fasta"sv, fasta);
        register_enabler("matching"sv, hi_name_matching);
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
