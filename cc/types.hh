#pragma once

#include "acmacs-base/named-type.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        using clade_t = named_string_t<struct clade_tag>;
        using clades_t = named_vector_t<clade_t>;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
