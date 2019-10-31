#pragma once

#include "acmacs-base/named-type.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    using seq_id_t = acmacs::named_string_t<struct acmacs_seqdb_SeqId_tag>;

    seq_id_t make_seq_id(std::string_view designation);
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
