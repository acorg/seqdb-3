#pragma once

#include <string_view>
#include <vector>


// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        struct SeqdbEntry;

        void parse(std::string_view source, std::vector<SeqdbEntry>& entries_);

    } // namespace v3
} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
