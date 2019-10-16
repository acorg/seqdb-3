#pragma once

#include <stdexcept>

namespace acmacs::seqdb::inline v3
{
    class error : public std::runtime_error
    {
      public:
        using std::runtime_error::runtime_error;
    };

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
