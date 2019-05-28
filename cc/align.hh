#pragma once

#include <tuple>
#include <optional>
#include <string_view>

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        std::optional<std::tuple<int, std::string_view>> align(std::string_view amino_acids, std::string_view type_subtype_hint, std::string_view debug_name);

        class Aligner
        {
          public:
            Aligner() = default;

            void update(std::string_view amino_acids, int shift, std::string_view type_subtype);
            std::optional<std::tuple<int, std::string_view>> align(std::string_view amino_acids, std::string_view type_subtype_hint, std::string_view debug_name) const;

          private:
        };

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
