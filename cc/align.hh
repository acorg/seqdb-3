#pragma once

#include <tuple>
#include <map>
#include <array>
#include <optional>
#include <string_view>

#include "acmacs-base/range-v3.hh"

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
            constexpr static size_t max_sequence_length{1000};
            constexpr static size_t number_of_symbols{128};
            constexpr static size_t table_size{number_of_symbols * max_sequence_length};

            struct table_t
            {
                table_t() { ranges::fill(data, 1); }

                // shift is non-positive!
                void update(std::string_view amino_acids, int shift);

                std::array<int, table_size> data;
            };

            std::map<std::string, table_t> tables_;
        };

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
