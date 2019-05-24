#pragma once

#include <string>

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        class sequence_t
        {
          public:
            sequence_t() = default;
            sequence_t(std::string_view source) { import(source); }
            void translate();

            const std::string_view aa() const { return aa_; }
            const std::string_view nuc() const { return nuc_; }

          private:
            std::string aa_;
            std::string nuc_;
            size_t nu_translation_offset_{0};

            void import(std::string_view source);
        };
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
