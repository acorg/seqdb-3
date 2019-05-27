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
            bool align(std::string_view type_subtype_hint, std::string_view debug_name); // returns if aligining succeeded

            const std::string_view aa() const { return aa_; }
            const std::string_view aa_aligned() const { return {aa_.data() + shift_aa_}; }
            const std::string_view nuc() const { return nuc_; }
            const std::string_view nuc_aligned() const { return {nuc_.data() + shift_nuc_}; }
            const std::string_view type_subtype() const { return type_subtype_; }
            const std::string_view lineage() const { return lineage_; }

            template <typename Int> void set_shift_aa(Int shift_aa) { shift_aa_ = static_cast<decltype(shift_aa_)>(shift_aa); shift_nuc_ = nuc_translation_offset_ + shift_aa_ * 3; }

          private:
            std::string aa_;
            std::string nuc_;
            int nuc_translation_offset_{0};
            std::string type_subtype_;
            std::string lineage_;
            int shift_nuc_{0};
            int shift_aa_{0};

            void import(std::string_view source);

            bool align_h3n2(std::string_view debug_name);
            bool align_any(std::string_view debug_name, std::string_view except = {});
        };
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
