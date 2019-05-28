#pragma once

#include <string>

#include "acmacs-base/date.hh"
#include "acmacs-virus/virus-name.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        class sequence_t
        {
            using shift_t = int;
            constexpr static const shift_t not_aligned = -99999;

          public:
            sequence_t() = default;

            void import(std::string_view source);
            void translate();
            bool align(std::string_view type_subtype_hint, std::string_view debug_name); // returns if aligining succeeded

            std::string_view aa() const { return aa_; }
            std::string_view nuc() const { return nuc_; }
            std::string_view type_subtype() const { return type_subtype_; }
            std::string_view lineage() const { return lineage_; }

            std::string aa_aligned() const
            {
                if (shift_aa_ > 0)
                    return aa_.substr(static_cast<size_t>(shift_aa_));
                else if (shift_aa_ == 0)
                    return aa_;
                else
                    return std::string(static_cast<size_t>(- shift_aa_), 'X') + aa_;
            }

            // returned shift is non-negative
            // if it is 0, aligned sequence is returned
            // if it is negative, sequence must be prepended with shift 'X's
            std::tuple<std::string_view, shift_t> aa_shifted() const
            {
                if (shift_aa_ >= 0)
                    return {std::string_view(aa_).substr(static_cast<size_t>(shift_aa_)), 0};
                else
                    return {std::string_view(aa_), shift_aa_};
            }

            std::string nuc_aligned() const
            {
                if (shift_nuc_ > 0)
                    return nuc_.substr(static_cast<size_t>(shift_nuc_));
                else if (shift_nuc_ == 0)
                    return nuc_;
                else
                    return std::string(static_cast<size_t>(- shift_nuc_), '-') + nuc_;
            }

            const acmacs::virus::virus_name_t& name() const { return name_; }
            // const acmacs::virus::host_t& host() const { return host_; }
            std::string_view annotations() const { return annotations_; }
            std::string_view lab_id() const { return lab_id_; }
            std::string_view lab() const { return lab_; }

            constexpr bool aligned() const { return shift_aa_ != not_aligned; }
            bool translated() const { return !aa_.empty(); }

            void set_shift(int shift_aa, std::optional<std::string_view> type_subtype = std::nullopt)
            {
                shift_aa_ = shift_aa;
                shift_nuc_ = nuc_translation_offset_ + shift_aa_ * 3;
                if (type_subtype.has_value())
                    type_subtype_ = *type_subtype;
            }

            void date(const Date& a_date) { date_ = a_date; }
            void passage(acmacs::virus::Passage&& a_passage) { passage_ = std::move(a_passage); }
            void reassortant(const acmacs::virus::Reassortant& a_reassortant) { reassortant_ = a_reassortant; }
            void lab_id(std::string&& a_lab_id) { lab_id_ = std::move(a_lab_id); }
            void lab_id(const std::string& a_lab_id) { lab_id_ = a_lab_id; }
            void lab(std::string_view a_lab) { lab_ = a_lab; }
            void name(acmacs::virus::virus_name_t&& a_name) { name_ = std::move(a_name); }
            // void host(acmacs::virus::host_t&& a_host) { host_ = std::move(a_host); }
            void annotations(std::string&& a_annotations) { annotations_ = std::move(a_annotations); }
            void remove_annotations() { annotations_.clear(); }
            // void (const & a_) { _ = a_; }

          private:
            acmacs::virus::virus_name_t name_;
            // acmacs::virus::host_t host_;
            Date date_;
            acmacs::virus::Reassortant reassortant_;
            acmacs::virus::Passage passage_;
            std::string annotations_;
            std::string lab_id_;
            std::string lab_;
            std::string aa_;
            std::string nuc_;
            int nuc_translation_offset_{0};
            shift_t shift_nuc_{not_aligned};
            shift_t shift_aa_{not_aligned};
            std::string type_subtype_; // by alignment
            std::string lineage_;      // by deletion detection
        };

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
