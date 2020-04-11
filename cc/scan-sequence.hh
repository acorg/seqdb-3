#pragma once

#include <tuple>

#include "acmacs-base/date.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-strip.hh"
#include "acmacs-base/string-from-chars.hh"
#include "acmacs-virus/virus-name.hh"
#include "seqdb-3/types.hh"
#include "seqdb-3/sequence.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace scan
        {
            struct deletions_insertions_t
            {
                struct pos_num_t
                {
                    pos0_t pos;
                    size_t num;
                };

                std::vector<pos_num_t> deletions, insertions;

                bool empty() const { return deletions.empty() && insertions.empty(); }

                // returns {pos deleted, adjusted pos}
                std::pair<bool, pos0_t> aa_apply_deletions(pos0_t pos) const
                {
                    for (const auto& pos_num : deletions) {
                        if (pos_num.pos <= pos) {
                            if ((pos_num.pos + pos_num.num) > pos)
                                return {true, pos};
                            else
                                pos -= pos_num.num;
                        }
                        else
                            break;
                    }
                    return {false, pos};
                }

                // returns {pos deleted, adjusted pos}
                std::pair<bool, pos0_t> nuc_apply_deletions(pos0_t pos) const
                {
                    const auto [pos_deleted, adjusted_pos_aa] = aa_apply_deletions(pos.nuc_to_aa());
                    return {pos_deleted, pos0_t{(*adjusted_pos_aa * 3) + pos.nuc_offset()}};
                }

            }; // struct deletions_insertions_t

            std::string format_aa(const std::vector<deletions_insertions_t::pos_num_t>& pos_num, std::string_view sequence, char deletion_symbol = '-');
            // std::string format_nuc(const std::vector<deletions_insertions_t::pos_num_t>& pos_num, std::string_view sequence, char deletion_symbol = '-');
            std::string format(const deletions_insertions_t& deletions);

            inline std::string format_date(const date::year_month_day& a_date)
            {
                if (a_date.year().ok() && date::get_year(a_date) > 0)
                    return date::display(a_date, date::allow_incomplete::yes);
                else
                    return std::string{"0000-00-00"};
            }

            static const auto empty_month_or_day = [](const auto& date) { return date.find("-00") != std::string::npos; };
            static const auto not_empty_month_or_day = [](const auto& date) { return date.find("-00") == std::string::npos; };

            // ----------------------------------------------------------------------

            class sequence_t
            {
              public:
                sequence_t() = default;
                static sequence_t from_aligned_aa(const acmacs::virus::name_t& name, std::string_view source);

                using shift_t = named_size_t<struct shift_t_tag>;
                constexpr static const shift_t not_aligned{99999};

                void import(std::string_view source);
                void translate();

                std::string_view aa() const { return aa_; }
                constexpr auto aa_shift() const { return shift_aa_; }
                std::string_view nuc() const { return nuc_; }
                constexpr auto nuc_shift() const { return shift_nuc_; }
                constexpr const acmacs::virus::type_subtype_t& type_subtype() const { return type_subtype_; }
                constexpr const acmacs::virus::lineage_t& lineage() const { return lineage_; }
                std::optional<std::string> date() const
                {
                    if (dates_.empty())
                        return std::nullopt;
                    else
                        return dates_.front();
                }
                constexpr const auto& dates() const { return dates_; }
                std::string date_simulated() const noexcept; // returns either stored date or date inferred from name

                size_t year() const
                {
                    if (dates_.empty()) {
                        if (auto yr = acmacs::virus::year(name_); yr.has_value())
                            return *yr;
                        else
                            return 0;
                    }
                    else {
                        try {
                            return acmacs::string::from_chars<size_t>(dates_.front().substr(0, 4));
                        }
                        catch (std::exception&) {
                            fmt::print(stderr, "WARNING: cannot read year from {}\n", dates_.front());
                            return 0;
                        }
                    }
                }

                // aligned, deletions inserted
                std::string aa_format() const;
                std::string nuc_format() const;
                // NOT aligned, deletions inserted
                std::string aa_format_not_aligned() const;
                std::string nuc_format_not_aligned() const;

                // without applying deletions
                std::string_view aa_aligned() const
                {
                    std::string_view aa{aa_};
                    aa.remove_prefix(*shift_aa_);
                    return aa;
                }

                auto aa_aligned_length() const { return aa_.size() - *shift_aa_; }

                std::string_view aa_aligned_substr(size_t pos, size_t num) const { return std::string_view(aa_.data() + *shift_aa_ + pos, num); }

                // pos is 0 based
                // returns '-' if deleted
                // returns '\0' if beyond the end of sequence or before the beginning
                char aa_at_pos(pos0_t pos0) const
                {
                    if (const auto [deleted, pos_with_deletions] = deletions_.aa_apply_deletions(pos0); deleted)
                        return '-';
                    else if ((pos_with_deletions + *shift_aa_) > pos0_t{aa_.size()})
                        return 0;
                    else
                        return aa_[*pos_with_deletions + *shift_aa_];
                }

                // pos is 1 based
                char aa_at_pos(acmacs::seqdb::pos1_t pos) const { return aa_at_pos(pos0_t{*pos - 1}); }

                // pos is 0 based
                // returns '-' if deleted
                // returns '\0' if beyond the end of sequence or before the beginning
                char nuc_at_pos(pos0_t pos0) const
                {
                    const auto [deleted, pos_with_deletions] = deletions_.nuc_apply_deletions(pos0);
                    if (deleted)
                        return '-';
                    else if ((pos_with_deletions + *shift_nuc_) > pos0_t{nuc_.size()})
                        return 0;
                    else
                        return nuc_[*pos_with_deletions + *shift_nuc_];
                }

                char nuc_at_pos(acmacs::seqdb::pos1_t pos) const { return nuc_at_pos(pos0_t{*pos - 1}); }

                size_t aa_number_of_X() const
                {
                    if (aa_.empty())
                        throw std::runtime_error("internal in sequence_t::aa_number_of_X");
                    return static_cast<size_t>(std::count(aa_.begin() + static_cast<ssize_t>(*shift_aa_), aa_.end(), 'X'));
                }

                size_t aa_number_of_not_X() const { return aa_.size() - *shift_aa_ - aa_number_of_X(); }

                // without applying deletions
                std::string_view nuc_aligned() const
                {
                    std::string_view nuc{nuc_};
                    nuc.remove_prefix(*shift_nuc_);
                    return nuc;
                }

                constexpr const acmacs::virus::name_t& name() const { return name_; }
                std::string_view annotations() const { return annotations_; }
                std::string_view hash() const { return hash_; }
                // const acmacs::virus::host_t& host() const { return host_; }
                constexpr const acmacs::virus::Reassortant& reassortant() const { return reassortant_; }
                constexpr const auto& passages() const { return passages_; }
                acmacs::virus::Passage passage() const { return passages_.empty() ? acmacs::virus::Passage{} : passages_.front(); }
                std::string full_name() const;
                constexpr auto shift_aa() const { return shift_aa_; }
                constexpr auto shift_nuc() const { return shift_nuc_; }
                constexpr const clades_t& clades() const { return clades_; }
                std::string_view country() const { return country_; }
                std::string_view continent() const { return continent_; }
                constexpr const auto& hi_names() const { return hi_names_; }

                constexpr bool aligned() const { return shift_aa_ != not_aligned; }
                bool translated() const { return !aa_.empty(); }

                void set_shift(int shift_aa, std::optional<acmacs::virus::type_subtype_t> type_subtype = std::nullopt);

                void add_date(const std::string& date);
                void add_date(const date::year_month_day& a_date) { add_date(format_date(a_date)); }
                void remove_dates() { dates_.clear(); }
                void add_passage(acmacs::virus::Passage&& a_passage) { passages_.add(std::move(a_passage)); }
                void reassortant(const acmacs::virus::Reassortant& a_reassortant) { reassortant_ = a_reassortant; }
                void name(const acmacs::virus::name_t& a_name) { name_ = a_name; }
                void name(acmacs::virus::name_t&& a_name) { name_ = std::move(a_name); }
                // void host(acmacs::virus::host_t&& a_host) { host_ = std::move(a_host); }
                void annotations(std::string&& a_annotations) { annotations_ = std::move(a_annotations); }
                void remove_annotations() { annotations_.clear(); }
                void lineage(const acmacs::virus::lineage_t& lin) { lineage_ = lin; }
                void add_clade(const clade_t& clade) { clades_.add(clade); }
                void country(std::string_view country) { country_.assign(country); }
                void country(std::string&& country) { country_ = std::move(country); }
                void continent(std::string_view continent) { continent_.assign(continent); }
                void continent(std::string&& continent) { continent_ = std::move(continent); }
                void add_hi_name(const std::string& hi_name) { hi_names_.add(hi_name); }
                void update_subtype(const acmacs::virus::type_subtype_t& subtype);

                void add_lab_id(const acmacs::uppercase& lab, const acmacs::uppercase& lab_id);
                void add_lab_id(const acmacs::uppercase& lab);
                bool lab_in(std::initializer_list<std::string_view> labs) const;
                const auto& lab_ids() const { return lab_ids_; }

                constexpr deletions_insertions_t& deletions() { return deletions_; }
                constexpr const deletions_insertions_t& deletions() const { return deletions_; }

                // std::string name_with_annotations() const
                // {
                //     if (annotations_.empty())
                //         return *name_;
                //     else
                //         return fmt::format("{} {}", name_, annotations_);
                // }

                void add_isolate_id(std::string_view src) { isolate_id_.add(std::string{acmacs::string::strip(src)}); }
                void add_submitter(std::string_view src) { submitters_.add(std::string{acmacs::string::strip(src)}); }
                void add_sample_id_by_sample_provider(std::string_view src) { sample_id_by_sample_provider_.add(std::string{acmacs::string::strip(src)}); }
                void add_gisaid_last_modified(std::string_view src) { gisaid_last_modified_.add(std::string{acmacs::string::strip(src)}); }
                void add_originating_lab(std::string_view src) { originating_lab_.add(std::string{acmacs::string::strip(src)}); }
                void add_gisaid_segment(std::string_view src) { gisaid_segment_.add(std::string{acmacs::string::strip(src)}); }
                void add_gisaid_segment_number(std::string_view src) { gisaid_segment_number_.add(std::string{acmacs::string::strip(src)}); }
                void add_gisaid_identifier(std::string_view src) { gisaid_identifier_.add(std::string{acmacs::string::strip(src)}); }
                void add_gisaid_dna_accession_no(std::string_view src) { gisaid_dna_accession_no_.add(std::string{acmacs::string::strip(src)}); }
                void add_gisaid_dna_insdc(std::string_view src) { gisaid_dna_insdc_.add(std::string{acmacs::string::strip(src)}); }

                const auto& isolate_id() const { return isolate_id_; }
                const auto& submitters() const { return submitters_; }
                const auto& sample_id_by_sample_provider() const { return sample_id_by_sample_provider_; }
                const auto& gisaid_last_modified() const { return gisaid_last_modified_; }
                const auto& originating_lab() const { return originating_lab_; }
                const auto& gisaid_segment() const { return gisaid_segment_; }
                const auto& gisaid_segment_number() const { return gisaid_segment_number_; }
                const auto& gisaid_identifier() const { return gisaid_identifier_; }
                const auto& gisaid_dna_accession_no() const { return gisaid_dna_accession_no_; }
                const auto& gisaid_dna_insdc() const { return gisaid_dna_insdc_; }

                void merge_from(const sequence_t& source);

              private:
                acmacs::virus::name_t name_;
                // acmacs::virus::host_t host_;
                std::string country_;
                std::string continent_;
                flat_set_t<std::string> dates_;
                acmacs::virus::Reassortant reassortant_;
                flat_set_t<acmacs::virus::Passage> passages_;
                flat_set_t<std::string> hi_names_;
                std::string annotations_;
                std::map<acmacs::uppercase, std::set<acmacs::uppercase>> lab_ids_;
                std::string aa_;
                std::string nuc_;
                std::string hash_;
                int nuc_translation_offset_{0};
                shift_t shift_nuc_{not_aligned};
                shift_t shift_aa_{not_aligned};
                acmacs::virus::type_subtype_t type_subtype_; // by alignment
                deletions_insertions_t deletions_;
                acmacs::virus::lineage_t lineage_; // by deletion detection
                flat_set_t<std::string> isolate_id_;
                flat_set_t<std::string> submitters_;
                flat_set_t<std::string> sample_id_by_sample_provider_;
                flat_set_t<std::string> gisaid_last_modified_;
                flat_set_t<std::string> originating_lab_;
                flat_set_t<std::string> gisaid_segment_;
                flat_set_t<std::string> gisaid_segment_number_; // also ncbi segment number
                flat_set_t<std::string> gisaid_identifier_;
                flat_set_t<std::string> gisaid_dna_accession_no_;
                flat_set_t<std::string> gisaid_dna_insdc_;
                clades_t clades_;
                void aa_trim_absent(); // remove leading and trailing X and - from aa

            }; // class sequence_t

            inline auto designation(const sequence_t& seq)
            {
                return std::tuple{seq.name(), seq.annotations(), seq.reassortant()};
            }


        } // namespace scan

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
