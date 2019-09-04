#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/uppercase.hh"
#include "acmacs-base/flat-map.hh"

// ----------------------------------------------------------------------

namespace acmacs::chart { class Antigens; class Chart; class ChartModify; class PointIndexList; }

namespace acmacs::seqdb
{
    inline namespace v3
    {
        struct SeqdbSeq;
        struct SeqdbEntry;
        struct ref;
        class subset;

        using seq_id_index_t = flat_map_t<std::string, ref>;
        using hi_name_index_t = flat_map_t<std::string_view, ref>;

        class Seqdb
        {
          public:
            static const Seqdb& get();
            bool empty() const { return entries_.empty(); }

            subset all() const;
            subset select_by_seq_id(std::string_view seq_id) const;
            subset select_by_seq_id(const std::vector<std::string_view>& seq_ids) const;
            subset select_by_name(std::string_view name) const;
            subset select_by_name(const std::vector<std::string_view>& names) const;
            subset select_by_regex(std::string_view re) const;
            ref find_hi_name(std::string_view full_name) const;
            const seq_id_index_t& seq_id_index() const;
            const hi_name_index_t& hi_name_index() const;

            // returned subset contains elements for each antigen, i.e. it may contain empty ref's
            subset match(const acmacs::chart::Antigens& aAntigens, std::string_view aChartVirusType = {}) const;

            using aas_indexes_t = std::map<std::string, std::vector<size_t>>;
            aas_indexes_t aa_at_pos1_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions1) const;

            using clade_t = std::string_view;
            using clades_t = std::vector<clade_t>;
            enum class clades_for_name_inclusive { no /* only common clades for matching sequences */, yes /* all possible clades */ };
            clades_t clades_for_name(std::string_view name, clades_for_name_inclusive inclusive = clades_for_name_inclusive::no) const;

            void add_clades(acmacs::chart::ChartModify& chart) const;

            // returns subset where each entry corresponds to the entry in seq_ids
            subset find_by_seq_ids(const std::vector<std::string_view>& seq_ids) const;

            // returns json with data for ace-view/2018 sequences_of_chart command
            std::string sequences_of_chart_for_ace_view_1(const acmacs::chart::Chart& chart) const;
            // returns sequences in the fasta format
            std::string sequences_of_chart_as_fasta(const acmacs::chart::Chart& chart) const;

          private:
            std::string json_text_;
            std::vector<SeqdbEntry> entries_;
            mutable seq_id_index_t seq_id_index_;
            mutable hi_name_index_t hi_name_index_;

            Seqdb(std::string_view filename);
            // Seqdb(std::string&& source);

            void select_by_name(std::string_view name, subset& subs) const;

        };

        void setup(std::string_view filename);
        inline const Seqdb& get() { return Seqdb::get(); }

        // ----------------------------------------------------------------------

        struct export_options
        {
            enum class format { fasta_aa, fasta_nuc };

            format e_format{format::fasta_nuc};
            size_t e_wrap_at{0};
            bool e_aligned{true};
            bool e_most_common_length{false};
            std::string e_name_format{"{seq_id}"};

            export_options& fasta(bool nucs)
            {
                e_format = nucs ? format::fasta_nuc : format::fasta_aa;
                return *this;
            }
            export_options& wrap(size_t wrap_at)
            {
                e_wrap_at = wrap_at;
                return *this;
            }
            export_options& no_wrap()
            {
                e_wrap_at = 0;
                return *this;
            }
            export_options& aligned(bool aligned = true)
            {
                e_aligned = aligned;
                return *this;
            }
            export_options& most_common_length(bool most_common_length = true)
            {
                e_most_common_length = most_common_length;
                return *this;
            }
            export_options& name_format(std::string_view name_format)
            {
                e_name_format = name_format;
                return *this;
            }
        };

        // ----------------------------------------------------------------------

        struct aa_at_pos1_t
        {
            size_t pos1;
            char aa;
        };

        using list_aa_at_pos1_t = std::vector<aa_at_pos1_t>;

        list_aa_at_pos1_t parse_list_aa_at_pos1(std::string_view source); // space or comma separated list, e.g. "183P 141E"

        // ----------------------------------------------------------------------

        struct SeqdbSeq
        {
            using lab_ids_t = std::vector<std::string_view>;
            using labs_t = std::vector<std::pair<std::string_view, lab_ids_t>>;

            std::string_view amino_acids;
            std::string_view a_shift;
            std::string_view nucs;
            std::string_view n_shift;
            std::string_view annotations;
            std::vector<std::string_view> reassortants;
            std::vector<std::string_view> passages;
            std::vector<std::string_view> clades;
            std::vector<std::string_view> hi_names;
            labs_t lab_ids;

            bool has_lab(std::string_view lab) const
            {
                return std::any_of(std::begin(lab_ids), std::end(lab_ids), [lab](const auto& en) { return en.first == lab; });
            }
            bool has_clade(std::string_view clade) const { return std::find(std::begin(clades), std::end(clades), clade) != std::end(clades); }
            bool has_reassortant(std::string_view reassortant) const { return std::find(std::begin(reassortants), std::end(reassortants), reassortant) != std::end(reassortants); }
            bool match(const list_aa_at_pos1_t& aa_at_pos1) const;

            size_t aa_nuc_shift(std::string_view shift_s) const
            {
                if (shift_s.empty())
                    return 0;
                const auto shift = ::string::from_chars<int>(shift_s);
                if (shift > 0)
                    throw std::runtime_error(fmt::format("unsupported shift {}, hi_names: {}", shift, hi_names));
                return static_cast<size_t>(-shift);
            }

            std::string_view aa_aligned(size_t length = std::string_view::npos) const { return amino_acids.substr(aa_nuc_shift(a_shift), length); }
            std::string_view nuc_aligned(size_t length = std::string_view::npos) const { return nucs.substr(aa_nuc_shift(n_shift), length); }
            size_t aa_aligned_length() const { return amino_acids.size() - aa_nuc_shift(a_shift); }
            size_t nuc_aligned_length() const { return nucs.size() - aa_nuc_shift(n_shift); }

            char aa_at_pos0(size_t pos0) const
            {
                const auto aligned = aa_aligned();
                if (pos0 < aligned.size())
                    return aligned[pos0];
                else
                    return '?';
            }

            char aa_at_pos1(size_t pos1) const { return aa_at_pos0(pos1 - 1); }

            std::string_view lab() const { return lab_ids.empty() ? std::string_view{} : lab_ids.front().first; }
            std::string_view lab_id() const { return (lab_ids.empty() || lab_ids.front().second.empty()) ? std::string_view{} : lab_ids.front().second.front(); }
            std::string_view passage() const { return passages.empty() ? std::string_view{} : passages.front(); }
            std::string designation() const { return ::string::join(" ", {annotations, ::string::join(" ", reassortants), passage()}); }
        };

        struct SeqdbEntry
        {
            std::string_view name;
            std::string_view continent;
            std::string_view country;
            std::vector<std::string_view> dates;
            std::string_view lineage;
            std::string_view virus_type;
            std::vector<SeqdbSeq> seqs;

            std::string_view host() const;
            bool date_within(std::string_view start, std::string_view end) const { return !dates.empty() && (start.empty() || dates.front() >= start) && (end.empty() || dates.front() < end); }
            std::string_view date() const { return dates.empty() ? name.substr(name.size() - 4) : dates.front(); }
        };

        // ----------------------------------------------------------------------

        struct ref
        {
            const SeqdbEntry* entry;
            size_t seq_index;
            bool to_be_removed{false};  // for subsetting at random
            size_t group_no{0};         // for group_by_hamming_distance
            size_t hamming_distance{0}; // for group_by_hamming_distance

            ref() : entry{nullptr}, seq_index{static_cast<size_t>(-1)} {}
            ref(const SeqdbEntry* a_entry, size_t a_index) : entry{a_entry}, seq_index{a_index} {}

            constexpr bool operator==(const ref& rhs) const { return entry == rhs.entry && seq_index == rhs.seq_index; }
            constexpr bool operator!=(const ref& rhs) const { return !operator==(rhs); }
            constexpr explicit operator bool() const { return entry != nullptr; }

            const auto& seq() const { return entry->seqs[seq_index]; }
            std::string seq_id() const;
            std::string full_name() const { return ::string::join(" ", {entry->name, ::string::join(" ", seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()}); }
            std::string full_name_with_date() const
            {
                return fmt::format("{} [{}]", ::string::join(" ", {entry->name, ::string::join(" ", seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()}),
                                   entry->date());
            }
            std::string hi_name_or_full_name() const
            {
                if (seq().hi_names.empty())
                    return full_name();
                else
                    return std::string{seq().hi_names.front()};
            }
            bool has_lab(std::string_view lab) const { return seq().has_lab(lab); }
            bool has_clade(std::string_view clade) const { return seq().has_clade(clade); }
            bool has_hi_names() const { return !seq().hi_names.empty(); }
            bool match(const list_aa_at_pos1_t& aa_at_pos1) const { return seq().match(aa_at_pos1); }
        };

        class subset
        {
          public:
            using refs_t = std::vector<ref>;
            using amino_acid_at_pos0_t = std::tuple<size_t, char, bool>; // pos (0-based), aa, equal/not-equal
            enum class sorting { none, name_asc, name_desc, date_asc, date_desc };

            auto empty() const { return refs_.empty(); }
            auto size() const { return refs_.size(); }
            auto begin() const { return refs_.begin(); }
            auto end() const { return refs_.end(); }
            auto begin() { return refs_.begin(); }
            auto end() { return refs_.end(); }
            const auto& operator[](size_t index) const { return refs_.at(index); }
            const auto& front() const { return refs_.front(); }

            subset& multiple_dates(bool do_filter = true);
            subset& subtype(const acmacs::uppercase& virus_type);
            subset& lineage(const acmacs::uppercase& lineage);
            subset& lab(const acmacs::uppercase& lab);
            subset& whocc_lab(bool do_filter = true);
            subset& host(const acmacs::uppercase& host);
            subset& dates(std::string_view start, std::string_view end);
            subset& continent(const acmacs::uppercase& continent);
            subset& country(const acmacs::uppercase& country);
            subset& clade(const acmacs::uppercase& clade);
            subset& recent(size_t recent);
            subset& recent_matched(const std::vector<size_t>& recent_matched);
            subset& random(size_t random);
            subset& group_by_hamming_distance(size_t dist_threshold, size_t output_size);  // Eu's algorithm 2019-07-23
            subset& subset_by_hamming_distance_random(bool do_subset, size_t output_size); // davipatti algorithm 2019-07-23
            subset& remove_nuc_duplicates(bool do_remove, bool keep_hi_matched);
            subset& with_hi_name(bool with_hi_name);
            subset& aa_at_pos(const std::vector<amino_acid_at_pos0_t>& aa_at_pos0);
            subset& names_matching_regex(const std::vector<std::string_view>& regex_list);
            subset& names_matching_regex(std::string_view re) { return names_matching_regex(std::vector<std::string_view>{re}); }
            subset& prepend_single_matching(std::string_view re, const Seqdb& seqdb);
            subset& nuc_hamming_distance_to_base(size_t threshold, bool do_filter = true);
            subset& sort(sorting srt);
            subset& export_sequences(std::string_view filename, const export_options& options);
            subset& print(std::string_view name_format, bool do_print = true);
            subset& report_stat(bool do_report = true);

            subset& append(const subset& another);

            // returns new subset, this subset is not modified
            enum class matched_only { no, yes };
            subset filter_by_indexes(const acmacs::chart::PointIndexList& indexes, enum matched_only matched_only = matched_only::yes) const;

          private:
            refs_t refs_;

            subset() = default;
            subset(size_t size) : refs_(size) {}

            void sort_by_name_asc()
            {
                std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.seq_id() < e2.seq_id(); });
            }
            void sort_by_name_desc()
            {
                std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.seq_id() > e2.seq_id(); });
            }
            void sort_by_date_recent_first()
            {
                std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.entry->date() > e2.entry->date(); });
            }
            void sort_by_date_oldest_first()
            {
                std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.entry->date() < e2.entry->date(); });
            }
            void sort_by_hamming_distance()
            {
                std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.hamming_distance < e2.hamming_distance; });
            }

            refs_t::iterator most_recent_with_hi_name();
            void remove_marked()
            {
                refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.to_be_removed; }), std::end(refs_));
            }
            void set_remove_marker(bool marker)
            {
                std::for_each(std::begin(refs_), std::end(refs_), [marker](auto& en) { en.to_be_removed = marker; });
            }

            void resize(size_t size) { refs_.resize(size); }

            using collected_entry_t = std::pair<std::string, std::string>; // {seq_id, sequence}
            using collected_t = std::vector<collected_entry_t>;
            collected_t export_collect(const export_options& options) const;
            std::string export_fasta(const collected_t& entries, const export_options& options);
            std::string make_name(std::string_view name_format, const ref& entry) const;

            friend class Seqdb;
        };


    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
