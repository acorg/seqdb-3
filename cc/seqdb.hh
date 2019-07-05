#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/uppercase.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        struct SeqdbSeq;
        struct SeqdbEntry;
        class subset;

        class Seqdb
        {
          public:
            Seqdb(const std::string& filename);
            // Seqdb(std::string&& source);

            subset all() const;
            subset select_by_name(std::string_view name) const;
            subset select_by_regex(std::string_view re) const;

          private:
            std::string json_text_;
            std::vector<SeqdbEntry> entries_;
        };

        // ----------------------------------------------------------------------

        struct export_options
        {
            enum class format { fasta_aa, fasta_nuc };

            format e_format{format::fasta_nuc};
            size_t e_wrap_at{0};
            bool   e_aligned{true};

            export_options& fasta(bool nucs) { e_format = nucs ? format::fasta_nuc : format::fasta_aa; return *this; }
            export_options& wrap(size_t wrap_at) { e_wrap_at = wrap_at; return *this; }
            export_options& no_wrap() { e_wrap_at = 0; return *this; }
            export_options& aligned(bool aligned = true) { e_aligned = aligned; return *this; }

        };

        // ----------------------------------------------------------------------

        struct SeqdbSeq
        {
            using lab_ids_t = std::vector<std::string_view>;
            using labs_t = std::vector<std::pair<std::string_view, lab_ids_t>>;

            std::string_view amino_acids;
            std::string_view a_shift;
            std::string_view nucs;
            std::string_view n_shift;
            std::vector<std::string_view> reassortants;
            std::vector<std::string_view> passages;
            std::vector<std::string_view> clades;
            std::vector<std::string_view> hi_names;
            labs_t lab_ids;

            bool has_lab(std::string_view lab) const { return std::any_of(std::begin(lab_ids), std::end(lab_ids), [lab](const auto& en) { return en.first == lab; }); }
            bool has_clade(std::string_view clade) const { return std::find(std::begin(clades), std::end(clades), clade) != std::end(clades); }

            size_t aa_nuc_shift(std::string_view shift_s) const
            {
                if (shift_s.empty())
                    return 0;
                const auto shift = ::string::from_chars<int>(shift_s);
                if (shift > 0)
                    throw std::runtime_error(fmt::format("unsupported shift {}, hi_names: {}", shift, hi_names));
                return static_cast<size_t>(- shift);
            }

            std::string_view aa_aligned() const { return amino_acids.substr(aa_nuc_shift(a_shift)); }
            std::string_view nuc_aligned() const { return nucs.substr(aa_nuc_shift(n_shift)); }

            char aa_at(size_t pos0) const
            {
                const auto aligned = aa_aligned();
                if (pos0 < aligned.size())
                    return aligned[pos0];
                else
                    return '?';
            }
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
            bool selected{false}; // for choosing at random

            ref(const SeqdbEntry* a_entry, size_t a_index) : entry{a_entry}, seq_index{a_index} {}

            constexpr bool operator==(const ref& rhs) const { return entry == rhs.entry && seq_index == rhs.seq_index; }
            constexpr bool operator!=(const ref& rhs) const { return !operator==(rhs); }

            const auto& seq() const { return entry->seqs[seq_index]; }
            std::string seq_id() const;
            std::string full_name() const { return ::string::join(" ", {entry->name, ::string::join(" ", seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()}); }
            bool has_lab(std::string_view lab) const { return seq().has_lab(lab); }
            bool has_clade(std::string_view clade) const { return seq().has_clade(clade); }
        };

        class subset
        {
          public:
            using refs_t = std::vector<ref>;
            using amino_acid_at_pos0_t = std::tuple<size_t, char, bool>; // pos (0-based), aa, equal/not-equal
            enum class print_options { details, seq_id };

            constexpr auto empty() const { return refs_.empty(); }
            constexpr auto size() const { return refs_.size(); }
            constexpr auto begin() const { return refs_.begin(); }
            constexpr auto end() const { return refs_.end(); }
            constexpr const auto& front() const { return refs_.front(); }

            subset& multiple_dates(bool do_filter = true);
            subset& subtype(const acmacs::uppercase& virus_type);
            subset& lineage(const acmacs::uppercase& lineage);
            subset& lab(const acmacs::uppercase& lab);
            subset& host(const acmacs::uppercase& host);
            subset& dates(std::string_view start, std::string_view end);
            subset& continent(const acmacs::uppercase& continent);
            subset& country(const acmacs::uppercase& country);
            subset& clade(const acmacs::uppercase& clade);
            subset& recent(size_t recent);
            subset& random(size_t random);
            subset& with_hi_name(bool with_hi_name);
            subset& aa_at_pos(const std::vector<amino_acid_at_pos0_t>& aa_at_pos0);
            subset& names_matching_regex(const std::vector<std::string_view>& regex_list);
            subset& names_matching_regex(std::string_view re) { return names_matching_regex(std::vector<std::string_view>{re}); }
            subset& prepend_single_matching(std::string_view re, const Seqdb& seqdb);
            subset& nuc_hamming_distance_to_base(size_t threshold, bool do_filter = true);
            subset& export_sequences(std::string_view filename, const export_options& options);
            subset& print(print_options po, bool do_print = true);

          private:
            refs_t refs_;

            subset() = default;

            void sort_by_date_recent_first() { std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.entry->date() > e2.entry->date(); }); }
            void export_fasta(const ref& entry, const export_options& options, std::string& output);

            friend class Seqdb;
        };


    } // namespace v3
} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
