#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/uppercase.hh"

// ----------------------------------------------------------------------

namespace seqdb
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

          private:
            std::string json_text_;
            std::vector<SeqdbEntry> entries_;
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
            char aa_at(size_t pos0) const
            {
                const auto shift = string::from_chars<int>(a_shift);
                const auto shifted = static_cast<decltype(shift)>(pos0) - shift;
                if (shifted < 0)
                    return 'X';
                else if (static_cast<size_t>(shifted) >= amino_acids.size())
                    return '?';
                else
                    return amino_acids[static_cast<size_t>(shifted)];
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

            const auto& seq() const { return entry->seqs[seq_index]; }
            std::string seq_id() const { return string::join(" ", {entry->name, string::join(" ", seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()}); }
            bool has_lab(std::string_view lab) const { return seq().has_lab(lab); }
            bool has_clade(std::string_view clade) const { return seq().has_clade(clade); }
        };

        class subset
        {
          public:
            using refs_t = std::vector<ref>;
            using amino_acid_at_pos0_t = std::tuple<size_t, char, bool>; // pos (0-based), aa, equal/not-equal

            auto empty() const { return refs_.empty(); }
            auto size() const { return refs_.size(); }
            auto begin() const { return refs_.begin(); }
            auto end() const { return refs_.end(); }

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

            subset& print();

          private:
            refs_t refs_;

            subset() = default;

            void sort_by_date_recent_first() { std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.entry->date() > e2.entry->date(); }); }

            friend class Seqdb;
        };


    } // namespace v3
} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
