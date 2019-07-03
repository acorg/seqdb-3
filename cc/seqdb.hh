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
        };

        // ----------------------------------------------------------------------

        struct ref
        {
            const SeqdbEntry* entry;
            size_t seq_index;

            ref(const SeqdbEntry* a_entry, size_t a_index) : entry{a_entry}, seq_index{a_index} {}
            std::string seq_id() const
            {
                const auto& seq = entry->seqs[seq_index];
                return string::join(" ", {entry->name, string::join(" ", seq.reassortants), seq.passages.empty() ? std::string_view{} : seq.passages.front()});
            }
        };

        class subset
        {
          public:
            using refs_t = std::vector<ref>;

            auto empty() const { return refs_.empty(); }
            auto size() const { return refs_.size(); }
            auto begin() const { return refs_.begin(); }
            auto end() const { return refs_.end(); }

            subset& multiple_dates();
            subset& subtype(const acmacs::uppercase& virus_type);

          private:
            refs_t refs_;

            subset() = default;

            friend class Seqdb;
        };


    } // namespace v3
} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
