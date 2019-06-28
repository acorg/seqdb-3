#pragma once

#include <string>
#include <string_view>
#include <vector>

// ----------------------------------------------------------------------

namespace seqdb
{
    inline namespace v3
    {

        struct SeqdbSeq
        {
            using lab_ids_t = std::vector<std::string_view>;
            using labs_t = std::vector<std::pair<std::string_view, lab_ids_t>>;

            std::string_view amino_acids;
            std::string_view a_shift;
            std::string_view nucs;
            std::string_view n_shift;
            std::vector<std::string_view> passages;
            std::vector<std::string_view> reassortants;
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

        class Seqdb
        {
          public:
            Seqdb(const std::string& filename);
            // Seqdb(std::string&& source);

            using seq_ref_t = std::pair<const SeqdbEntry*, size_t>;
            using refs_t = std::vector<seq_ref_t>;

            refs_t select_by_name(std::string_view name) const;

          private:
            std::string json_text_;
            std::vector<SeqdbEntry> entries_;
        };

    } // namespace v3
} // namespace seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
