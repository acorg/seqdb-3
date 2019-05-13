#pragma once

#include "acmacs-base/date.hh"
#include "acmacs-chart-2/passage.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        struct source_ref_t
        {
            std::string_view filename;
            size_t line_no = 0;
        };

        inline std::ostream& operator<<(std::ostream& out, const source_ref_t& source_ref) { return out << source_ref.filename << ':' << source_ref.line_no; }
        
        class FastaEntry
        {
          public:
            std::string raw_name;
            std::string name;
            Date date;
            acmacs::chart::Passage passage;
            std::string lab_id;
            std::string lab;
            std::string virus_type;
            std::string sequence;
            source_ref_t source_ref;

            FastaEntry() = default;
            FastaEntry(std::string_view rn, std::string_view seq, source_ref_t a_source_ref) : raw_name(rn), sequence(seq), source_ref(a_source_ref) {}

          private:
            bool normalize_sequence();
            bool parse();
            bool name_gisaid_spaces(std::string_view source);
            bool name_gisaid_underscores(std::string_view source);
            void parse_date(std::string_view source);
            friend std::vector<FastaEntry> fasta_scan(std::string_view filename, std::string_view data);
        };

        std::ostream& operator<<(std::ostream& out, const FastaEntry& entry);

        std::vector<FastaEntry> fasta_scan(std::string_view filename);
        std::vector<FastaEntry> fasta_scan(std::string_view filename, std::string_view data);
    }
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
