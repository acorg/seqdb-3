#pragma once

#include "acmacs-base/date.hh"
#include "acmacs-chart-2/passage.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        struct FastaEntry
        {
            std::string raw_name;
            std::string name;
            Date date;
            acmacs::chart::Passage passage;
            std::string lab_id;
            std::string lab;
            std::string virus_type;
            std::string sequence;

            FastaEntry() = default;
            FastaEntry(std::string_view rn, std::string seq) : raw_name(rn), sequence(seq) {}
            void parse_raw_name();
        };

        inline std::ostream& operator<<(std::ostream& out, const FastaEntry& entry) { return out << entry.raw_name; }
        
        std::vector<FastaEntry> fasta_scan(std::string_view data);
    }
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
