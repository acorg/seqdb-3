#pragma once

#include <tuple>
#include <optional>

#include "acmacs-base/date.hh"
#include "acmacs-virus/virus-name.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace fasta
        {
            struct scan_error : public std::runtime_error
            {
                using std::runtime_error::runtime_error;
            };

            struct scan_input_t
            {
                std::string_view::const_iterator first;
                std::string_view::const_iterator last;
                size_t line_no = 1;
                size_t name_line_no = 1;
                constexpr bool done() const { return first == last; }
            };

            struct scan_output_t
            {
                std::string_view name;
                std::string_view sequence;
            };

            struct sequence_t
            {
                std::string fasta_name;
                std::string raw_name;
                acmacs::virus::virus_name_t name;
                Date date;
                acmacs::virus::Reassortant reassortant;
                acmacs::virus::Passage passage;
                std::string annotations;
                std::string lab_id;
                std::string lab;
                std::string virus_type;
                std::string lineage;
                std::string sequence;
            };

            // ----------------------------------------------------------------------

            std::tuple<scan_input_t, scan_output_t> scan(scan_input_t input);

            std::optional<sequence_t> name_gisaid_spaces(std::string_view name, std::string_view filename, size_t line_no);
            std::optional<sequence_t> name_gisaid_underscores(std::string_view name, std::string_view filename, size_t line_no);
            std::optional<sequence_t> name_plain(std::string_view name, std::string_view filename, size_t line_no);

            // returns error and warning messages
            std::vector<acmacs::virus::v2::parse_result_t::message_t> normalize_name(sequence_t& source);

            std::string normalize_sequence(std::string_view raw_sequence);

        } // namespace fasta

        // struct source_ref_t
        // {
        //     std::string_view filename;
        //     size_t line_no = 0;
        // };

        // inline std::ostream& operator<<(std::ostream& out, const source_ref_t& source_ref) { return out << source_ref.filename << ':' << source_ref.line_no; }

        // class FastaEntry
        // {
        //   public:
        //     std::string raw_name;
        //     std::string name;
        //     Date date;
        //     acmacs::virus::Passage passage;
        //     std::string lab_id;
        //     std::string lab;
        //     std::string virus_type;
        //     std::string sequence;
        //     source_ref_t source_ref;

        //     FastaEntry() = default;
        //     FastaEntry(std::string_view rn, std::string_view seq, source_ref_t a_source_ref) : raw_name(rn), sequence(seq), source_ref(a_source_ref) {}

        //   private:
        //     bool normalize_sequence();
        //     bool parse();
        //     bool name_gisaid_spaces(std::string_view source);
        //     bool name_gisaid_underscores(std::string_view source);
        //     void parse_date(std::string_view source);
        //     friend std::vector<FastaEntry> fasta_scan(std::string_view filename, std::string_view data);
        // };

        // std::ostream& operator<<(std::ostream& out, const FastaEntry& entry);

        // std::vector<FastaEntry> fasta_scan(std::string_view filename);
        // std::vector<FastaEntry> fasta_scan(std::string_view filename, std::string_view data);
    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
