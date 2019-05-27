#pragma once

#include <tuple>
#include <optional>
#include <vector>

#include "acmacs-base/date.hh"
#include "acmacs-virus/virus-name.hh"
#include "seqdb-3/sequence.hh"

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

            struct scan_options_t
            {
                size_t remove_too_short_nucs{100}; // remove nuc sequences shorter than this (if value 1000, sequence of length 1000 is kept)
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
                acmacs::virus::host_t host;
                Date date;
                acmacs::virus::Reassortant reassortant;
                acmacs::virus::Passage passage;
                std::string annotations;
                std::string lab_id;
                std::string lab;
                std::string type_subtype;
                std::string lineage;
                seqdb::sequence_t sequence;
            };

            struct scan_result_t
            {
                sequence_t seq;
                bool aligned{false};
                std::vector<acmacs::virus::v2::parse_result_t::message_t> messages;
                std::string filename;
                size_t line_no; // of the sequence name in filename
            };

            struct hint_t
            {
                std::string lab;
                std::string subtype;
                std::string lineage;
            };

            // ----------------------------------------------------------------------

            std::vector<scan_result_t> scan(const std::vector<std::string_view>& filenames, const scan_options_t& options);
            std::vector<std::reference_wrapper<scan_result_t>> aligned(std::vector<scan_result_t>& source);

            // ----------------------------------------------------------------------
            // details

            std::tuple<scan_input_t, scan_output_t> scan(scan_input_t input);

            std::optional<sequence_t> name_gisaid_spaces(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);
            std::optional<sequence_t> name_gisaid_underscores(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);
            std::optional<sequence_t> name_plain(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);

            using messages_t = std::vector<acmacs::virus::v2::parse_result_t::message_t>;

            // returns error and warning messages
            messages_t normalize_name(sequence_t& source);

            std::optional<acmacs::seqdb::sequence_t> import_sequence(std::string_view raw_sequence, const scan_options_t& options);

        } // namespace fasta

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
