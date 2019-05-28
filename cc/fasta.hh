#pragma once

#include <tuple>
#include <optional>
#include <vector>

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

            struct data_t
            {
                std::string entry_name;
                std::string name;
                std::string type_subtype;
                std::string lineage;
                std::string passage;
                std::string filename;
                size_t line_no; // of the sequence name in filename
                std::vector<acmacs::virus::parse_result_t::message_t> messages;
            };

            struct scan_result_t
            {
                data_t fasta;
                seqdb::sequence_t sequence;
            };

            struct hint_t
            {
                std::string lab;
                std::string subtype;
                std::string lineage;
            };

            // ----------------------------------------------------------------------

            std::vector<scan_result_t> scan(const std::vector<std::string_view>& filenames, const scan_options_t& options);
            // std::vector<std::reference_wrapper<scan_result_t>> aligned(std::vector<scan_result_t>& source);

            // ----------------------------------------------------------------------
            // details

            std::tuple<scan_input_t, scan_output_t> scan(scan_input_t input);

            std::optional<scan_result_t> name_gisaid_spaces(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);
            std::optional<scan_result_t> name_gisaid_underscores(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);
            std::optional<scan_result_t> name_plain(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);

            using messages_t = std::vector<acmacs::virus::v2::parse_result_t::message_t>;

            // returns error and warning messages
            messages_t normalize_name(scan_result_t& source);

            bool import_sequence(std::string_view raw_sequence, seqdb::sequence_t& sequence_data, const scan_options_t& options);

        } // namespace fasta

    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
