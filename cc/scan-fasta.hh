#pragma once

#include <tuple>
#include <optional>
#include <vector>

#include "seqdb-3/scan-sequence.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace scan
        {
            namespace fasta
            {
                struct data_t
                {
                    std::string entry_name;
                    std::string name;
                    acmacs::virus::type_subtype_t type_subtype;
                    acmacs::virus::lineage_t lineage;
                    std::string passage;
                    std::string filename;
                    size_t line_no; // of the sequence name in filename
                    std::vector<acmacs::virus::parse_result_t::message_t> messages;
                };

                struct scan_result_t
                {
                    data_t fasta;
                    sequence_t sequence;
                };

                constexpr const auto is_aligned = [](const scan_result_t& sc) { return sc.sequence.aligned(); };
                constexpr const auto isnot_aligned = [](const scan_result_t& sc) { return !sc.sequence.aligned(); };
                constexpr const auto is_translated = [](const scan_result_t& sc) { return sc.sequence.translated(); };
                constexpr const auto is_different_type_subtype = [](const scan_result_t& sc) { return sc.fasta.type_subtype != sc.sequence.type_subtype(); };
                constexpr const auto is_different_type_subtype_ignore_h0 = [](const scan_result_t& sc) {
                    return sc.fasta.type_subtype != sc.sequence.type_subtype() && (sc.fasta.type_subtype.h_or_b() != "H0" || sc.sequence.type_subtype().type() != 'A');
                };

                // ----------------------------------------------------------------------

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

                struct hint_t
                {
                    std::string lab;
                    std::string subtype;
                    std::string lineage;
                };

                // ----------------------------------------------------------------------

                std::vector<scan_result_t> scan(const std::vector<std::string_view>& filenames, const scan_options_t& options);
                void sort_by_date(std::vector<fasta::scan_result_t>& sequences) noexcept;
                void sort_by_name(std::vector<fasta::scan_result_t>& sequences) noexcept;

                std::string report_false_positive(const std::vector<scan_result_t>& sequences,
                                                  size_t sequence_cutoff = std::string::npos); // report aligned having type_subtype that differs from provided with fasta
                std::string report_not_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_prefix, size_t sequence_cutoff = std::string::npos);
                std::string report_aa(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff = std::string::npos);
                std::string report_aa_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff = std::string::npos);

                // ----------------------------------------------------------------------
                // details

                std::tuple<scan_input_t, scan_output_t> scan(scan_input_t input);

                std::optional<scan_result_t> name_gisaid_spaces(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);
                std::optional<scan_result_t> name_gisaid_underscores(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);
                std::optional<scan_result_t> name_plain(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no);

                using messages_t = std::vector<acmacs::virus::v2::parse_result_t::message_t>;

                // returns error and warning messages
                messages_t normalize_name(scan_result_t& source);

                bool import_sequence(std::string_view raw_sequence, sequence_t& sequence_data, const scan_options_t& options);

            } // namespace fasta
        }     // namespace scan
    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
