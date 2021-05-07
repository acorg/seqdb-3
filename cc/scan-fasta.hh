#pragma once

#include <optional>
#include <vector>

#include "acmacs-base/log.hh"
#include "acmacs-virus/virus-name-normalize.hh"
#include "seqdb-3/scan-sequence.hh"

// ----------------------------------------------------------------------

namespace acmacs::messages::inline v1
{
    namespace key
    {
        constexpr static inline std::string_view ncbi_unrecognized{"ncbi-unrecognized"};
        constexpr static inline std::string_view ncbi_unrecognized_subtype{"ncbi-unrecognized-subtype"};
        constexpr static inline std::string_view ncbi_unrecognized_fna_name{"ncbi-unrecognized-fna-name"};
        constexpr static inline std::string_view ncbi_dat_fna_name_difference{"ncbi-dat-fna-name-difference"};

        constexpr static inline std::string_view gisaid_invalid_date{"gisaid-invalid-date"};
        constexpr static inline std::string_view gisaid_invalid_subtype{"gisaid-invalid-subtype"};

        constexpr static inline std::string_view fasta_no_year_at_the_end_of_name{"fasta-no-year-at-the-end-of-name"};
        constexpr static inline std::string_view fasta_name_contains_annotations{"fasta-name-contains-annotations"};
        constexpr static inline std::string_view fasta_country_name_mismatch{"fasta-country-name-mismatch"};

    }
}

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
                    acmacs::uppercase passage;
                    std::string country; // ncbi
                    std::string filename;
                    size_t line_no; // of the sequence name in filename
                    acmacs::messages::messages_t messages;
                };

                struct master_ref_t
                {
                    acmacs::virus::name_t name;
                    std::string hash;
                    // std::string annotations;
                    // acmacs::virus::Reassortant reassortant;
                    // acmacs::uppercase passage;
                };

                struct scan_result_t
                {
                    data_t fasta;
                    sequence_t sequence;
                    std::optional<master_ref_t> reference; // master entry with identical nuc sequence and nuc_shift
                    bool remove{false};
                    acmacs::virus::name::parsed_fields_t name_fields; // for merging dat and fan names for ncbi
                };

                struct scan_results_t
                {
                    std::vector<scan_result_t> results;
                    acmacs::messages::messages_t messages;

                    void merge(scan_results_t&& source);
                };

                constexpr const auto is_aligned = [](const scan_result_t& sc) { return sc.sequence.aligned(); };
                constexpr const auto isnot_aligned = [](const scan_result_t& sc) { return !sc.sequence.aligned(); };
                constexpr const auto is_translated = [](const scan_result_t& sc) { return sc.sequence.translated(); };
                constexpr const auto is_different_type_subtype = [](const scan_result_t& sc) { return !(sc.fasta.type_subtype == sc.sequence.type_subtype()); };
                constexpr const auto is_different_type_subtype_ignore_h0 = [](const scan_result_t& sc) {
                    const auto f_hb = sc.fasta.type_subtype.h_or_b(), s_hb = sc.sequence.type_subtype().h_or_b();
                    return f_hb != s_hb && (f_hb != "H0" || s_hb == "B");
                };

                // ----------------------------------------------------------------------

                struct scan_error : public std::runtime_error
                {
                    using std::runtime_error::runtime_error;
                };

                struct manually_excluded : public std::runtime_error
                {
                    manually_excluded(std::string_view msg) : std::runtime_error{std::string{msg}} {}
                };

                enum class scan_name_adjustments { none, gisaid, ncbi };
                enum class print_names { no, yes };

                struct scan_options_t
                {
                    scan_options_t(debug a_dbg, scan_name_adjustments na = scan_name_adjustments::none, print_names pn = print_names::no) : dbg{a_dbg}, name_adjustements{na}, prnt_names{pn} {}
                    scan_options_t(size_t a_remove_too_short_nucs, debug a_dbg = debug::no, scan_name_adjustments na = scan_name_adjustments::none)
                        : remove_too_short_nucs{a_remove_too_short_nucs}, dbg{a_dbg}, name_adjustements{na} {}
                    size_t remove_too_short_nucs{100}; // remove nuc sequences shorter than this (if value 1000, sequence of length 1000 is kept)
                    debug dbg;
                    scan_name_adjustments name_adjustements{scan_name_adjustments::none};
                    print_names prnt_names{print_names::no};
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
                    acmacs::uppercase lab;
                    acmacs::uppercase subtype;
                    acmacs::uppercase lineage;
                };

                // ----------------------------------------------------------------------

                scan_results_t scan(const std::vector<std::string_view>& filenames, const scan_options_t& options);
                scan_results_t scan_ncbi(std::string_view directory, const scan_options_t& options);

                inline void sort_by_date(std::vector<fasta::scan_result_t>& sequences) noexcept
                {
                    std::sort(std::begin(sequences), std::end(sequences), [](const auto& e1, const auto& e2) { return e1.sequence.date_simulated() < e2.sequence.date_simulated(); });
                }

                inline void sort_by_name(std::vector<fasta::scan_result_t>& sequences) noexcept
                {
                    std::sort(std::begin(sequences), std::end(sequences), [](const auto& e1, const auto& e2) { return designation(e1.sequence) < designation(e2.sequence); });
                }

                void remove_without_names(std::vector<fasta::scan_result_t>& sequences);
                void merge_duplicates(std::vector<fasta::scan_result_t>& sequences);

                struct min_max_dates_t
                {
                    std::string min_isolation_date, max_isolation_date, min_submission_date, max_submission_date;
                };

                min_max_dates_t min_max_dates(const std::vector<scan_result_t>& sequences);

                std::string report_false_positive(const std::vector<scan_result_t>& sequences, size_t sequence_cutoff = std::string::npos); // report aligned having type_subtype that differs from provided with fasta
                std::string report_not_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_prefix, size_t sequence_cutoff = std::string::npos);
                std::string report_aa(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff = std::string::npos);
                std::string report_aa_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff = std::string::npos);

                // ----------------------------------------------------------------------
                // details

                std::tuple<scan_input_t, scan_output_t> scan(scan_input_t input);

                // returns error and warning messages
                acmacs::messages::messages_t normalize_name(scan_result_t& source, debug dbg, scan_name_adjustments name_adjustements, print_names prnt_names);
                void fix_gisaid_name(scan_result_t& source, acmacs::messages::messages_t& messages, debug dbg);
                std::string fix_ncbi_name(std::string_view source, acmacs::messages::messages_t& messages, debug dbg);
                // date::year_month_day parse_date(const acmacs::uppercase& source, std::string_view filename, size_t line_no);
                bool import_sequence(std::string_view raw_sequence, sequence_t& sequence_data, const scan_options_t& options);

            } // namespace fasta
        }     // namespace scan
    } // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
