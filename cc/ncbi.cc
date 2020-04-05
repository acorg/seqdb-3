#include "acmacs-base/read-file.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

using cursor_t = decltype(std::string{}.cbegin());

inline cursor_t token(cursor_t cursor, cursor_t end)
{
    for (; cursor != end; ++cursor) {
        switch (*cursor) {
            case '\t':
            case '\n':
                return cursor;
        }
    }
    return cursor;
};

// ----------------------------------------------------------------------

inline std::string_view string_view(cursor_t begin, cursor_t end) { return std::string_view(&*begin, static_cast<size_t>(end - begin)); }

enum class na_field : int { genbank_accession = 0, host, segment_no, subtype, country, date, sequence_length, virus_name, age, gender, completeness };
constexpr inline na_field& operator++(na_field& fld) { return fld = static_cast<na_field>(static_cast<int>(fld) + 1); }

// ----------------------------------------------------------------------

inline acmacs::virus::type_subtype_t parse_subtype(const acmacs::uppercase& source, std::string_view filename, size_t line_no)
{
    if (source.empty())
        AD_WARNING("no subtype @@ {}:{}", filename, line_no);
    if (source.size() >= 5 && source->substr(0, 5) == "MIXED") {
        if (source.size() == 5)
            return acmacs::virus::type_subtype_t{};
        if (source.size() == 6 || source[6] == '.' || source[6] == ',') {
            AD_WARNING("unrecognized subtype: \"{}\" @@ {}:{}", source, filename, line_no);
            return acmacs::virus::type_subtype_t{};
        }
        return acmacs::virus::type_subtype_t{fmt::format("A({})", source->substr(6))};
    }
    if (source.size() > 7 && source->substr(0, 5) == "MIXED,")
        return acmacs::virus::type_subtype_t{fmt::format("A({})", source->substr(6))};
    if (source.size() < 2 || source.size() > 6) // H3 N3 H13 H3N2 H9N1,9
        AD_WARNING("unrecognized subtype: \"{}\" @@ {}:{}", source, filename, line_no);
    return acmacs::virus::type_subtype_t{fmt::format("A({})", source)};
}

// ----------------------------------------------------------------------

inline std::string fix_country(std::string_view source)
{
    using pp = std::pair<std::string_view, std::string_view>;
    using namespace std::string_view_literals;
    static std::array country_mapping{
        pp{"USA"sv,                              "UNITED STATES OF AMERICA"sv},
        pp{"DEMOCRATIC REPUBLIC OF THE CONGO"sv, "CONGO DEMOCRATIC REPUBLIC"sv},
        pp{"VIET NAM"sv,                         "VIETNAM"sv},
        pp{"COTE D'IVOIRE"sv,                    "IVORY COAST"sv},
        pp{"LAB"sv,                              ""sv}, // error in ncbi database?
        // pp{sv,                                sv},
    };

    for (const auto& [from, to] : country_mapping) {
        if (source == from)
            return std::string(to);
    }

    return std::string(source);
}

// ----------------------------------------------------------------------

inline std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> influenza_na_read_entry(cursor_t& cur, cursor_t end, std::string_view filename, size_t line_no)
{
    acmacs::seqdb::v3::scan::fasta::scan_result_t result;
    result.fasta.filename = filename;
    result.fasta.line_no = line_no;

    na_field field{na_field::genbank_accession};
    cursor_t segment_number;
    for (auto tok_beg = cur, tok_end = token(tok_beg, end); tok_end != end; tok_beg = std::next(tok_end), tok_end = token(tok_beg, end), ++field) {
        if (tok_beg != tok_end) {
            switch (field) {
                case na_field::genbank_accession:
                    result.sequence.add_sample_id_by_sample_provider(string_view(tok_beg, tok_end));
                    break;
                case na_field::segment_no:
                    segment_number = tok_beg;
                    result.sequence.add_gisaid_segment_number(string_view(tok_beg, tok_end));
                    break;
                case na_field::virus_name:
                    result.fasta.name = string_view(tok_beg, tok_end);
                    break;
                case na_field::subtype:
                    result.fasta.type_subtype = parse_subtype(string_view(tok_beg, tok_end), filename, line_no);
                    break;
                case na_field::date:
                    break;
                case na_field::country:
                    result.fasta.country = fix_country(string::upper(tok_beg, tok_end));
                    break;
                case na_field::sequence_length:
                case na_field::age:
                case na_field::host:
                case na_field::gender:
                case na_field::completeness:
                    break;
            }
        }

        const auto make_result = [&]() -> std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> {
            if (*segment_number == '4') // interested in segment 4 (HA) only
                return result;
            else
                return std::nullopt;
        };

        if (tok_end == end) {
            cur = tok_end;
            return make_result();
        }
        else if (*tok_end == '\n') {
            cur = std::next(tok_end);
            return make_result();
        }
    }
    return std::nullopt;
}

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t acmacs::seqdb::v3::scan::fasta::scan_ncbi(const std::string_view directory, const scan_options_t& options)
{
    const auto filename = fmt::format("{}/influenza_na.dat.xz", directory);
    const std::string influenza_na_dat = acmacs::file::read(filename);
    AD_DEBUG("influenza_na_dat: {}", influenza_na_dat.size());

    scan_results_t results;
    auto cur = std::begin(influenza_na_dat);
    const auto end = std::end(influenza_na_dat);
    for (size_t line_no = 1; cur != end; ++line_no) {
        if (auto scan_result = influenza_na_read_entry(cur, end, filename, line_no); scan_result.has_value()) {
            auto messages = normalize_name(*scan_result, options.dbg, scan_name_adjustments::ncbi);
            // fmt::print("{:4d} {:8s} \"{}\" {} {}\n", line_no, *res->fasta.type_subtype, res->fasta.name, res->fasta.country, res->sequence.sample_id_by_sample_provider());
            results.results.push_back(std::move(*scan_result));
            std::move(std::begin(messages), std::end(messages), std::back_inserter(results.messages));
        }
    }

    // parse name
    // find sequence

    return results;

} // acmacs::seqdb::v3::scan::fasta::scan_ncbi

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
