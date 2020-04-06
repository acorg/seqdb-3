#include <map>
#include "acmacs-base/regex.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/string-split.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

using cursor_t = decltype(std::string{}.cbegin());

inline cursor_t dat_token(cursor_t cursor, cursor_t end)
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

static acmacs::virus::type_subtype_t parse_subtype(const acmacs::uppercase& source, std::string_view filename, size_t line_no);
static std::string fix_country(std::string_view source);
static std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> read_influenza_na_dat_entry(cursor_t& cur, cursor_t end, std::string_view filename, size_t line_no);
static acmacs::seqdb::v3::scan::fasta::scan_results_t read_influenza_na_dat(const std::string_view directory, const acmacs::seqdb::v3::scan::fasta::scan_options_t& options);
static void read_influenza_fna(acmacs::seqdb::v3::scan::fasta::scan_results_t& results, const std::string_view directory, const acmacs::seqdb::v3::scan::fasta::scan_options_t& options);
static date::year_month_day parse_date(std::string_view source, std::string_view filename, size_t line_no);

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t acmacs::seqdb::v3::scan::fasta::scan_ncbi(const std::string_view directory, const scan_options_t& options)
{
    scan_results_t results = read_influenza_na_dat(directory, options);
    read_influenza_fna(results, directory, options);

    // remove entries with empty sequences
    results.results.erase(std::remove_if(std::begin(results.results), std::end(results.results), [](const auto& en) { return en.sequence.nuc().empty(); }), std::end(results.results));

    AD_INFO("{} ncbi sequences found in {}", results.results.size(), directory);

    return results;

} // acmacs::seqdb::v3::scan::fasta::scan_ncbi

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::fix_ncbi_name(std::string_view source, debug dbg)
{
    using namespace acmacs::regex;

#include "acmacs-base/global-constructors-push.hh"

#define RE_VIRUS_NAME "([^\\(\\)]+)"
    // "(H1N1)", "(MIXED,H1N1)", "(MIXED.H1N1)", "(MIXED)", "(H1N1)(H1N1)" "(HxNx)"
#define RE_MIXED_AND_SUBTYPE "(?:\\((?:MIXED|(?:MIXED[\\.,])?(?:[HNX\\d]+))\\))*"
#define RE_CLONE "(?:\\((CLONE)=([A-Z\\d]+)\\))?"

    static const std::array fix_data{
         // allow text at the end, e.g. "segment 4 hemagglutinin (HA) gene, complete cds" found in influenza.fna
        look_replace2_t{std::regex("^INFLUENZA A VIRUS \\(A/REASSORTANT/(?:A/)?" RE_VIRUS_NAME "\\(([^\\(\\)]+) X PUERTO RICO/8/1934\\)" RE_MIXED_AND_SUBTYPE "\\)", std::regex::icase), "A/$2", " $1"},
        look_replace2_t{std::regex("^INFLUENZA A VIRUS \\(A/REASSORTANT/(?:A/)?" RE_VIRUS_NAME "\\(([^\\(\\)]+)\\)" RE_MIXED_AND_SUBTYPE "\\)", std::regex::icase), "A/$2", " $1"},
        look_replace2_t{std::regex("^INFLUENZA A VIRUS \\(" RE_VIRUS_NAME RE_CLONE RE_MIXED_AND_SUBTYPE "\\)", std::regex::icase), "$1", " $2 $3"},
        look_replace2_t{std::regex("^INFLUENZA A VIRUS STRAIN " RE_VIRUS_NAME RE_MIXED_AND_SUBTYPE " SEGMENT ", std::regex::icase), "$1", ""},
        look_replace2_t{std::regex("^INFLUENZA A VIRUS " RE_VIRUS_NAME RE_MIXED_AND_SUBTYPE " [A-Z]+ GENE ", std::regex::icase), "$1", ""},
        look_replace2_t{std::regex("^INFLUENZA A VIRUS \\(" RE_MIXED_AND_SUBTYPE "\\) SEGMENT \\d ISOLATE " RE_VIRUS_NAME RE_MIXED_AND_SUBTYPE " CRNA SEQUENCE", std::regex::icase), "$1", ""},
        // name absent
        look_replace2_t{std::regex("^(?:CDNA ENCODING HA OF INFLUENZA TYPE A|SEQUENCE \\d+ FROM PATENT [^ ]+)$", std::regex::icase), " ", ""},
    };

#include "acmacs-base/diagnostics-pop.hh"

    const auto remove_like_at_end = [](std::string_view text) {
        if (text.size() > 5 && ::string::upper(text.substr(text.size() - 5)) == "-LIKE")
            text.remove_suffix(5);
        return text;
    };

    if (const auto [r1, r2] = scan_replace2(source, fix_data); !r1.empty()) {
        const std::string fixed = ::string::strip(::string::replace(fmt::format("{}{}", remove_like_at_end(r1), r2), '_', ' '));
        AD_DEBUG_IF(dbg, "\"{}\" -> \"{}\"", source, fixed);
        return fixed;
    }
    if (::string::upper(source) == "INFLUENZA A VIRUS")
        return {};

    AD_WARNING("fix_ncbi_name: unable to fix: \"{}\"", source);
    return std::string(source);

} // acmacs::seqdb::v3::scan::fasta::fix_ncbi_name

// ----------------------------------------------------------------------

acmacs::virus::type_subtype_t parse_subtype(const acmacs::uppercase& source, std::string_view filename, size_t line_no)
{
    using namespace acmacs::regex;

#include "acmacs-base/global-constructors-push.hh"

    static const std::array fix_data{
         // allow text at the end, e.g. "segment 4 hemagglutinin (HA) gene, complete cds" found in influenza.fna
        look_replace_t{std::regex("^H\\d{1,2}(?:N\\d)?$", std::regex::icase), "A($0)"},
        look_replace_t{std::regex("^(H\\d{1,2})NX$", std::regex::icase), "A($1)"},
        look_replace_t{std::regex("^(H\\d{1,2})N\\d,\\d", std::regex::icase), "A($1)"},
        look_replace_t{std::regex("^HXNX$", std::regex::icase), "A"},
        look_replace_t{std::regex("^N\\d$", std::regex::icase), "A"},
        look_replace_t{std::regex("^MIXED[\\.,](H\\d{1,2})$", std::regex::icase), "A($1)"},
        look_replace_t{std::regex("^MIXED[\\.,]N\\d$", std::regex::icase), "A"},
        look_replace_t{std::regex("^MIXED$", std::regex::icase), "A"},
    };

#include "acmacs-base/diagnostics-pop.hh"

    if (const auto res = scan_replace(source, fix_data); !res.empty()) {
        if (res == "A")
            return acmacs::virus::type_subtype_t{};
        else
            return acmacs::virus::type_subtype_t{res};
    }

    AD_WARNING("unrecognized subtype: \"{}\" @@ {}:{}", source, filename, line_no);
    return acmacs::virus::type_subtype_t{};
}

// ----------------------------------------------------------------------

std::string fix_country(std::string_view source)
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

std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> read_influenza_na_dat_entry(cursor_t& cur, cursor_t end, std::string_view filename, size_t line_no)
{
    acmacs::seqdb::v3::scan::fasta::scan_result_t result;
    result.fasta.filename = filename;
    result.fasta.line_no = line_no;

    na_field field{na_field::genbank_accession};
    cursor_t segment_number;
    for (auto tok_beg = cur, tok_end = dat_token(tok_beg, end); tok_end != end; tok_beg = std::next(tok_end), tok_end = dat_token(tok_beg, end), ++field) {
        if (tok_beg != tok_end) {
            const auto token = string_view(tok_beg, tok_end);
            switch (field) {
                case na_field::genbank_accession:
                    result.sequence.add_sample_id_by_sample_provider(token);
                    break;
                case na_field::segment_no:
                    segment_number = tok_beg;
                    result.sequence.add_gisaid_segment_number(token);
                    break;
                case na_field::virus_name:
                    result.fasta.name = token;
                    break;
                case na_field::subtype:
                    result.fasta.type_subtype = parse_subtype(token, filename, line_no);
                    break;
                case na_field::date:
                    if (const auto dt = parse_date(token, filename, line_no); date::year_ok(dt))
                        result.sequence.add_date(acmacs::seqdb::scan::format_date(dt));
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

date::year_month_day parse_date(std::string_view source, std::string_view filename, size_t line_no)
{
    date::year_month_day result = date::invalid_date();
    auto ok{false};

    switch (source.size()) {
      case 4:                   // year
          result = date::year_from_string(source)/0/0;
          ok = date::year_ok(result);
          break;
      case 7:                   // year/month
          if (source[4] == '/') {
              result = date::year_from_string(source.substr(0, 4)) / date::month_from_string(source.substr(5)) / 0;
              ok = date::year_ok(result) && date::month_ok(result);
          }
          else
              ok = source == "UNKNOWN";
          break;
      case 9:
          ok = source.substr(0, 4) == "NON/";
          break;
      case 10:                  // year/month/day
          result = date::from_string(source, "%Y/%m/%d");
          ok = result.ok();
          break;
      case 3:
          ok = source == "NON";
          break;
      case 0:
          ok = true;
          break;
    }
    if (!ok)
        AD_ERROR("cannot parse date: [{}] @@ {}:{}", source, filename, line_no);
    return result;

} // parse_date

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t read_influenza_na_dat(const std::string_view directory, const acmacs::seqdb::v3::scan::fasta::scan_options_t& options)
{
    using namespace acmacs::seqdb::v3::scan::fasta;

    scan_results_t results;

    const auto filename_dat = fmt::format("{}/influenza_na.dat.xz", directory);
    const std::string influenza_na_dat = acmacs::file::read(filename_dat);
    // AD_DEBUG("influenza_na_dat: {}", influenza_na_dat.size());

    auto cur = std::begin(influenza_na_dat);
    const auto end = std::end(influenza_na_dat);
    for (size_t line_no = 1; cur != end; ++line_no) {
        if (auto scan_result = read_influenza_na_dat_entry(cur, end, filename_dat, line_no); scan_result.has_value()) {
            auto messages = normalize_name(*scan_result, options.dbg, scan_name_adjustments::ncbi);
            // fmt::print("{:4d} {:8s} \"{}\" {} {}\n", line_no, *res->fasta.type_subtype, res->fasta.name, res->fasta.country, res->sequence.sample_id_by_sample_provider());
            if (scan_result->fasta.type_subtype.empty() && !scan_result->sequence.name().empty())
                scan_result->fasta.type_subtype = acmacs::virus::v2::type_subtype_t{std::string(1, scan_result->sequence.name()->front())};

            results.results.push_back(std::move(*scan_result));
            std::move(std::begin(messages), std::end(messages), std::back_inserter(results.messages));
        }
    }
    AD_INFO("{} HA entries found in \"{}\"", results.results.size(), filename_dat);

    return results;
}

// ----------------------------------------------------------------------

void read_influenza_fna(acmacs::seqdb::v3::scan::fasta::scan_results_t& results, const std::string_view directory, const acmacs::seqdb::v3::scan::fasta::scan_options_t& options)
{
    using namespace acmacs::seqdb::v3::scan::fasta;

    std::map<std::string, scan_result_t*, std::less<>> ncbi_id_to_entry;
    for (auto& en : results.results)
        ncbi_id_to_entry[en.sequence.sample_id_by_sample_provider().front()] = &en;

    const auto filename_fna = fmt::format("{}/influenza.fna.xz", directory);
    const std::string influenza_fna_s = acmacs::file::read(filename_fna);
    const std::string_view influenza_fna(influenza_fna_s);
    // AD_DEBUG("influenza_fna: {}", influenza_fna.size());

    scan_input_t file_input{influenza_fna.begin(), influenza_fna.end()};
    while (!file_input.done()) {
        scan_output_t sequence_ref;
        std::tie(file_input, sequence_ref) = scan(file_input);
        if (const auto fields = acmacs::string::split(sequence_ref.name, "|"); fields.size() == 5) {
            if (const auto found = ncbi_id_to_entry.find(fields[3]); found != ncbi_id_to_entry.end()) {
                if (import_sequence(sequence_ref.sequence, found->second->sequence, options)) {
                    // merge names from dat and fna
                    scan_result_t result_for_name_in_fna{*found->second};
                    result_for_name_in_fna.fasta.name = fields[4];
                    results.merge(normalize_name(result_for_name_in_fna, options.dbg, scan_name_adjustments::ncbi));
                    if (!result_for_name_in_fna.sequence.name().empty()) {
                        if (found->second->sequence.name().empty())
                            found->second->sequence.name(result_for_name_in_fna.sequence.name());
                        else if (result_for_name_in_fna.sequence.name() != found->second->sequence.name())
                            results.messages.push_back({{"ncbi-dat-fna-name-difference", fmt::format("dat:\"{}\" fna:\"{}\"", found->second->sequence.name(), result_for_name_in_fna.sequence.name())}, result_for_name_in_fna.fasta.filename, result_for_name_in_fna.fasta.line_no});
                    }
                }
            }
        }
        else
            AD_WARNING("unrecognized ncbi fna name: \"{}\"", sequence_ref.name);
    }

}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
