#include <map>
#include "acmacs-base/regex.hh"
#include "acmacs-base/timeit.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/string-strip.hh"
#include "acmacs-base/string-compare.hh"
#include "acmacs-virus/defines.hh"
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

static acmacs::virus::type_subtype_t parse_subtype(const acmacs::uppercase& source, acmacs::messages::messages_t& messages, std::string_view filename, size_t line_no);
static std::string fix_country(std::string_view source);
static std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> read_influenza_na_dat_entry(cursor_t& cur, cursor_t end, acmacs::messages::messages_t& messages, std::string_view filename, size_t line_no);
static acmacs::seqdb::v3::scan::fasta::scan_results_t read_influenza_na_dat(const std::string_view directory, const acmacs::seqdb::v3::scan::fasta::scan_options_t& options);
static void read_influenza_fna(acmacs::seqdb::v3::scan::fasta::scan_results_t& results, const std::string_view directory, const acmacs::seqdb::v3::scan::fasta::scan_options_t& options);
static date::year_month_day parse_date(std::string_view source, std::string_view filename, size_t line_no);
static std::string fix_ncbi_name_influenza_a(std::string_view source, acmacs::messages::messages_t& messages, acmacs::debug dbg);
static std::string fix_ncbi_name_influenza_b(std::string_view source, acmacs::messages::messages_t& messages, acmacs::debug dbg);
static std::string fix_ncbi_name_rest(std::string_view source, acmacs::messages::messages_t& messages, acmacs::debug dbg);
static void fix_ncbi_name_remove_meaningless(std::string& source);

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t acmacs::seqdb::v3::scan::fasta::scan_ncbi(const std::string_view directory, const scan_options_t& options)
{
    Timeit timeit("scan_ncbi: ");
    Timeit timeit_na_dat("scan_ncbi (read na.dat): ");

    scan_results_t results = read_influenza_na_dat(directory, options);
    timeit_na_dat.report();
    Timeit timeit_fna("scan_ncbi (read fna): ");
    read_influenza_fna(results, directory, options);
    timeit_fna.report();

    // remove entries with empty sequences
    results.results.erase(std::remove_if(std::begin(results.results), std::end(results.results), [](const auto& en) { return en.sequence.nuc().empty(); }), std::end(results.results));

    AD_INFO("{} ncbi sequences found in {}", results.results.size(), directory);

    return results;

} // acmacs::seqdb::v3::scan::fasta::scan_ncbi

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::fix_ncbi_name(std::string_view source, acmacs::messages::messages_t& messages, debug dbg)
{

    static const std::string_view prefix_a{"INFLUENZA A VIRUS"};
    static const std::string_view prefix_b{"INFLUENZA B VIRUS "};
    static const std::string_view prefix_cdna_a{"CDNA ENCODING HA OF INFLUENZA TYPE A "};

    std::string fixed;
    if (acmacs::string::startswith_ignore_case(source, prefix_a))
        fixed = fix_ncbi_name_influenza_a(source.substr(prefix_a.size()), messages, dbg);
    else if (acmacs::string::startswith_ignore_case(source, prefix_b))
        fixed = fix_ncbi_name_influenza_b(source.substr(prefix_b.size()), messages, dbg);
    else if (acmacs::string::startswith_ignore_case(source, prefix_cdna_a))
        fixed.assign(source.substr(prefix_cdna_a.size()));
    else
        fixed = fix_ncbi_name_rest(source, messages, dbg);
    // fix_ncbi_name_remove_meaningless(fixed);
    ::string::replace_in_place(fixed, '_', ' ');
    acmacs::string::strip_in_place(fixed);
    return fixed;

} // acmacs::seqdb::v3::scan::fasta::fix_ncbi_name

// ----------------------------------------------------------------------

namespace acmacs::string
{
    // if the first symbol of source is '(', returns the segement inside parentheses (including nested parentheses).
    // otherwise returns an empty segment
    // if there is no matching ')' returns the whole source except initial '('
    inline std::string_view prefix_in_parentheses(std::string_view source)
    {
        if (source.empty() || source[0] != '(')
            return {};
        int paren_level{1};
        for (auto cc = std::next(std::begin(source)); cc != std::end(source); ++cc) {
            switch (*cc) {
              case '(':
                  ++paren_level;
                  break;
              case ')':
                  --paren_level;
                  if (paren_level == 0)
                      return source.substr(1, static_cast<size_t>(cc - std::begin(source) - 1));
                  break;
            }
        }
        return source.substr(1);
    }
}

// ----------------------------------------------------------------------

#define NCBI_VIRUS_NAME                                                                                                                                                                                  \
    "("                                                                                                                                                                                                \
    "[^\\(\\)]+"                                                                                                                                                                                       \
    "|"                                                                                                                                                                                                \
    "[AB]/(?:[^/]+/){2,3}\\d{4}"                                                                                                                                                                       \
    ")"

    // "(H1N1)", "(MIXED,H1N1)", "(MIXED.H1N1)", "(MIXED)", "(H1N1)(H1N1)" "(HxNx)"
#define NCBI_SUBTYPE "H\\d{1,2}/?(?:N(?:\\d{1,2}|-|\\?)V?)?\\b"
#define NCBI_SUBTYPE_OPTIONAL "(?:" NCBI_SUBTYPE "\\s+)?"
#define NCBI_MIXED_AND_SUBTYPE "(?:\\((?:MIXED|(?:MIXED[\\.,])?" NCBI_SUBTYPE ")\\))*"
#define NCBI_REASSORTANT_IN_PAREN "(?:\\((X-\\d+[A-Z]*)\\))?"
#define NCBI_PASSAGED_IN_PAREN "(?:\\(([A-Z]+[\\-\\s]PASSAGED?)\\))?"
#define NCBI_CLONE "(?:\\((CLONE)=([A-Z\\d]+)\\))?"

// ----------------------------------------------------------------------

std::string fix_ncbi_name_influenza_a(std::string_view source, acmacs::messages::messages_t& messages, acmacs::debug dbg)
{
    if (source.empty())
        return {};

    if (source[0] == ' ')
        source.remove_prefix(1);

    if (auto prefix = acmacs::string::prefix_in_parentheses(source); !prefix.empty())
        source = prefix;

    using namespace acmacs::regex;
#include "acmacs-base/global-constructors-push.hh"
    static const std::array fix_data_1{
        look_replace_t{std::regex("^" FLU_A_SUBTYPE "\\s(?:STRAIN\\s) A/", std::regex::icase), {"A($1$2$2$4)/$'"}},
    };

    // garbage to remove
    static const std::array fix_data_2{
        look_replace_t{std::regex("^\\(?" FLU_A_SUBTYPE "\\)?$", std::regex::icase), {"$`$'"}},
        look_replace_t{std::regex("(?:"
                                  "ha?emagglutinin (?:(?:precursor *)?HA\\d region (?:\\(HA\\) *)?)?gene, partial cds"
                                  "|"
                                  "segment \\d ha?emagglutinin \\(HA\\) gene, partial cds"
                                  "|"
                                  "HA (?:partial *)?gene for Ha?emagglutinin, (?:genomic RNA, strain|complete cds)"
                                  "|"
                                  "genomic RNA for ha?emagglutinin \\(ha gene\\) strain"
                                  "|"
                                  "partial HA gene for Ha?emagglutinin(?:, genomic RNA| subunit HA1,) strain"
                                  "|"
                                  "strain "
                                  ")", std::regex::icase), {"$`$'"}},
    };
#include "acmacs-base/diagnostics-pop.hh"

    std::string result;
    if (const auto res = scan_replace(source, fix_data_1); res.has_value())
        result.assign(res->front());
    else
        result.assign(source);

    while (true) {
        if (const auto res = scan_replace(result, fix_data_2); res.has_value())
            result.assign(res->front());
        else
            break;
    }

    return result;

//     using namespace acmacs::regex;
// #include "acmacs-base/global-constructors-push.hh"
//     static const std::array fix_data{
//         // allow text at the end, e.g. "segment 4 hemagglutinin (HA) gene, complete cds" found in influenza.fna
//         look_replace_t{std::regex("^ *\\(A/REASSORTANT/(?:A/)?" NCBI_VIRUS_NAME "\\(([^\\(\\)]+) X PUERTO RICO/8/1934\\)" NCBI_MIXED_AND_SUBTYPE "\\)", std::regex::icase), {"A/$2 $1"}},
//         look_replace_t{std::regex("^ *\\(A/REASSORTANT/(?:A/)?" NCBI_VIRUS_NAME "\\(([^\\(\\)]+)\\)" NCBI_MIXED_AND_SUBTYPE "\\)", std::regex::icase), {"A/$2 $1"}},
//         look_replace_t{std::regex("^ *\\(A/((?:NYMC )?X-[\\dA-Z]+)\\((?:A/)?" NCBI_VIRUS_NAME "\\)\\(([^\\(\\)]+)\\)\\)", std::regex::icase), {"A/$2$3 $1"}},
//         look_replace_t{std::regex("^ *\\(" NCBI_VIRUS_NAME NCBI_CLONE NCBI_REASSORTANT_IN_PAREN NCBI_PASSAGED_IN_PAREN NCBI_MIXED_AND_SUBTYPE "\\)", std::regex::icase), {"$1 $2 $3 $4 $5"}},
//         look_replace_t{std::regex("^ *\\(" NCBI_MIXED_AND_SUBTYPE "\\) SEGMENT \\d ISOLATE " NCBI_VIRUS_NAME NCBI_MIXED_AND_SUBTYPE " CRNA SEQUENCE", std::regex::icase), {"$1"}},
//         look_replace_t{std::regex("^ STRAIN " NCBI_VIRUS_NAME NCBI_MIXED_AND_SUBTYPE " SEGMENT ", std::regex::icase), {"$1"}},
//         look_replace_t{std::regex("^[A-Z\\s\\d,]* STRAIN[\\s:]+" NCBI_VIRUS_NAME "\\s*" NCBI_MIXED_AND_SUBTYPE, std::regex::icase), {"$1$2"}},
//         look_replace_t{std::regex("^\\s*" NCBI_SUBTYPE "$", std::regex::icase), {""}},
//         look_replace_t{std::regex("^\\s*" NCBI_SUBTYPE_OPTIONAL NCBI_VIRUS_NAME "\\s*" NCBI_MIXED_AND_SUBTYPE, std::regex::icase), {"$1$2"}}, // " [A-Z]+ GENE "
//     };
// #include "acmacs-base/diagnostics-pop.hh"

//     if (source.empty()) {
//         return {};
//     }
//     else if (const auto res = scan_replace(source, fix_data); res.has_value()) {
//         AD_DEBUG_IF(dbg, "\"{}\" -> \"{}\"", source, res->front());
//         return res->front();
//     }
//     else {
//         messages.emplace_back(acmacs::messages::key::ncbi_influenza_a_not_fixed, source, MESSAGE_CODE_POSITION);
//         return std::string{source};
//     }

}  // fix_ncbi_name_influenza_a

// ----------------------------------------------------------------------

std::string fix_ncbi_name_influenza_b(std::string_view source, acmacs::messages::messages_t& messages, acmacs::debug dbg)
{
    using namespace acmacs::regex;
#include "acmacs-base/global-constructors-push.hh"
    static const std::array fix_data{
        look_replace_t{std::regex("^ *\\(" NCBI_VIRUS_NAME "\\)", std::regex::icase), {"$1"}},
        look_replace_t{std::regex("^ *\\(B/REASSORTANT/(NYMC BX-[\\dA-Z]+)\\((?:LEE/1940|PANAMA/45/1990) X " NCBI_VIRUS_NAME "\\)\\)", std::regex::icase), {"B/$2 $1"}},
        look_replace_t{std::regex("^ *\\(B/REASSORTANT/(NYMC BX-[\\dA-Z]+)\\(" NCBI_VIRUS_NAME "\\)\\)", std::regex::icase), {"B/$2 $1"}},
        look_replace_t{std::regex("^[A-Z\\s\\d,]* STRAIN[\\s:]+" NCBI_VIRUS_NAME, std::regex::icase), {"$1"}},
        look_replace_t{std::regex("\\s+" NCBI_VIRUS_NAME, std::regex::icase), {"$1"}},
    };
#include "acmacs-base/diagnostics-pop.hh"

    if (source.empty()) {
        return {};
    }
    else if (const auto res = scan_replace(source, fix_data); res.has_value()) {
        AD_DEBUG_IF(dbg, "\"{}\" -> \"{}\"", source, res->front());
        return res->front();
    }
    else {
        messages.emplace_back(acmacs::messages::key::ncbi_influenza_b_not_fixed, source, MESSAGE_CODE_POSITION);
        return std::string{source};
    }

} // fix_ncbi_name_influenza_b

// ----------------------------------------------------------------------

std::string fix_ncbi_name_rest(std::string_view source, acmacs::messages::messages_t& messages, acmacs::debug dbg)
{
    using namespace acmacs::regex;

#include "acmacs-base/global-constructors-push.hh"
    static const std::array fix_data{
        look_replace_t{std::regex("^(?:SEQUENCE \\d+ FROM PATENT [^ ]+|unidentified influenza virus.*|(?:Low temperature-adaptable )?Equine influenza virus(?: H3N8)?)$", std::regex::icase), {""}},
        look_replace_t{std::regex("^Influenza\\s+" NCBI_VIRUS_NAME "[A-Z,\\s\\-]*" NCBI_MIXED_AND_SUBTYPE, std::regex::icase), {"$1$2"}},
    };
#include "acmacs-base/diagnostics-pop.hh"

    if (const auto res = scan_replace(source, fix_data); res.has_value()) {
        AD_DEBUG_IF(dbg, "\"{}\" -> \"{}\"", source, res->front());
        return res->front();
    }
    else {
        messages.emplace_back(acmacs::messages::key::ncbi_not_fixed, source, MESSAGE_CODE_POSITION);
        return std::string{source};
    }

} // fix_ncbi_name_rest

// ----------------------------------------------------------------------

void fix_ncbi_name_remove_meaningless(std::string& source)
{
#include "acmacs-base/global-constructors-push.hh"
    static const std::regex re_meaningless("\\s*(?:"
                                           "genomic"
                                           "|"
                                           "gene"
                                           "|"
                                           "RNA"
                                           "|"
                                           "(?:for\\s*)(?:pre)?ha?emagglutinin(?:,?\\s*HA\\d?\\s*DOMAIN|\\s*MRNA)?"
                                           "|"
                                           "(?:precursor\\s*)HA1 region"
                                           "|"
                                           ",?\\s*partial\\s*cds"
                                           "|"
                                           "segment\\s*\\d"
                                           "|"
                                           "\\d\\s*SUBUNIT"
                                           ")",
                                           std::regex::icase);
#include "acmacs-base/diagnostics-pop.hh"

    std::smatch match_meaningless;
    while (std::regex_search(source, match_meaningless, re_meaningless)) {
        source.erase(static_cast<size_t>(match_meaningless.position(0)), static_cast<size_t>(match_meaningless.length(0)));
    }

} // fix_ncbi_name_remove_meaningless

// ----------------------------------------------------------------------

acmacs::virus::type_subtype_t parse_subtype(const acmacs::uppercase& source, acmacs::messages::messages_t& messages, std::string_view filename, size_t line_no)
{
    using namespace acmacs::regex;

#include "acmacs-base/global-constructors-push.hh"

    static const std::array fix_data{
         // allow text at the end, e.g. "segment 4 hemagglutinin (HA) gene, complete cds" found in influenza.fna
        look_replace_t{std::regex("^H\\d{1,2}(?:N\\d{1,2}V?)?(?:NSB)?$", std::regex::icase), {"A($0)"}},
        look_replace_t{std::regex("^(H\\d{1,2})N[X\\-\\?]$", std::regex::icase), {"A($1)"}},
        look_replace_t{std::regex("^(H\\d{1,2})N\\d{1,2}[/,]N?\\d{1,2}$", std::regex::icase), {"A($1)"}},
        look_replace_t{std::regex("^(H\\d{1,2})N\\d{1,2},H\\d{1,2}$", std::regex::icase), {"A"}},
        look_replace_t{std::regex("^(H\\d{1,2})N$", std::regex::icase), {"A($1)"}},
        look_replace_t{std::regex("^H[X\\?I]N[X\\d]$", std::regex::icase), {"A"}},
        look_replace_t{std::regex("^N\\d{1,2}$", std::regex::icase), {"A"}},
        look_replace_t{std::regex("^MIXED[\\.,] *(H\\d{1,2})$", std::regex::icase), {"A($1)"}},
        look_replace_t{std::regex("^MIXED[\\.,] *N\\d{1,2}$", std::regex::icase), {""}},
        look_replace_t{std::regex("^MIXED$", std::regex::icase), {""}},
        look_replace_t{std::regex("^(H\\d{1,2}),MIXED$", std::regex::icase), {"A($1)"}},
        look_replace_t{std::regex("^UNKNOWN$", std::regex::icase), {""}},
    };

#include "acmacs-base/diagnostics-pop.hh"

    if (const auto res = scan_replace(source, fix_data); res.has_value()) {
        return acmacs::virus::type_subtype_t{res->front()};
    }

    messages.emplace_back(acmacs::messages::key::ncbi_unrecognized_subtype, source, acmacs::messages::position_t{filename, line_no}, MESSAGE_CODE_POSITION);
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
        pp{"COTE DIVOIRE"sv,                     "IVORY COAST"sv},
        pp{"COOK ISLANDS"sv,                     "NEW ZEALAND"sv},
        pp{"HONG KONG"sv,                        "CHINA"sv},
        pp{"GREENLAND"sv,                        "DENMARK"sv},
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

std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> read_influenza_na_dat_entry(cursor_t& cur, cursor_t end, acmacs::messages::messages_t& messages, std::string_view filename, size_t line_no)
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
                    result.fasta.type_subtype = parse_subtype(token, messages, filename, line_no);
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
            if (*segment_number != '4') // interested in segment 4 (HA) only
                return std::nullopt;
            if (result.fasta.name.size() > 17 && result.fasta.name[10] == 'C' && ::string::upper(result.fasta.name.substr(0, 17)) == "INFLUENZA C VIRUS")
                return std::nullopt;
            return result;
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
              ok = ::string::upper(source) == "UNKNOWN";
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
        if (auto scan_result = read_influenza_na_dat_entry(cur, end, results.messages, filename_dat, line_no); scan_result.has_value()) {
            auto messages = normalize_name(*scan_result, options.dbg, scan_name_adjustments::ncbi, options.prnt_names);
            // fmt::print("{:4d} {:8s} \"{}\" {} {}\n", line_no, *res->fasta.type_subtype, res->fasta.name, res->fasta.country, res->sequence.sample_id_by_sample_provider());
            if (scan_result->fasta.type_subtype.empty() && !scan_result->sequence.name().empty())
                scan_result->fasta.type_subtype = acmacs::virus::v2::type_subtype_t{std::string(1, scan_result->sequence.name()->front())};

            results.results.push_back(std::move(*scan_result));
            acmacs::messages::move_and_add_source(results.messages, std::move(messages), acmacs::messages::position_t{filename_dat, line_no});
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

    const auto fna_dat_difference_report = [&](std::string_view fna, std::string_view dat, const auto& file_pos) {
        if (fna == dat)
            return;
        if (acmacs::string::startswith(fna, "A(") && acmacs::string::startswith(dat, "A/") && fna.substr(fna.find(')') + 1) == dat.substr(1))
            return;
        results.messages.emplace_back(acmacs::messages::key::ncbi_dat_fna_name_difference, fmt::format("dat:\"{}\" fna:\"{}\"", dat, fna), file_pos, MESSAGE_CODE_POSITION);
    };

    scan_input_t file_input{influenza_fna.begin(), influenza_fna.end()};
    while (!file_input.done()) {
        scan_output_t sequence_ref;
        std::tie(file_input, sequence_ref) = scan(file_input);
        const acmacs::messages::position_t file_pos{filename_fna, file_input.name_line_no};
        if (const auto fields = acmacs::string::split(sequence_ref.name, "|"); fields.size() == 5) {
            if (const auto found = ncbi_id_to_entry.find(fields[3]); found != ncbi_id_to_entry.end()) {
                if (import_sequence(sequence_ref.sequence, found->second->sequence, options)) {
                    // merge names from dat and fna
                    scan_result_t result_for_name_in_fna{*found->second};
                    result_for_name_in_fna.fasta.name = fields[4];
                    acmacs::messages::move_and_add_source(results.messages, normalize_name(result_for_name_in_fna, options.dbg, scan_name_adjustments::ncbi, options.prnt_names),
                                                          acmacs::messages::position_t{filename_fna, file_input.name_line_no});
                    if (!result_for_name_in_fna.sequence.name().empty()) {
                        if (found->second->sequence.name().empty())
                            found->second->sequence.name(result_for_name_in_fna.sequence.name());
                        else
                            fna_dat_difference_report(result_for_name_in_fna.sequence.name(), found->second->sequence.name(), file_pos);
                    }
                }
            }
        }
        else
            results.messages.emplace_back(acmacs::messages::key::ncbi_unrecognized_fna_name, sequence_ref.name, file_pos, MESSAGE_CODE_POSITION);
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
