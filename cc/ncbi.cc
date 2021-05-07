#include <map>
#include "acmacs-base/regex.hh"
#include "acmacs-base/timeit.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/string-strip.hh"
#include "acmacs-base/string-compare.hh"
#include "acmacs-base/bits.hh"
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
static void merge_dat_fna_names(acmacs::seqdb::v3::scan::fasta::scan_result_t& dat_result, acmacs::messages::messages_t& messages, std::string_view fna_name, const acmacs::messages::position_t& fna_pos);

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t acmacs::seqdb::v3::scan::fasta::scan_ncbi(const std::string_view directory, const scan_options_t& options)
{
    Timeit timeit_scan_ncbi("scan_ncbi: ");

    scan_results_t results = timeit("scan_ncbi (read na.dat)", [&]() { return read_influenza_na_dat(directory, options); });
    timeit("scan_ncbi (read fna)", [&]() { read_influenza_fna(results, directory, options); });

    // remove entries with empty sequences
    results.results.erase(std::remove_if(std::begin(results.results), std::end(results.results), [](const auto& en) { return en.sequence.nuc().empty(); }), std::end(results.results));

    AD_INFO("{} ncbi sequences found in {}", results.results.size(), directory);

    return results;

} // acmacs::seqdb::v3::scan::fasta::scan_ncbi

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::fix_ncbi_name(std::string_view source, acmacs::messages::messages_t& messages, debug /*dbg*/)
{
#include "acmacs-base/global-constructors-push.hh"
    static const std::regex re_prefix_influenza_ab_virus{"^Influenza [AB] virus *", acmacs::regex::icase};
    static const std::regex re_influenza_ab_find_1{".*(?:strain|isolate|H\\d+N\\d+)[\\s:]([AB]/[\\w\\s/\\-\\(\\)]+/\\d+(?:\\s*\\(H\\dN\\d\\))?)", acmacs::regex::icase};
    static const std::regex re_influenza_ab_find_2{"([AB]/[\\w\\s/\\-\\(\\)]+/\\d+(?:\\s*\\(H\\dN\\d\\))?)", acmacs::regex::icase};
    static const std::regex re_find_3{"[\\s\\(]([AB]/[\\w\\s/\\-]+/\\d+(?:\\s*\\(H\\dN\\d\\))?)", acmacs::regex::icase};

    static const std::regex re_influenza_ab_meaningless{"^("
                                                        "H\\d+N\\d+"
                                                        "|"
                                                        "\\w\\w gene for ha?emagglutinin, complete cds"
                                                        "|"
                                                        "ha?emagglutinin (\\([^\\)]+\\) )gene, (complete|partial) cds"
                                                        "|"
                                                        "segment \\d gene for ha?emagglutinin, genomic RNA, strain clone \\w+( \\(H\\d+N\\d+\\))?"
                                                        "|"
                                                        "PX[\\w\\-]+ segment \\d ha?emagglutinin mRNA, (complete|partial) cds"
                                                        ")$",
                                                        acmacs::regex::icase | std::regex::nosubs};
    static const std::regex re_meaningless{"^("
                                           "(Low temperature-adaptable )?equine influenza virus( H\\d+N\\d+)?"
                                           "|"
                                           "Influenza virus type [AB] hemagglutinin gene, \\d'' end"
                                           "|"
                                           "unidentified influenza virus.*"
                                           "|"
                                           "cDNA encoding HA of influenza type [AB]"
                                           "|"
                                           "Sequence \\d+ from Patent \\w+"
                                           "|"
                                           "MULTI PLASMID SYSTEM FOR THE PRODUCTION OF INFLUENZA VIRUS"
                                           "|"
                                           "Recombinant infectious laryngotracheitis virus vaccine"
                                           "|"
                                           "UNVERIFIED.*"
                                           ")$",
                                           acmacs::regex::icase | std::regex::nosubs};
#include "acmacs-base/diagnostics-pop.hh"

    if (std::cmatch match_prefix_influenza_ab_virus; acmacs::regex::search(source, match_prefix_influenza_ab_virus, re_prefix_influenza_ab_virus)) {
        if (auto prefix = acmacs::string::prefix_in_parentheses(source.substr(static_cast<size_t>(match_prefix_influenza_ab_virus.length(0)))); !prefix.empty()) {
            return std::string{acmacs::string::remove_prefix_ignore_case(prefix, "STRAIN ")};
        }
        else if (source.length() == static_cast<size_t>(match_prefix_influenza_ab_virus.length(0)) ||
                 acmacs::regex::search(source, match_prefix_influenza_ab_virus.length(0), re_influenza_ab_meaningless))
            return std::string{}; // no name
        else if (std::cmatch match_influenza_ab_find_1; acmacs::regex::search(source, match_prefix_influenza_ab_virus.length(0), match_influenza_ab_find_1, re_influenza_ab_find_1))
            return match_influenza_ab_find_1.str(1); // fmt::print(rest, "--1 {} -- {}\n", match_influenza_ab_find_1.str(1), source);
        else if (std::cmatch match_influenza_ab_find_2; acmacs::regex::search(source, match_prefix_influenza_ab_virus.length(0), match_influenza_ab_find_2, re_influenza_ab_find_2))
            return match_influenza_ab_find_2.str(1); // fmt::print(rest, "--2 {} -- {}\n", match_influenza_ab_find_2.str(1), source);
        else
            messages.emplace_back(acmacs::messages::key::ncbi_unrecognized, source, MESSAGE_CODE_POSITION);
    }
    else if (std::regex_search(std::begin(source), std::end(source), re_meaningless))
        return std::string{}; // no name
    else if (std::cmatch match_find_3; acmacs::regex::search(source, match_prefix_influenza_ab_virus.length(0), match_find_3, re_find_3))
        return match_find_3.str(1); // fmt::print(rest, "--3 {} -- {}\n", match_find_3.str(1), source);
    else
        messages.emplace_back(acmacs::messages::key::ncbi_unrecognized, source, MESSAGE_CODE_POSITION);
    return std::string{};

} // acmacs::seqdb::v3::scan::fasta::fix_ncbi_name

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
        AD_ERROR("cannot parse date: [{}] (size: {}) @@ {}:{}", source, source.size(), filename, line_no);
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

void merge_dat_fna_names(acmacs::seqdb::v3::scan::fasta::scan_result_t& dat_result, acmacs::messages::messages_t& messages, std::string_view fna_name, const acmacs::messages::position_t& fna_pos)
{
    using namespace acmacs::seqdb::v3::scan::fasta;
    scan_result_t fna_result{dat_result};

    constexpr const std::string_view unknown{"UNKNOWN"};

    const auto add_message = [&]() {
        messages.emplace_back(acmacs::messages::key::ncbi_dat_fna_name_difference, fmt::format("dat:\"{}\" fna:\"{}\"", dat_result.name_fields.full_name(), fna_result.name_fields.full_name()),
                              fna_pos, MESSAGE_CODE_POSITION);
    };

    fna_result.fasta.name = fna_name;
    auto fna_name_messages = normalize_name(fna_result, acmacs::debug::no, scan_name_adjustments::ncbi, print_names::no);

    const auto use_fna = [&]() {
        dat_result.sequence.name(fna_result.sequence.name());
        acmacs::messages::move_and_add_source(messages, std::move(fna_name_messages), fna_pos);
    };

    const auto& fnaf = fna_result.name_fields;
    const auto& datf = dat_result.name_fields;
    enum dat_fna_diff : size_t {
        dat_fna_same = 0,
        subtype_diff = 0b1000000,
        host_diff = 0b100000,
        location_diff = 0b10000,
        isolation_diff = 0b1000,
        year_diff = 0b100,
        reassortant_diff = 0xb10,
        extra_diff = 0b1
    }; // subtype, host, location, isolation, year, reassortant, extra
    const auto dat_fna_diff{acmacs::bits::from_bool(!(datf.subtype == fnaf.subtype), !(datf.host == fnaf.host), !(datf.location == fnaf.location), !(datf.isolation == fnaf.isolation),
                                                    !(datf.year == fnaf.year), !(datf.reassortant == fnaf.reassortant), !(datf.extra == fnaf.extra))};

    // AD_DEBUG("++  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
    switch (acmacs::bits::from_bool(dat_result.name_fields.good(), fna_result.name_fields.good(), dat_result.name_fields.good_but_no_country(), fna_result.name_fields.good_but_no_country())) {
        case 0: // bot are bad
            // use longest?
            if (dat_fna_diff != dat_fna_same)
                add_message(); // AD_DEBUG("BB  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
            break;
        case 0b1100: // both are good
        case 0b1111: // both are not so good (both locations are unknown)
            switch (dat_fna_diff) {
                case dat_fna_same: // use dat, do nothing
                    break;
                case subtype_diff: // use longer subtype
                    if (datf.subtype.size() < fnaf.subtype.size())
                        use_fna();
                    else if (datf.subtype.size() == fnaf.subtype.size())
                        add_message(); // AD_DEBUG("S+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case host_diff: // use longer host in case it is long enough
                    if (datf.host.size() < fnaf.host.size() && fnaf.host.size() > 3)
                        use_fna();
                    else if (datf.host.size() > fnaf.host.size() && datf.host.size() > 3)
                        ; // use dat
                    else
                        add_message(); // AD_DEBUG("H+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case location_diff: // use longer location
                    if (datf.location.size() < fnaf.location.size())
                        use_fna();
                    else if (datf.location.size() == fnaf.location.size())
                        add_message(); // AD_DEBUG("L+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case isolation_diff: // use longer isolation
                    if (datf.isolation.size() < fnaf.isolation.size())
                        use_fna();
                    else if (datf.isolation.size() == fnaf.isolation.size())
                        add_message(); // AD_DEBUG("I+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case year_diff: // check date
                case year_diff | isolation_diff:
                    if (const auto date = dat_result.sequence.date(); date.has_value()) {
                        if (const auto year = date->substr(0, 4); fnaf.year == year)
                            use_fna();
                        else if (datf.year != year)
                            add_message(); // AD_DEBUG("Y+  {:70s}    {:70s} date:{} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), *date, fna_pos.filename,
                                           // fna_pos.line_no);
                    }
                    else
                        add_message(); // AD_DEBUG("Y+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case reassortant_diff:
                    add_message(); // AD_DEBUG("R+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case extra_diff:
                    add_message(); // AD_DEBUG("E+  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
                case subtype_diff | host_diff:
                    switch (acmacs::bits::from_bool(datf.subtype.size() == fnaf.subtype.size(), datf.subtype.size() < fnaf.subtype.size(), datf.host.size() == fnaf.host.size(),
                                                    datf.host.size() < fnaf.host.size())) {
                        case 0b0000: // both in fnaf are longer
                            use_fna();
                            break;
                        case 0b0101: // both in datf are longer
                            break;
                        default:
                            add_message(); // AD_DEBUG("SH+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                    }
                    break;
                case subtype_diff | location_diff:
                    switch (acmacs::bits::from_bool(datf.subtype.size() == fnaf.subtype.size(), datf.subtype.size() < fnaf.subtype.size(), datf.location.size() == fnaf.location.size(),
                                                    datf.location.size() < fnaf.location.size())) {
                        case 0b0000: // both in fnaf are longer
                            use_fna();
                            break;
                        case 0b0101: // both in datf are longer
                            break;
                        default:
                            add_message(); // AD_DEBUG("HL+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                    }
                    break;
                case host_diff | location_diff:
                    switch (acmacs::bits::from_bool(datf.host.size() == fnaf.host.size(), datf.host.size() < fnaf.host.size(), datf.location.size() == fnaf.location.size(),
                                                    datf.location.size() < fnaf.location.size() && fnaf.location != unknown)) {
                        case 0b0000: // both in fnaf are longer and fnaf data is not unknown
                            use_fna();
                            break;
                        case 0b0101: // both in datf are longer
                            if (datf.location == unknown)
                                add_message(); // AD_DEBUG("HL+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                        default:
                            add_message(); // AD_DEBUG("HL+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                    }
                    break;
                case host_diff | isolation_diff:
                    switch (acmacs::bits::from_bool(datf.host.size() == fnaf.host.size(), datf.host.size() < fnaf.host.size(), datf.isolation.size() == fnaf.isolation.size(),
                                                    datf.isolation.size() < fnaf.isolation.size() && fnaf.isolation != unknown)) {
                        case 0b0000: // both in fnaf are longer and fnaf data is not unknown
                            use_fna();
                            break;
                        case 0b0101: // both in datf are longer
                            if (datf.isolation == unknown)
                                add_message(); // AD_DEBUG("HI+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                        default:
                            add_message(); // AD_DEBUG("HI+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                    }
                    break;
                case location_diff | isolation_diff:
                    switch (acmacs::bits::from_bool(datf.location.size() == fnaf.location.size(), datf.location.size() < fnaf.location.size() && fnaf.location != unknown,
                                                    datf.isolation.size() == fnaf.isolation.size(), datf.isolation.size() < fnaf.isolation.size() && fnaf.isolation != unknown)) {
                        case 0b0000: // both in fnaf are longer and fnaf data is not unknown
                            use_fna();
                            break;
                        case 0b0101: // both in datf are longer
                            if (datf.location == unknown || datf.isolation == unknown)
                                add_message(); // AD_DEBUG("LI+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                        default:
                            add_message(); // AD_DEBUG("LI+ {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                            break;
                    }
                    break;
                default:
                    add_message(); // AD_DEBUG("++  {:70s}    {:70s} @@ {}:{}", dat_result.name_fields.full_name(), fna_result.name_fields.full_name(), fna_pos.filename, fna_pos.line_no);
                    break;
            }
            break;
        case 0b1000: // dat is good
        case 0b1101: // dat is better
            break;
        case 0b0100: // fna is good
        case 0b1110: // fna is better
            use_fna();
            break;
    }

} // merge_dat_fna_names

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
        const acmacs::messages::position_t fna_pos{filename_fna, file_input.name_line_no};
        if (const auto fields_fna = acmacs::string::split(sequence_ref.name, "|"); fields_fna.size() == 5) {
            if (const auto found = ncbi_id_to_entry.find(fields_fna[3]); found != ncbi_id_to_entry.end() && import_sequence(sequence_ref.sequence, found->second->sequence, options))
                merge_dat_fna_names(*found->second, results.messages, fields_fna[4], fna_pos);
        }
        else
            results.messages.emplace_back(acmacs::messages::key::ncbi_unrecognized_fna_name, sequence_ref.name, fna_pos, MESSAGE_CODE_POSITION);
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
