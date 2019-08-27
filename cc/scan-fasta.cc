#include <map>
#include <regex>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/algorithm.hh"

#include "acmacs-base/string-split.hh"
#include "acmacs-base/read-file.hh"
#include "locationdb/locdb.hh"
#include "acmacs-virus/virus-name.hh"
#include "seqdb-3/scan-fasta.hh"
#include "seqdb-3/scan-align.hh"
#include "seqdb-3/hamming-distance.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb
{
    inline namespace v3
    {
        namespace scan
        {
            namespace fasta
            {
                static date::year_month_day parse_date(const acmacs::uppercase& source, std::string_view filename, size_t line_no);
                static std::string_view parse_lab(const acmacs::uppercase& source, std::string_view filename, size_t line_no);
                static acmacs::virus::type_subtype_t parse_subtype(const acmacs::uppercase& source, std::string_view filename, size_t line_no);
                static std::string_view parse_lineage(const acmacs::uppercase& source, std::string_view filename, size_t line_no);
                static acmacs::seqdb::v3::scan::fasta::hint_t find_hints(std::string_view filename);
                static acmacs::uppercase fix_passage(const acmacs::uppercase& passage);

            } // namespace fasta
        }     // namespace scan
    }         // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::scan::fasta::scan_result_t> acmacs::seqdb::v3::scan::fasta::scan(const std::vector<std::string_view>& filenames, const scan_options_t& options)
{
    using namespace fmt::literals;

    get_locdb(); // load locbd outside of threading code, it is not thread safe

    std::vector<std::vector<scan_result_t>> sequences_per_file(filenames.size());
#pragma omp parallel for default(shared) schedule(static, 4)
    for (size_t f_no = 0; f_no < filenames.size(); ++f_no) {
        const auto& filename = filenames[f_no];
        const auto hints = find_hints(filename);
        try {
            const std::string file_data_s = acmacs::file::read(filename);
            const std::string_view file_data = file_data_s;
            scan_input_t file_input{std::begin(file_data), std::end(file_data)};
            while (!file_input.done()) {
                scan_output_t sequence_ref;
                std::tie(file_input, sequence_ref) = scan(file_input);

                std::optional<scan_result_t> scan_result;
                for (auto parser : {&name_gisaid_fields, &name_gisaid_spaces, &name_gisaid_underscores, &name_plain}) {
                    scan_result = (*parser)(sequence_ref.name, hints, filename, file_input.name_line_no);
                    if (scan_result.has_value())
                        break;
                }
                if (scan_result.has_value()) {
                    auto messages = normalize_name(*scan_result, options.dbg);
                    if (import_sequence(sequence_ref.sequence, scan_result->sequence, options)) {
                        if (!scan_result->sequence.reassortant().empty()  // dates for reassortants in gisaid are irrelevant
                            || scan_result->sequence.lab_in({"NIBSC"})) { // dates provided by NIBSC cannot be trusted, they seem to be put date when they made reassortant
                            scan_result->sequence.remove_dates();
                        }
                        sequences_per_file[f_no].push_back(std::move(*scan_result));
                    }
                }
                else
                    fmt::print(stderr, "WARNING: {}:{}: unable to parse fasta name: {}\n", filename, file_input.name_line_no, sequence_ref.name);
            }
        }
        catch (std::exception& err) {
            throw scan_error(fmt::format("{}: {}", filename, err));
        }
    }

    std::vector<scan_result_t> all_sequences;
    for (auto& per_file : sequences_per_file)
        std::move(per_file.begin(), per_file.end(), std::back_inserter(all_sequences));

    return all_sequences;

} // acmacs::seqdb::v3::scan::fasta::scan

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::fasta::merge_duplicates(std::vector<fasta::scan_result_t>& sequences)
{
    sort_by_name(sequences);

    const auto merge = [](auto first, auto last) {
        if ((last - first) > 1 && !first->sequence.nuc().empty()) {
            std::vector<decltype(first)> entries(static_cast<size_t>(last - first));
            std::generate(std::begin(entries), std::end(entries), [it = first]() mutable { return it++; });
            std::sort(std::begin(entries), std::end(entries), [](const auto& it1, const auto& it2) { return it1->sequence.nuc() < it2->sequence.nuc(); });
            auto master = entries.front();
            for (auto it = std::next(entries.begin()); it != entries.end(); ++it) {
                if (master->sequence.nuc() == (*it)->sequence.nuc()) {
                    master->sequence.merge_from((*it)->sequence);
                    (*it)->remove = true;
                }
                else
                    master = *it;
            }
        }
    };

    auto start = sequences.begin();
    for (auto it = std::next(start); it != sequences.end(); ++it) {
        if (designation(it->sequence) != designation(start->sequence)) {
            merge(start, it);
            start = it;
        }
    }
    merge(start, sequences.end());

    sequences.erase(std::remove_if(std::begin(sequences), std::end(sequences), [](const auto& sr) { return sr.remove; }), std::end(sequences));

} // acmacs::seqdb::v3::scan::fasta::merge_duplicates

// ----------------------------------------------------------------------

std::tuple<acmacs::seqdb::v3::scan::fasta::scan_input_t, acmacs::seqdb::v3::scan::fasta::scan_output_t> acmacs::seqdb::v3::scan::fasta::scan(scan_input_t input)
{
    for (; !input.done() && (*input.first == '\r' || *input.first == '\n'); ++input.first) {
        if (*input.first == '\n')
            ++input.line_no;
    }
    if (input.done())
        return {input, {}};

    if (*input.first != '>')
        throw scan_error(::string::concat(':', input.line_no, ": '>' expected"));
    const auto name_start = ++input.first;
    for (; !input.done() && *input.first != '\n'; ++input.first);
    if (input.done())
        throw scan_error(::string::concat(':', input.line_no, ": unexpected end of input"));
    input.name_line_no = input.line_no;
    ++input.line_no;
    auto name_size = static_cast<size_t>(input.first - name_start);
    if (*(input.first - 1) == '\r')
        --name_size;
    const std::string_view name(name_start, name_size);
    const auto seq_start = ++input.first;

    bool eol = false;
    for (; !input.done(); ++input.first) {
        switch (*input.first) {
          case '>':
              if (eol)
                  return {input, {name, std::string_view(seq_start, static_cast<size_t>(input.first - seq_start))}};
              else
                  throw scan_error(::string::concat(':', input.line_no, ": unexpected '>'"));
              // break;
          case '\n':
              ++input.line_no;
              eol = true;
              break;
          default:
              eol = false;
              break;
        }
    }
    return {input, {name, std::string_view(seq_start, static_cast<size_t>(input.first - seq_start))}};

} // acmacs::seqdb::v3::scan::fasta::scan

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> acmacs::seqdb::v3::scan::fasta::name_gisaid_fields(std::string_view name, const hint_t& /*hints*/, std::string_view filename, size_t line_no)
{
    // Isolate name|a=Isolate ID|b=Type|c=Passage details/history|d=Lineage|e=Collection date|f=Submitter|g=Sample ID by sample provider|
    // h=Sample ID by submitting lab|i=Last modified|j=Originating lab|k=Submitting lab|l=Segment|m=Segment number|n=Identifier|o=DNA Accession no.|p=DNA INSDC|

    auto fields = acmacs::string::split(name, "_|_");
    if (fields.size() != 18 || fields[1].substr(0, 2) != "a=" || !fields.back().empty()) {
        fmt::print(stderr, "DEBUG: not name_gisaid_fields: {} [{}] [{}]\n", fields.size(), fields[1].substr(0, 2), fields.back());
        return std::nullopt;
    }

    scan_result_t result;
    result.fasta.entry_name = name;
    result.fasta.name = fields[0];
    result.fasta.filename = filename;
    result.fasta.line_no = line_no;

    std::string_view lab, lab_id;

    for (auto it = std::next(std::begin(fields)); it != std::prev(std::end(fields)); ++it) {
        if (it->at(1) != '=')
            throw scan_error(fmt::format("line:{} field:{}: unrecognized", line_no, it - std::begin(fields)));
        switch(it->at(0)) {
          case 'a':
              result.sequence.add_isolate_id(it->substr(2));
              break;
          case 'b':
              result.fasta.type_subtype = parse_subtype(it->substr(2), filename, line_no);
              break;
          case 'c':
              result.fasta.passage = it->substr(2);
              break;
          case 'd':
              result.fasta.lineage = acmacs::virus::lineage_t{parse_lineage(it->substr(2), filename, line_no)};
              break;
          case 'e':
              result.sequence.add_date(seqdb::scan::format_date(parse_date(it->substr(2), filename, line_no)));
              break;
          case 'f':
              result.sequence.add_submitter(it->substr(2));
              break;
          case 'g':
              result.sequence.add_sample_id_by_sample_provider(it->substr(2));
              break;
          case 'h':
              lab_id = it->substr(2);
              break;
          case 'i':
              result.sequence.add_gisaid_last_modified(seqdb::scan::format_date(parse_date(it->substr(2), filename, line_no)));
              break;
          case 'j':
              result.sequence.add_originating_lab(it->substr(2));
              break;
          case 'k':
              lab = it->substr(2);
              break;
          case 'l':
              result.sequence.add_gisaid_segment(it->substr(2));
              break;
          case 'm':
              result.sequence.add_gisaid_segment_number(it->substr(2));
              break;
          case 'n':
              result.sequence.add_gisaid_identifier(it->substr(2));
              break;
          case 'o':
              result.sequence.add_gisaid_dna_accession_no(it->substr(2));
              break;
          case 'p':
              result.sequence.add_gisaid_dna_insdc(it->substr(2));
              break;
          default:
            throw scan_error(fmt::format("line:{} field:{}: unrecognized", line_no, it - std::begin(fields)));
        }
    }

    result.sequence.add_lab_id(parse_lab(lab, filename, line_no), lab_id);

    return result;

} // acmacs::seqdb::v3::scan::fasta::name_gisaid_fields

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> acmacs::seqdb::v3::scan::fasta::name_gisaid_spaces(std::string_view name, const hint_t& /*hints*/, std::string_view filename, size_t line_no)
{
    // name | date | passage | lab_id | lab | subtype | lineage

    auto fields = acmacs::string::split(name, " | ");
    if (fields.size() < 2)
        return std::nullopt;
    if (auto& last_field = fields.back(); last_field.back() == '|') {
        last_field.remove_suffix(1);
        while (!last_field.empty() && last_field.back() == ' ')
            last_field.remove_suffix(1);
    }

    scan_result_t result;
    result.fasta.entry_name = name;
    result.fasta.name = fields[0];
    result.fasta.filename = filename;
    result.fasta.line_no = line_no;
    result.sequence.add_date(seqdb::scan::format_date(parse_date(::string::strip(fields[1]), filename, line_no)));
    if (fields.size() > 2)
        result.fasta.passage = ::string::strip(fields[2]);
    if (fields.size() > 4)
        result.sequence.add_lab_id(parse_lab(::string::strip(fields[4]), filename, line_no), ::string::strip(fields[3]));
    if (fields.size() > 5)
        result.fasta.type_subtype = parse_subtype(::string::strip(fields[5]), filename, line_no);
    if (fields.size() > 6)
        result.fasta.lineage = acmacs::virus::lineage_t{parse_lineage(::string::strip(fields[6]), filename, line_no)};

    if (!result.fasta.lineage.empty() && result.fasta.lineage != acmacs::virus::lineage_t{"UNKNOWN"})
        result.sequence.lineage(result.fasta.lineage);

    return result;

} // acmacs::seqdb::v3::scan::fasta::name_gisaid_spaces

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> acmacs::seqdb::v3::scan::fasta::name_gisaid_underscores(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no)
{
    const auto fields = acmacs::string::split(name, "_|_");
    if (fields.size() < 2)
        return std::nullopt;
    std::string source_without_underscores(name);
    ::string::replace(source_without_underscores, '_', ' ');
    return name_gisaid_spaces(source_without_underscores, hints, filename, line_no);

} // acmacs::seqdb::v3::scan::fasta::name_gisaid_underscores

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::scan::fasta::scan_result_t> acmacs::seqdb::v3::scan::fasta::name_plain(std::string_view name, const hint_t& hints, std::string_view /*filename*/, size_t /*line_no*/)
{
    scan_result_t result;
    result.fasta.entry_name = result.fasta.name = name;
    result.sequence.add_lab_id(hints.lab);
    result.fasta.type_subtype = acmacs::virus::type_subtype_t{hints.subtype};
    result.fasta.lineage = acmacs::virus::lineage_t{hints.lineage};
    return result;

} // acmacs::seqdb::v3::scan::fasta::name_plain

// ----------------------------------------------------------------------

#include "acmacs-base/global-constructors-push.hh"
static const std::regex re_valid_annotations{
    "^("
        "\\((?:[\\d\\-ABC]+"
            "|VS\\d+"
            "|SU\\d+"
            "|\\d\\d/\\d\\d\\d"
            "|CNIC-\\w+"
            "|TR-\\d+"
            ")\\)"
        "|[BCD]-?\\d\\.\\d"
        "|CDC\\d+A"
        ")"
}; // Crick stuff from gisaid and HI, C1.4, CDC19A, NIBSC

static const std::regex re_empty_annotations_if_just{"^[\\(\\)_\\-\\s,\\.]+$"};

static const std::regex re_name_ends_with_year{"/(19\\d\\d|20[0-2]\\d)$"};

#include "acmacs-base/diagnostics-pop.hh"

acmacs::seqdb::v3::scan::fasta::messages_t acmacs::seqdb::v3::scan::fasta::normalize_name(acmacs::seqdb::v3::scan::fasta::scan_result_t& source, debug dbg)
{
    if (dbg == debug::yes)
        fmt::print(stderr, "DEBUG: source.fasta.name: {}\n", source.fasta.name);

    fix_gisaid_name(source);

    auto result = acmacs::virus::parse_name(source.fasta.name);
    source.sequence.name(std::move(result.name));
    if (source.sequence.year() >= 2016 && !std::regex_search(*source.sequence.name(), re_name_ends_with_year))
        fmt::print(stderr, "WARNING: no year at the end of name: {} {}:{}\n", source.sequence.name(), source.fasta.filename, source.fasta.line_no);
    // if (auto name_year = acmacs::virus::year(source.sequence.name()); !name_year || (!source.sequence.dates().empty() && *name_year != ::string::from_chars<size_t>(source.sequence.dates().front().substr(0, 4))))
    //     fmt::print(stderr, "WARNING: no year in the name or year in the name does not correspond to the date: {} and {}, fasta name: {}\n", source.sequence.name(), source.sequence.dates(), source.fasta.name);
    if (dbg == debug::yes)
        fmt::print(stderr, "DEBUG: source.sequence.name: {}\n", source.sequence.name());

    // source.sequence.host(std::move(result.host));
    source.sequence.country(std::move(result.country));
    source.sequence.continent(std::move(result.continent));
    source.sequence.reassortant(result.reassortant);
    source.sequence.annotations(std::move(result.extra));

    const auto [passage, passage_extra] = acmacs::virus::parse_passage(fix_passage(source.fasta.passage), acmacs::virus::passage_only::yes);
    if (!passage_extra.empty()) {
        if (passage.empty()) {
            result.messages.emplace_back(acmacs::virus::parse_result_t::message_t::unrecognized_passage, passage_extra);
            source.sequence.add_passage(acmacs::virus::Passage{passage_extra});
        }
        else {
            source.sequence.add_passage(acmacs::virus::Passage{passage});
            source.sequence.annotations(::string::join(" ", {source.sequence.annotations(), passage_extra}));
        }
    }
    else if (!passage.empty())
        source.sequence.add_passage(acmacs::virus::Passage{passage});

    // adjust subtype
    // parse lineage

    // if (!result.passage.empty())
    //     result.messages.emplace_back("name field contains passage", result.passage);

    if (const auto annotations = source.sequence.annotations(); !annotations.empty()) {
        if (std::regex_match(std::begin(annotations), std::end(annotations), re_empty_annotations_if_just))
            source.sequence.remove_annotations();
        else if (!std::regex_match(std::begin(annotations), std::end(annotations), re_valid_annotations))
            result.messages.emplace_back("fasta name contains annotations", annotations);
    }
    return result.messages;

} // acmacs::seqdb::v3::scan::fasta::normalize_name

// ----------------------------------------------------------------------

#include "acmacs-base/global-constructors-push.hh"

static const std::regex re_CSISP_name{"/[\\d_]+(_)(20\\d\\d)\\d\\d\\d\\d$"};
static const std::regex re_year_at_end_of_name{"(19\\d\\d|20[0-2]\\d)$"};
// static const std::regex re_year_3000{"/(30)([0-2]\\d)$"};

#include "acmacs-base/diagnostics-pop.hh"

void acmacs::seqdb::v3::scan::fasta::fix_gisaid_name(scan_result_t& source)
{
    // CSISP has names with the isolation date: A/Valencia/07_0435_20171111 -> A/Valencia/07_0435/2017
    if (std::smatch match_CSISP_name; std::regex_search(source.fasta.name, match_CSISP_name, re_CSISP_name)) {
        // fmt::print("INFO: {}\n", source.fasta.name);
        source.fasta.name = ::string::concat(source.fasta.name.substr(0, static_cast<size_t>(match_CSISP_name.position(1))), '/', match_CSISP_name.str(2));
        // fmt::print("INFO: {}\n", source.fasta.name);
    }
    else if (std::smatch match_year_at_end_of_name; source.fasta.name.size() > 4 && source.fasta.name[source.fasta.name.size() - 5] != '/' && std::regex_search(source.fasta.name, match_year_at_end_of_name, re_year_at_end_of_name)) {
        // A/Iasi/2416022019
        source.fasta.name = ::string::concat(source.fasta.name.substr(0, static_cast<size_t>(match_year_at_end_of_name.position(1))), '/', match_year_at_end_of_name.str(1));
    }
    // else if (std::smatch match_year_3000; std::regex_search(source.fasta.name, match_year_3000, re_year_3000)) {
    //     // A/OMSK/3296/3018
    //     source.fasta.name = ::string::concat(source.fasta.name.substr(0, static_cast<size_t>(match_year_at_end_of_name.position(1))), "/20", match_year_at_end_of_name.str(2));
    // }

} // acmacs::seqdb::v3::scan::fasta::fix_gisaid_name

// ----------------------------------------------------------------------

acmacs::uppercase acmacs::seqdb::v3::scan::fasta::fix_passage(const acmacs::uppercase& passage)
{
    const std::array to_remove{
        std::string_view{"PASSAGE DETAILS:"},
        std::string_view{"PASSAGE HISTORY:"},
        std::string_view{"PASSAGE:"},
        std::string_view{"YAMAGATA LINEAGE;"},
        std::string_view{"YAMAGATA LINEAGE"},
        std::string_view{"VICTORIA LINEAGE;"},
        std::string_view{"VICTORIA LINEAGE"},
        std::string_view{"LINEAGE: SWL;"},
        std::string_view{"LINEAGE: A(H1N1)PDM09"},
        std::string_view{"LINEAGE:"},
    };

    std::string result{*passage};
    for (const auto& en : to_remove) {
        if (const auto found = result.find(en); found != std::string::npos)
            result.erase(found, en.size());
    }
    return ranges::to<std::string>(
        result
        | ranges::views::trim([](char cc) { return std::isspace(cc); }) // remove leading and trailing spaces
        | ranges::views::adjacent_filter([](char first, char second) { return !std::isspace(first) || !std::isspace(second); }) // collapse spaces
                      );

} // acmacs::seqdb::v3::scan::fasta::fix_passage

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::scan::fasta::import_sequence(std::string_view raw_sequence, sequence_t& sequence_data, const scan_options_t& options)
{
    std::string sequence(raw_sequence);
    sequence.erase(std::remove_if(std::begin(sequence), std::end(sequence), [](char c) { return c == '\n' || c == '\r'; }), sequence.end());
    if (sequence.size() < options.remove_too_short_nucs)
        return false;
    sequence_data.import(sequence);
    return true;

} // acmacs::seqdb::v3::scan::fasta::import_sequence

// ----------------------------------------------------------------------

date::year_month_day acmacs::seqdb::v3::scan::fasta::parse_date(const acmacs::uppercase& src, std::string_view filename, size_t line_no)
{
    date::year_month_day result;
    const std::string_view source = src;

    const auto month_and_day_unknown = [source,&result]() -> bool {
        if (source.size() > 25 && source.substr(4) == " (MONTH AND DAY UNKNOWN)") {
            result = date::year_from_string(source.substr(0, 4))/0/0;
            return true;
        }
        else
            return false;
    };

    const auto day_unknown = [source,&result]() -> bool {
        if (source.size() > 15 && source.substr(7) == " (DAY UNKNOWN)") {
            result = date::year_from_string(source.substr(0, 4))/date::month_from_string(source.substr(5, 2))/0;
            return true;
        }
        else
            return false;
    };

    const auto extract_date = [source,&result]() -> bool {
        result = date::from_string(source, date::allow_incomplete::no, date::throw_on_error::no);
        return result.ok();
    };

    if (!source.empty() && !month_and_day_unknown() && !day_unknown() && !extract_date())
        fmt::print(stderr, "ERROR: {}:{}: cannot parse date: [{}]\n", filename, line_no, source);
    return result;

} // acmacs::seqdb::v3::scan::fasta::parse_date

// ----------------------------------------------------------------------

#include "acmacs-base/global-constructors-push.hh"
static const std::map<std::string_view, std::string_view> sLabs{
    {"CENTERS FOR DISEASE CONTROL AND PREVENTION", "CDC"},
    {"CRICK WORLDWIDE INFLUENZA CENTRE", "Crick"},
    {"NATIONAL INSTITUTE FOR MEDICAL RESEARCH", "Crick"},
    {"NATIONAL INSTITUTE OF INFECTIOUS DISEASES (NIID)", "NIID"},
    {"WHO COLLABORATING CENTRE FOR REFERENCE AND RESEARCH ON INFLUENZA", "VIDRL"},
    {"ERASMUS MEDICAL CENTER", "EMC"},
    {"WHO CHINESE NATIONAL INFLUENZA CENTER", "CNIC"},
    {"NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL (NIBSC)", "NIBSC"},
};
#include "acmacs-base/diagnostics-pop.hh"

std::string_view acmacs::seqdb::v3::scan::fasta::parse_lab(const acmacs::uppercase& source, std::string_view /*filename*/, size_t /*line_no*/)
{
    if (const auto found = sLabs.find(source); found != sLabs.end())
        return found->second;
    return source;

} // acmacs::seqdb::v3::scan::fasta::parse_lab

// ----------------------------------------------------------------------

acmacs::virus::type_subtype_t acmacs::seqdb::v3::scan::fasta::parse_subtype(const acmacs::uppercase& source, std::string_view filename, size_t line_no)
{
    if (source.empty())
        fmt::print(stderr, "WARNING: {}:{}: no subtype\n", filename, line_no, source);
    if (source.size() >= 8 && source->front() == 'A') {
        if (source[5] != '0' && source[7] == '0') // H3N0
            return acmacs::virus::type_subtype_t{fmt::format("A({})", source->substr(4, 2))};
        else
            return acmacs::virus::type_subtype_t{fmt::format("A({})", source->substr(4))};
    }
    else if (!source.empty() && source->front() == 'B')
        return acmacs::virus::type_subtype_t{"B"};
    return {};

} // acmacs::seqdb::v3::scan::fasta::parse_subtype

// ----------------------------------------------------------------------

std::string_view acmacs::seqdb::v3::scan::fasta::parse_lineage(const acmacs::uppercase& source, std::string_view /*filename*/, size_t /*line_no*/)
{
    return source;

} // acmacs::seqdb::v3::scan::fasta::parse_lineage

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::hint_t acmacs::seqdb::v3::scan::fasta::find_hints(std::string_view filename)
{
    const auto stem = fs::path{filename}.stem().stem().string();
    const auto fields = acmacs::string::split(stem, "-");
    hint_t hints;
    hints.lab = fields[0];
    if (fields[1] == "h1pdm" || fields[1] == "h1seas" || fields[1] == "h1")
        hints.subtype = "A(H1N1)";
    else if (fields[1] == "h3")
        hints.subtype = "A(H3N2)";
    else if (fields[1] == "b" && fields[0] == "niid") {
        hints.subtype = "B";
        if (fields.size() >= 4) {
            if (fields[3] == "vic")
                hints.lineage = "VICTORIA";
            else if (fields[3] == "yam")
                hints.lineage = "YAMAGATA";
        }
    }
    return hints;

} // acmacs::seqdb::v3::scan::fasta::find_hints

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_false_positive(const std::vector<scan_result_t>& sequences, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::views::filter(is_aligned) | ranges::views::filter(is_different_type_subtype_ignore_h0))
        fmt::format_to(out, "detected:{} | fasta:{} | {} -- {}:{}\n{}\n", sc.sequence.type_subtype(), sc.fasta.type_subtype, sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no, sc.sequence.aa().substr(0, sequence_cutoff));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_false_positive

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_not_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::views::filter([type_subtype_infix](const auto& sc) { return sc.fasta.type_subtype->find(type_subtype_infix) != std::string::npos; }) | ranges::views::filter(isnot_aligned)) {
        // fmt::format_to(out, "{} -- {}:{}\n{}\n", sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no, sc.sequence.aa().substr(0, sequence_cutoff));
        fmt::format_to(out, "{} ::: {} ::: {}:{}\n", sc.sequence.aa().substr(0, sequence_cutoff), sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no);
    }
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_not_aligned

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_aa(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::views::filter([type_subtype_infix](const auto& sc) { return sc.fasta.type_subtype->find(type_subtype_infix) != std::string::npos; }) | ranges::views::filter(is_translated))
        fmt::format_to(out, "{}\n{}\n", sc.fasta.entry_name, sc.sequence.aa().substr(0, sequence_cutoff));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_aa

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_aa_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::views::filter([type_subtype_infix](const auto& sc) { return sc.fasta.type_subtype->find(type_subtype_infix) != std::string::npos; }) | ranges::views::filter(is_aligned)) {
        const auto seq = sc.sequence.aa_aligned();
        fmt::format_to(out, "{} [{}]\n{}\n", sc.sequence.full_name(), seq.size(), seq.substr(0, sequence_cutoff));
    }
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_aa_aligned

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
