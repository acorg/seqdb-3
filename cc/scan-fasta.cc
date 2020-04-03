#include <map>
#include <regex>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/algorithm.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/filesystem.hh"
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

                inline void add_message(messages_t& target, messages_t&& to_add)
                {
                    std::move(std::begin(to_add), std::end(to_add), std::back_inserter(target));
                }

            } // namespace fasta
        }     // namespace scan
    }         // namespace v3
} // namespace acmacs::seqdb

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::fasta::scan_results_t::merge(scan_results_t&& source)
{
    using diff_t = decltype(results.begin() - results.begin());

    const auto pos = static_cast<diff_t>(results.size());
    results.resize(results.size() + source.results.size());
    std::move(std::begin(source.results), std::end(source.results), std::next(std::begin(results), pos));

    const auto pos_messages = static_cast<diff_t>(messages.size());
    messages.resize(messages.size() + source.messages.size());
    std::move(std::begin(source.messages), std::end(source.messages), std::next(std::begin(messages), pos_messages));

} // acmacs::seqdb::v3::scan::fasta::scan_results_t::merge

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t acmacs::seqdb::v3::scan::fasta::scan(const std::vector<std::string_view>& filenames, const scan_options_t& options)
{
    using namespace std::string_view_literals;
    using namespace fmt::literals;

    get_locdb(); // load locbd outside of threading code, it is not thread safe

    std::vector<std::vector<scan_result_t>> sequences_per_file(filenames.size());
    std::vector<messages_t> messages_per_file(filenames.size());
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

                try {
                std::optional<scan_result_t> scan_result;
                for (auto parser : {&name_gisaid_fields, &name_gisaid_spaces, &name_gisaid_underscores, &name_plain}) {
                    scan_result = (*parser)(sequence_ref.name, hints, filename, file_input.name_line_no);
                    if (scan_result.has_value())
                        break;
                }
                if (scan_result.has_value()) {
                    auto messages = normalize_name(*scan_result, options.dbg, options.name_adjustements);
                    if (import_sequence(sequence_ref.sequence, scan_result->sequence, options)) {
                        if (!scan_result->sequence.reassortant().empty()  // dates for reassortants in gisaid are irrelevant
                            || scan_result->sequence.lab_in({"NIBSC"})) { // dates provided by NIBSC cannot be trusted, they seem to be put date when they made reassortant
                            scan_result->sequence.remove_dates();
                        }
                        if (scan_result->fasta.type_subtype.h_or_b() == "B" && scan_result->fasta.lineage.empty())
                            messages.push_back({message_t{"invalid-lineage", fmt::format("no lineage for \"{}\"", scan_result->fasta.name)}, scan_result->fasta.filename, scan_result->fasta.line_no});
                        sequences_per_file[f_no].push_back(std::move(*scan_result));
                        add_message(messages_per_file[f_no], std::move(messages));
                    }
                }
                else
                    fmt::print(stderr, "WARNING: {}:{}: unable to parse fasta name: {}\n", filename, file_input.name_line_no, sequence_ref.name);
                }
                catch (manually_excluded& /*msg*/) {
                    // fmt::print("INFO: manually excluded: {} {} {}:{}\n", sequence_ref.name.substr(0, sequence_ref.name.find("_|_")), msg, filename, file_input.name_line_no);
                }
            }
        }
        catch (scan_error& err) {
            fmt::print(stderr, "{}\n", err);
        }
        catch (std::exception& err) {
            fmt::print(stderr, "{}: error: {}\n", filename, err);
        }
    }

    std::vector<scan_result_t> all_sequences;
    for (auto& per_file : sequences_per_file)
        std::move(per_file.begin(), per_file.end(), std::back_inserter(all_sequences));

    for (auto it = std::next(std::begin(messages_per_file)); it != std::end(messages_per_file); ++it)
        add_message(messages_per_file.front(), std::move(*it));

    return {all_sequences, messages_per_file.front()};

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
    // Isolate name_|_a=Isolate ID_|_b=Type_|_c=Passage details/history_|_d=Lineage_|_e=Collection date_|_f=Submitter_|_g=Sample ID by sample provider_|_
    // h=Sample ID by submitting lab_|_i=Last modified_|_j=Originating lab_|_k=Submitting lab_|_l=Segment_|_m=Segment number_|_n=Identifier_|_o=DNA Accession no._|_p=DNA INSDC_|_
    // x=<exclusion-reason>_|_

    auto fields = acmacs::string::split(name, "_|_");
    if ((fields.size() != 18 && fields.size() != 19) || fields[1].substr(0, 2) != "a=" || !fields.back().empty()) {
        if (fields.size() > 1)
            fmt::print(stderr, "{}:{}: warning: name_gisaid_fields: unexpected number of fields: {}: {}\n", filename, line_no, fields.size(), name);
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
            throw scan_error(fmt::format("{}:{}: error: field {} unrecognized: {}", filename, line_no, it - std::begin(fields), *it));
        if (it->size() > 2) { // non-empty value
            switch (it->at(0)) {
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
                    result.sequence.add_submitter(::string::strip(it->substr(2)));
                    break;
                case 'g':
                    result.sequence.add_sample_id_by_sample_provider(::string::strip(it->substr(2)));
                    break;
                case 'h':
                    lab_id = it->substr(2);
                    break;
                case 'i':
                    result.sequence.add_gisaid_last_modified(seqdb::scan::format_date(parse_date(it->substr(2), filename, line_no)));
                    break;
                case 'j':
                    result.sequence.add_originating_lab(::string::strip(it->substr(2)));
                    break;
                case 'k':
                    lab = it->substr(2);
                    break;
                case 'l':
                    result.sequence.add_gisaid_segment(::string::strip(it->substr(2)));
                    break;
                case 'm':
                    result.sequence.add_gisaid_segment_number(::string::strip(it->substr(2)));
                    break;
                case 'n':
                    result.sequence.add_gisaid_identifier(::string::strip(it->substr(2)));
                    break;
                case 'o':
                    result.sequence.add_gisaid_dna_accession_no(::string::strip(it->substr(2)));
                    break;
                case 'p':
                    result.sequence.add_gisaid_dna_insdc(::string::strip(it->substr(2)));
                    break;
                case 'x': // manually excluded
                    throw manually_excluded{it->substr(2)};
                default:
                    throw scan_error(fmt::format("{}:{}: error: field {} unrecognized: {}", filename, line_no, it - std::begin(fields), *it));
            }
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

acmacs::seqdb::v3::scan::fasta::messages_t acmacs::seqdb::v3::scan::fasta::normalize_name(acmacs::seqdb::v3::scan::fasta::scan_result_t& source, debug dbg, scan_name_adjustments name_adjustements)
{

    switch (name_adjustements) {
      case scan_name_adjustments::gisaid:
          fix_gisaid_name(source, dbg);
          break;
      case scan_name_adjustments::none:
          // if (dbg == debug::yes)
          //     fmt::print(stderr, "DEBUG: source.fasta.name: {}\n", source.fasta.name);
          break;
    }

    auto name_parse_result = acmacs::virus::parse_name(source.fasta.name);
    source.sequence.name(std::move(name_parse_result.name));
    if (source.sequence.year() >= 2016 && !std::regex_search(*source.sequence.name(), re_name_ends_with_year))
        fmt::print(stderr, "{}:{}: warning: no year at the end of name: {}\n", source.fasta.filename, source.fasta.line_no, source.sequence.name());
    // if (auto name_year = acmacs::virus::year(source.sequence.name()); !name_year || (!source.sequence.dates().empty() && *name_year != ::string::from_chars<size_t>(source.sequence.dates().front().substr(0, 4))))
    //     fmt::print(stderr, "WARNING: no year in the name or year in the name does not correspond to the date: {} and {}, fasta name: {}\n", source.sequence.name(), source.sequence.dates(), source.fasta.name);
    // if (dbg == debug::yes)
    //     fmt::print(stderr, "DEBUG: source.sequence.name: {}\n", source.sequence.name());

    // source.sequence.host(std::move(name_parse_result.host));
    source.sequence.country(std::move(name_parse_result.country));
    source.sequence.continent(std::move(name_parse_result.continent));
    source.sequence.reassortant(name_parse_result.reassortant);
    source.sequence.annotations(std::move(name_parse_result.extra));

    messages_t messages;
    for (const auto& msg : name_parse_result.messages)
        messages.push_back({msg, source.fasta.filename, source.fasta.line_no});

    const auto [passage, passage_extra] = acmacs::virus::parse_passage(fix_passage(source.fasta.passage), acmacs::virus::passage_only::yes);
    if (!passage_extra.empty()) {
        if (passage.empty()) {
            messages.push_back({message_t{acmacs::virus::parse_result_t::message_t::unrecognized_passage, passage_extra}, source.fasta.filename, source.fasta.line_no});
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
    //     messages.push_back({message_t{"name field contains passage", result.passage}, source.fasta.filename, source.fasta.line_no});

    if (const auto annotations = source.sequence.annotations(); !annotations.empty()) {
        if (std::regex_match(std::begin(annotations), std::end(annotations), re_empty_annotations_if_just))
            source.sequence.remove_annotations();
        else if (!std::regex_match(std::begin(annotations), std::end(annotations), re_valid_annotations))
            messages.push_back({message_t{"fasta name contains annotations", annotations}, source.fasta.filename, source.fasta.line_no});
    }
    return messages;

} // acmacs::seqdb::v3::scan::fasta::normalize_name

// ----------------------------------------------------------------------

#include "acmacs-base/global-constructors-push.hh"

static const std::regex re_subtype_at_the_end("\\(H[0-9]+(N[0-9]+)?\\)$", std::regex_constants::icase | std::regex_constants::ECMAScript);
static const std::regex re_artefact_at_the_beginning("^[^A-Z]([AB]/)", std::regex_constants::icase | std::regex_constants::ECMAScript);

static const std::regex re_CSISP_name{"/[\\d_]+(_)(20\\d\\d)\\d\\d\\d\\d$"};
static const std::regex re_year_at_end_of_name{"(19\\d\\d|20[0-2]\\d)$"};
// static const std::regex re_year_3000{"/(30)([0-2]\\d)$"};
static const std::regex re_HK_name{"/HK/"};
static const std::regex re_CRIE1_name{"/([0-9]+)/CRIE/"};
static const std::regex re_CRIE2_name{"([^0-9])/CRIE/([0-9]+)/([0-9]+)$"};
static const std::regex re_INCMNSZ_name("/INCMNSZ/([^/]+)/[A-Z][A-Z][A-Z](20[0-9][0-9])/H[0-9]+N[0-9]+", std::regex_constants::icase | std::regex_constants::ECMAScript);
static const std::regex re_CDC_LV_name{"/(19\\d\\d|20[0-2]\\d)[_\\-]?CDC[_\\-]?LV[_\\-]?(\\d+[A-Z]*)$", std::regex_constants::icase | std::regex_constants::ECMAScript};

#include "acmacs-base/diagnostics-pop.hh"

void acmacs::seqdb::v3::scan::fasta::fix_gisaid_name(scan_result_t& source, debug dbg)
{
    const std::string name_orig{dbg == debug::yes ? source.fasta.name : std::string{}};

    if (std::smatch match_subtype_at_the_end; std::regex_search(source.fasta.name, match_subtype_at_the_end, re_subtype_at_the_end))
        source.fasta.name.erase(static_cast<size_t>(match_subtype_at_the_end.position(0)));
    if (std::smatch match_artefact_at_the_beginning; std::regex_search(source.fasta.name, match_artefact_at_the_beginning, re_artefact_at_the_beginning))
        source.fasta.name.erase(0, static_cast<size_t>(match_artefact_at_the_beginning.position(1)));

    // '-' instead of '/'
    if ((source.fasta.name[0] == 'A' || source.fasta.name[0] == 'B') && source.fasta.name[1] == '-' && std::count(std::begin(source.fasta.name), std::end(source.fasta.name), '/') < 2 &&
        std::count(std::begin(source.fasta.name), std::end(source.fasta.name), '-') > 2)
        std::replace(std::begin(source.fasta.name), std::end(source.fasta.name), '-', '/');

    const std::string_view name{source.fasta.name};
    // CSISP has names with the isolation date: A/Valencia/07_0435_20171111 -> A/Valencia/07_0435/2017
    if (std::cmatch match_CSISP_name; std::regex_search(std::begin(name), std::end(name), match_CSISP_name, re_CSISP_name)) {
        // fmt::print("INFO: {}\n", source.fasta.name);
        source.fasta.name = fmt::format("{}/{}", name.substr(0, static_cast<size_t>(match_CSISP_name.position(1))), match_CSISP_name.str(2));
        // fmt::print("INFO: {}\n", source.fasta.name);
    }
    else if (std::cmatch match_CDC_LV_name; name.size() > 4 && std::regex_search(std::begin(name), std::end(name), match_CDC_LV_name, re_CDC_LV_name)) {
        // A/ABU DHABI/240/2018-CDC-LV23A
        source.fasta.name = fmt::format("{}{} CDC-LV{}", name.substr(0, static_cast<size_t>(match_CDC_LV_name.position(1))), match_CDC_LV_name.str(1), match_CDC_LV_name.str(2));
    }
    else if (std::cmatch match_year_at_end_of_name;
             name.size() > 4 && name[source.fasta.name.size() - 5] != '/' && std::regex_search(std::begin(name), std::end(name), match_year_at_end_of_name, re_year_at_end_of_name)) {
        // A/Iasi/2416022019
        source.fasta.name = fmt::format("{}/{}", name.substr(0, static_cast<size_t>(match_year_at_end_of_name.position(1))), match_year_at_end_of_name.str(1));
    }
    // else if (std::cmatch match_year_3000; std::regex_search(std::begin(name), std::end(name), match_year_3000, re_year_3000)) {
    //     // A/OMSK/3296/3018
    //     source.fasta.name = fmt::format("{}/20{}", name.substr(0, static_cast<size_t>(match_year_at_end_of_name.position(1))), match_year_at_end_of_name.str(2));
    // }
    else if (const auto hk_pos = source.fasta.name.find("/HK/"); hk_pos != std::string::npos) {
        source.fasta.name = fmt::format("{}/HONG KONG/{}", name.substr(0, hk_pos), name.substr(hk_pos + 4));
    }
    else if (std::cmatch match_crie1; name.size() > 6 && std::regex_search(std::begin(name), std::end(name), match_crie1, re_CRIE1_name)) {
        // A/Moscow/14/CRIE/2019
        source.fasta.name = fmt::format("{}/CRIE-{}/{}", name.substr(0, static_cast<size_t>(match_crie1.position(0))), match_crie1.str(1),
                                        name.substr(static_cast<size_t>(match_crie1.position(0) + match_crie1.length(0))));
    }
    else if (std::cmatch match_crie2; name.size() > 6 && std::regex_search(std::begin(name), std::end(name), match_crie2, re_CRIE2_name)) {
        // A/Moscow/14/CRIE/2019
        source.fasta.name = fmt::format("{}{}/CRIE-{}/{}", name.substr(0, static_cast<size_t>(match_crie2.position(0))), match_crie2.str(1), match_crie2.str(2), match_crie2.str(3));
    }
    else if (std::cmatch match_incmnsz; name.size() > 20 && std::regex_search(std::begin(name), std::end(name), match_incmnsz, re_INCMNSZ_name)) {
        // A/MEXICO/INCMNSZ/BMR090/JAN2020/H1N1PDM09 -> A/MEXICO/BMR090/2020
        source.fasta.name = fmt::format("{}/{}/{}", name.substr(0, static_cast<size_t>(match_incmnsz.position(0))), match_incmnsz.str(1), match_incmnsz.str(2));
    }
    if (dbg == debug::yes && name_orig != source.fasta.name)
        fmt::print(stderr, "DEBUG: \"{}\" -> \"{}\"\n", name_orig, source.fasta.name);

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
    const auto filter_subtype = [type_subtype=string::split(type_subtype_infix, ",")](const auto& sc) {
        for (const auto& ts : type_subtype) {
            if (ts == "ALL" || sc.fasta.type_subtype->find(ts) != std::string::npos)
                return true;
        }
        return false;
    };

    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::views::filter(filter_subtype) | ranges::views::filter(isnot_aligned)) {
        // fmt::format_to(out, "{} -- {}:{}\n{}\n", sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no, sc.sequence.aa().substr(0, sequence_cutoff));
        fmt::format_to(out, "{}:{}: warning: {} ::: {} \n", sc.fasta.filename, sc.fasta.line_no, sc.sequence.aa().substr(0, sequence_cutoff), sc.fasta.entry_name);
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

acmacs::seqdb::v3::scan::fasta::min_max_dates_t acmacs::seqdb::v3::scan::fasta::min_max_dates(const std::vector<scan_result_t>& sequences)
{
    min_max_dates_t result;
    for (const auto& sc : sequences) {
        const auto isolation_date = sc.sequence.date_simulated();
        if (result.min_isolation_date.empty() || isolation_date < result.min_isolation_date)
            result.min_isolation_date = isolation_date;
        if (result.max_isolation_date.empty() || isolation_date > result.max_isolation_date)
            result.max_isolation_date = isolation_date;
        if (const auto& submission_dates = sc.sequence.gisaid_last_modified(); !submission_dates.empty()) {
            if (result.min_submission_date.empty() || submission_dates.front() < result.min_submission_date)
                result.min_submission_date = submission_dates.front();
            if (result.max_submission_date.empty() || submission_dates.front() > result.max_submission_date)
                result.max_submission_date = submission_dates.front();
        }
    }
    return result;

} // acmacs::seqdb::v3::scan::fasta::min_max_dates

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
