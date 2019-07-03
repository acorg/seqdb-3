#include <map>

// #include <regex>
// #include <array>
// #include <set>

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
                for (auto parser : {&name_gisaid_spaces, &name_gisaid_underscores, &name_plain}) {
                    scan_result = (*parser)(sequence_ref.name, hints, filename, file_input.name_line_no);
                    if (scan_result.has_value())
                        break;
                }
                if (scan_result.has_value()) {
                    auto messages = normalize_name(*scan_result);
                    if (import_sequence(sequence_ref.sequence, scan_result->sequence, options))
                        sequences_per_file[f_no].push_back(std::move(*scan_result));
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

void acmacs::seqdb::v3::scan::fasta::remove_invalid_dates(std::vector<fasta::scan_result_t>& sequences)
{
    for (auto& en : sequences)
        en.sequence.remove_invalid_dates();

} // acmacs::seqdb::v3::scan::fasta::remove_invalid_dates

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
    const std::string_view name(name_start, static_cast<size_t>(input.first - name_start));
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

#include "acmacs-base/diagnostics-pop.hh"

acmacs::seqdb::v3::scan::fasta::messages_t acmacs::seqdb::v3::scan::fasta::normalize_name(acmacs::seqdb::v3::scan::fasta::scan_result_t& source)
{
    // std::cout << source.name << '\n';
    // return source;

    auto result = acmacs::virus::parse_name(source.fasta.name);
    source.sequence.name(std::move(result.name));
    // source.sequence.host(std::move(result.host));
    source.sequence.country(std::move(result.country));
    source.sequence.continent(std::move(result.continent));
    source.sequence.reassortant(result.reassortant);
    source.sequence.annotations(std::move(result.extra));

    const auto [passage, passage_extra] = acmacs::virus::parse_passage(source.fasta.passage, acmacs::virus::passage_only::yes);
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
        result = date::from_string(source, date::throw_on_error::no);
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
    if (source.size() >= 8 && source->front() == 'A')
        return acmacs::virus::type_subtype_t{fmt::format("A({})", source->substr(4))};
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
    for (const auto& sc : sequences | ranges::view::filter(is_aligned) | ranges::view::filter(is_different_type_subtype_ignore_h0))
        fmt::format_to(out, "detected:{} | fasta:{} | {} -- {}:{}\n{}\n", sc.sequence.type_subtype(), sc.fasta.type_subtype, sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no, sc.sequence.aa().substr(0, sequence_cutoff));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_false_positive

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_not_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::view::filter([type_subtype_infix](const auto& sc) { return sc.fasta.type_subtype->find(type_subtype_infix) != std::string::npos; }) | ranges::view::filter(isnot_aligned)) {
        // fmt::format_to(out, "{} -- {}:{}\n{}\n", sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no, sc.sequence.aa().substr(0, sequence_cutoff));
        fmt::format_to(out, "{} ::: {} ::: {}:{}\n", sc.sequence.aa().substr(0, sequence_cutoff), sc.fasta.entry_name, sc.fasta.filename, sc.fasta.line_no);
    }
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_not_aligned

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_aa(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::view::filter([type_subtype_infix](const auto& sc) { return sc.fasta.type_subtype->find(type_subtype_infix) != std::string::npos; }) | ranges::view::filter(is_translated))
        fmt::format_to(out, "{}\n{}\n", sc.fasta.entry_name, sc.sequence.aa().substr(0, sequence_cutoff));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_aa

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::fasta::report_aa_aligned(const std::vector<scan_result_t>& sequences, std::string_view type_subtype_infix, size_t sequence_cutoff)
{
    fmt::memory_buffer out;
    for (const auto& sc : sequences | ranges::view::filter([type_subtype_infix](const auto& sc) { return sc.fasta.type_subtype->find(type_subtype_infix) != std::string::npos; }) | ranges::view::filter(is_aligned)) {
        const auto seq = sc.sequence.aa_aligned();
        fmt::format_to(out, "{} [{}]\n{}\n", sc.sequence.full_name(), seq.size(), seq.substr(0, sequence_cutoff));
    }
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::fasta::report_aa_aligned

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
