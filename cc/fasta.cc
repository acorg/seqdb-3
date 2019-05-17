#include <iostream>
#include <algorithm>
#include <cctype>
#include <regex>

#include "acmacs-base/string-split.hh"
#include "acmacs-virus/virus-name.hh"
#include "seqdb-3/fasta.hh"

static Date parse_date(std::string_view source, std::string_view filename, size_t line_no);

// ----------------------------------------------------------------------

std::tuple<acmacs::seqdb::v3::fasta::scan_input_t, acmacs::seqdb::v3::fasta::scan_output_t> acmacs::seqdb::v3::fasta::scan(acmacs::seqdb::v3::fasta::scan_input_t input)
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

} // acmacs::seqdb::v3::fasta::scan

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::fasta::sequence_t> acmacs::seqdb::v3::fasta::name_gisaid_spaces(std::string_view name, std::string_view filename, size_t line_no)
{
    // name | date | passage | lab_id | lab | subtype | lineage

    const auto fields = acmacs::string::split(name, " | ");
    if (fields.size() < 2)
        return std::nullopt;

    acmacs::seqdb::v3::fasta::sequence_t result;
    result.raw_name = fields[0];
    result.date = parse_date(::string::upper(::string::strip(fields[1])), filename, line_no);
    if (fields.size() > 2)
        result.passage = acmacs::virus::Passage{::string::upper(::string::strip(fields[2]))};
    if (fields.size() > 3)
        result.lab_id = ::string::upper(::string::strip(fields[3]));
    if (fields.size() > 4)
        result.lab = ::string::upper(::string::strip(fields[4]));
    if (fields.size() > 5)
        result.virus_type = ::string::upper(::string::strip(fields[5]));
    if (fields.size() > 6)
        result.lineage = ::string::upper(::string::strip(fields[6]));
    result.fasta_name = name;
    return std::move(result);

} // acmacs::seqdb::v3::fasta::name_gisaid_spaces

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::fasta::sequence_t> acmacs::seqdb::v3::fasta::name_gisaid_underscores(std::string_view name, std::string_view filename, size_t line_no)
{
    const auto fields = acmacs::string::split(name, "_|_");
    if (fields.size() < 2)
        return std::nullopt;
    std::string source_without_underscores(name);
    ::string::replace(source_without_underscores, '_', ' ');
    return name_gisaid_spaces(source_without_underscores, filename, line_no);

} // acmacs::seqdb::v3::fasta::name_gisaid_underscores

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::fasta::sequence_t> acmacs::seqdb::v3::fasta::name_plain(std::string_view name, std::string_view /*filename*/, size_t /*line_no*/)
{
    acmacs::seqdb::v3::fasta::sequence_t result;
    result.raw_name = name;
    result.fasta_name = result.raw_name;
    return std::move(result);

} // acmacs::seqdb::v3::fasta::name_plain

// ----------------------------------------------------------------------

acmacs::seqdb::fasta::sequence_t& acmacs::seqdb::v3::fasta::normalize_name(acmacs::seqdb::v3::fasta::sequence_t& source, std::string_view filename, size_t line_no)
{
    // std::cout << source.name << '\n';
    // return source;

    acmacs::virus::Passage no_passage;
    std::tie(source.name, source.reassortant, no_passage, source.annotations) = acmacs::virus::parse_name(source.raw_name);

    // adjust subtype
    // parse passage
    // parse lineage

    if (!no_passage.empty())
        std::clog << "WARNING: " << filename << ':' << line_no << ": name field contains passage: \"" << no_passage << "\" name: \"" << source.fasta_name << '"' << std::endl;
    if (!source.annotations.empty()) {
#include "acmacs-base/global-constructors-push.hh"
        static const std::regex re_valid_annotations{"^\\(([\\d\\-]+|VS\\d+)\\)"}; // Crick stuff from gisaid and HI
#include "acmacs-base/diagnostics-pop.hh"
        if (!std::regex_match(source.annotations, re_valid_annotations))
            std::clog << "WARNING: " << filename << ':' << line_no << ": name contains annotations: \"" << source.annotations << "\" name: \"" << source.fasta_name << '"' << std::endl;
    }
    return source;

} // acmacs::seqdb::v3::fasta::normalize_name

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::fasta::normalize_sequence(std::string_view raw_sequence, std::string_view /*filename*/, size_t /*line_no*/)
{
    std::string sequence(raw_sequence);
    sequence.erase(std::remove_if(std::begin(sequence), std::end(sequence), [](char c) { return c == '\n' || c == '\r'; }), sequence.end());
    std::transform(std::begin(sequence), std::end(sequence), std::begin(sequence), [](char c) { return std::toupper(c); });
    return sequence;

} // acmacs::seqdb::v3::fasta::normalize_sequence

// ----------------------------------------------------------------------

Date parse_date(std::string_view source, std::string_view filename, size_t line_no)
{
    Date result;

    const auto month_and_day_unknown = [&]() -> bool {
        if (source.substr(4) == " (MONTH AND DAY UNKNOWN)") {
            return result.from_string(::string::concat(source.substr(0, 4), "-01-01"), acmacs::throw_on_error::no);
        }
        else
            return false;
    };

    const auto day_unknown = [&]() -> bool {
        if (source.substr(7) == " (DAY UNKNOWN)") {
            return result.from_string(::string::concat(source.substr(0, 7), "-01"), acmacs::throw_on_error::no);
        }
        else
            return false;
    };

    if (!result.from_string(source, acmacs::throw_on_error::no) && !month_and_day_unknown() && !day_unknown())
        std::cerr << "ERROR: " << filename << ':' << line_no << ": cannot parse date: [" << source << ']' << '\n';
    return result;

} // acmacs::seqdb::v3::FastaEntry::parse_date


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
