#include <algorithm>
#include <cctype>
#include <regex>
#include <map>
#include <array>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/read-file.hh"
#include "locationdb/locdb.hh"
#include "acmacs-virus/virus-name.hh"
#include "seqdb-3/fasta.hh"

static Date parse_date(std::string_view source, std::string_view filename, size_t line_no);
static std::string_view parse_lab(std::string_view source, std::string_view filename, size_t line_no);
static std::string parse_subtype(std::string_view source, std::string_view filename, size_t line_no);
static std::string_view parse_lineage(std::string_view source, std::string_view filename, size_t line_no);
static acmacs::seqdb::v3::fasta::hint_t find_hints(std::string_view filename);

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::fasta::scan_result_t> acmacs::seqdb::v3::fasta::scan(const std::vector<std::string_view>& filenames, const scan_options_t& options)
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

                std::optional<sequence_t> seq;
                for (auto parser : {&name_gisaid_spaces, &name_gisaid_underscores, &name_plain}) {
                    seq = (*parser)(sequence_ref.name, hints, filename, file_input.name_line_no);
                    if (seq.has_value())
                        break;
                }
                if (seq.has_value()) {
                    auto messages = normalize_name(*seq);
                    seq->sequence = normalize_sequence(sequence_ref.sequence, messages, options);
                    if (!seq->sequence.empty())
                        sequences_per_file[f_no].push_back({*seq, messages, std::string(filename), file_input.name_line_no});
                }
                else
                    fmt::print(stderr, "WARNING: {}:{}: unable to parse fasta name: {}\n", filename, file_input.name_line_no, sequence_ref.name);
            }
        }
        catch (std::exception& err) {
            throw scan_error("{}: {}"_format(filename, err));
        }
    }

    std::vector<scan_result_t> all_sequences;
    for (auto& per_file : sequences_per_file)
        std::move(per_file.begin(), per_file.end(), std::back_inserter(all_sequences));

    return all_sequences;

} // acmacs::seqdb::v3::fasta::scan

// ----------------------------------------------------------------------

std::tuple<acmacs::seqdb::v3::fasta::scan_input_t, acmacs::seqdb::v3::fasta::scan_output_t> acmacs::seqdb::v3::fasta::scan(scan_input_t input)
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

std::optional<acmacs::seqdb::v3::fasta::sequence_t> acmacs::seqdb::v3::fasta::name_gisaid_spaces(std::string_view name, const hint_t& /*hints*/, std::string_view filename, size_t line_no)
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

    acmacs::seqdb::v3::fasta::sequence_t result;
    result.raw_name = fields[0];
    result.date = parse_date(::string::upper(::string::strip(fields[1])), filename, line_no);
    if (fields.size() > 2)
        result.passage = acmacs::virus::Passage{::string::upper(::string::strip(fields[2]))};
    if (fields.size() > 3)
        result.lab_id = ::string::upper(::string::strip(fields[3]));
    if (fields.size() > 4)
        result.lab = parse_lab(::string::upper(::string::strip(fields[4])), filename, line_no);
    if (fields.size() > 5)
        result.type_subtype = parse_subtype(::string::upper(::string::strip(fields[5])), filename, line_no);
    if (fields.size() > 6)
        result.lineage = parse_lineage(::string::upper(::string::strip(fields[6])), filename, line_no);
    result.fasta_name = name;
    return std::move(result);

} // acmacs::seqdb::v3::fasta::name_gisaid_spaces

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::fasta::sequence_t> acmacs::seqdb::v3::fasta::name_gisaid_underscores(std::string_view name, const hint_t& hints, std::string_view filename, size_t line_no)
{
    const auto fields = acmacs::string::split(name, "_|_");
    if (fields.size() < 2)
        return std::nullopt;
    std::string source_without_underscores(name);
    ::string::replace(source_without_underscores, '_', ' ');
    return name_gisaid_spaces(source_without_underscores, hints, filename, line_no);

} // acmacs::seqdb::v3::fasta::name_gisaid_underscores

// ----------------------------------------------------------------------

std::optional<acmacs::seqdb::v3::fasta::sequence_t> acmacs::seqdb::v3::fasta::name_plain(std::string_view name, const hint_t& hints, std::string_view /*filename*/, size_t /*line_no*/)
{
    acmacs::seqdb::v3::fasta::sequence_t result;
    result.raw_name = name;
    result.fasta_name = result.raw_name;
    result.lab = hints.lab;
    result.type_subtype = hints.subtype;
    result.lineage = hints.lineage;
    return std::move(result);

} // acmacs::seqdb::v3::fasta::name_plain

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

acmacs::seqdb::v3::fasta::messages_t acmacs::seqdb::v3::fasta::normalize_name(acmacs::seqdb::v3::fasta::sequence_t& source)
{
    // std::cout << source.name << '\n';
    // return source;

    auto result = acmacs::virus::parse_name(source.raw_name);
    source.name = result.name;
    source.reassortant = result.reassortant;
    source.annotations = result.extra;

    std::string passage_extra;
    std::tie(source.passage, passage_extra) = acmacs::virus::parse_passage(*source.passage, acmacs::virus::passage_only::yes);
    if (!passage_extra.empty()) {
        if (source.passage.empty()) {
            result.messages.emplace_back(acmacs::virus::parse_result_t::message_t::unrecognized_passage, passage_extra);
            source.passage = acmacs::virus::Passage{passage_extra};
        }
        else
            source.annotations = ::string::join(" ", {source.annotations, passage_extra});
    }

    // adjust subtype
    // parse lineage

    // if (!result.passage.empty())
    //     result.messages.emplace_back("name field contains passage", result.passage);

    if (!source.annotations.empty() && std::regex_match(source.annotations, re_empty_annotations_if_just))
        source.annotations.clear();

    if (!source.annotations.empty()) {
        if (!std::regex_match(source.annotations, re_valid_annotations))
            result.messages.emplace_back("fasta name contains annotations", source.annotations);
    }
    return result.messages;

} // acmacs::seqdb::v3::fasta::normalize_name

// ----------------------------------------------------------------------

static inline std::vector<std::pair<char, size_t>> is_nuc_sequence(std::string_view seq)
{
    const std::string_view nucs{"ACGTUWSMKRYBDHVN-"}; // https://en.wikipedia.org/wiki/Nucleotide
    std::array<size_t, 256> symbols_found;
    std::fill(std::begin(symbols_found), std::end(symbols_found), 0UL);
    bool not_nuc = false;
    for (auto cc : seq) {
        ++symbols_found[static_cast<unsigned char>(cc)];
        if (nucs.find(cc) == std::string_view::npos)
            not_nuc = true;
    }
    std::vector<std::pair<char, size_t>> counts;
    if (not_nuc) {
        for (size_t ord = 0; ord < symbols_found.size(); ++ord)
            if (const auto count = symbols_found[ord]; count)
                counts.emplace_back(ord, count);
        std::sort(std::begin(counts), std::end(counts), [](const auto& e1, const auto& e2) { return e1.second > e2.second; });
    }
    return counts;
}

std::string acmacs::seqdb::v3::fasta::normalize_sequence(std::string_view raw_sequence, messages_t& messages, const scan_options_t& options)
{
    std::string sequence(raw_sequence);
    sequence.erase(std::remove_if(std::begin(sequence), std::end(sequence), [](char c) { return c == '\n' || c == '\r'; }), sequence.end());
    if (sequence.size() < options.remove_too_short_nucs) // return empty if too short
        return {};
    std::transform(std::begin(sequence), std::end(sequence), std::begin(sequence), [](char c) { return std::toupper(c); });
    if (const auto letters = is_nuc_sequence(sequence); !letters.empty())
        messages.emplace_back("not-nucs", acmacs::to_string(letters));
    return sequence;

} // acmacs::seqdb::v3::fasta::normalize_sequence

// ----------------------------------------------------------------------

Date parse_date(std::string_view source, std::string_view filename, size_t line_no)
{
    Date result;

    const auto month_and_day_unknown = [&]() -> bool {
        if (source.size() > 25 && source.substr(4) == " (MONTH AND DAY UNKNOWN)") {
            return result.from_string(::string::concat(source.substr(0, 4), "-01-01"), acmacs::throw_on_error::no);
        }
        else
            return false;
    };

    const auto day_unknown = [&]() -> bool {
        if (source.size() > 15 && source.substr(7) == " (DAY UNKNOWN)") {
            return result.from_string(::string::concat(source.substr(0, 7), "-01"), acmacs::throw_on_error::no);
        }
        else
            return false;
    };

    if (!source.empty() && !result.from_string(source, acmacs::throw_on_error::no) && !month_and_day_unknown() && !day_unknown())
        fmt::print(stderr, "ERROR: {}:{}: cannot parse date: [{}]\n", filename, line_no, source);
    return result;

} // acmacs::seqdb::v3::FastaEntry::parse_date

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

std::string_view parse_lab(std::string_view source, std::string_view /*filename*/, size_t /*line_no*/)
{
    if (const auto found = sLabs.find(source); found != sLabs.end())
        return found->second;
    return source;

} // parse_lab

// ----------------------------------------------------------------------

std::string parse_subtype(std::string_view source, std::string_view filename, size_t line_no)
{
    if (source.empty())
        fmt::print(stderr, "WARNING: {}:{}: no subtype\n", filename, line_no, source);
    if (source.size() >= 8 && source[0] == 'A')
        return fmt::format("A({})", source.substr(4));
    else if (!source.empty() && source[0] == 'B')
        return "B";
    return {};

} // parse_subtype

// ----------------------------------------------------------------------

std::string_view parse_lineage(std::string_view source, std::string_view /*filename*/, size_t /*line_no*/)
{
    return source;

} // parse_lineage

// ----------------------------------------------------------------------

acmacs::seqdb::v3::fasta::hint_t find_hints(std::string_view filename)
{
    const auto stem = fs::path{filename}.stem().stem().string();
    const auto fields = acmacs::string::split(stem, "-");
    acmacs::seqdb::fasta::hint_t hints;
    hints.lab = ::string::upper(fields[0]);
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
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
