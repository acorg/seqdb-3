#include <iostream>
#include <algorithm>
#include <cctype>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/string-split.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::FastaEntry::parse()
{
    if (!name_gisaid_spaces(raw_name) && !name_gisaid_underscores(raw_name)) {
        std::cerr << "ERROR: " << source_ref << ": unrecognized name: " << raw_name << '\n';
        return false;
    }
    return normalize_sequence();

} // acmacs::seqdb::v3::FastaEntry::parse_raw_name

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::FastaEntry::name_gisaid_spaces(std::string_view source)
{
    // name | date | passage | lab_id | lab | subtype |
    const auto fields = acmacs::string::split(source, " | ");
    if (fields.size() < 2)
        return false;

    name = ::string::upper(::string::strip(fields[0]));
    parse_date(::string::upper(::string::strip(fields[1])));

    return true;

} // acmacs::seqdb::v3::FastaEntry::name_gisaid_spaces

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::FastaEntry::name_gisaid_underscores(std::string_view source)
{
    const auto fields = acmacs::string::split(source, "_|_");
    if (fields.size() < 2)
        return false;
    std::string source_without_underscores(source);
    ::string::replace(source_without_underscores, '_', ' ');
    return name_gisaid_spaces(source_without_underscores);

} // acmacs::seqdb::v3::FastaEntry::name_gisaid_underscores

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::FastaEntry::parse_date(std::string_view source)
{
    const auto month_and_day_unknown = [source,this]() -> bool {
        if (source.substr(4) == " (MONTH AND DAY UNKNOWN)") {
            return this->date.from_string(::string::concat(source.substr(0, 4), "-01-01"), acmacs::throw_on_error::no);
        }
        else
            return false;
    };

    const auto day_unknown = [source,this]() -> bool {
        if (source.substr(7) == " (DAY UNKNOWN)") {
            return this->date.from_string(::string::concat(source.substr(0, 7), "-01"), acmacs::throw_on_error::no);
        }
        else
            return false;
    };

    if (!date.from_string(source, acmacs::throw_on_error::no) && !month_and_day_unknown() && !day_unknown())
        std::cerr << "ERROR" << source_ref << ": cannot parse date: [" << source << ']' << '\n';

} // acmacs::seqdb::v3::FastaEntry::parse_date

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::FastaEntry::normalize_sequence()
{
    sequence.erase(std::remove_if(std::begin(sequence), std::end(sequence), [](char c) { return c == '\n' || c == '\r'; }), sequence.end());
    std::transform(std::begin(sequence), std::end(sequence), std::begin(sequence), [](char c) { return std::toupper(c); });
    return true;
    
} // acmacs::seqdb::v3::FastaEntry::normalize_sequence

// ----------------------------------------------------------------------

std::ostream& acmacs::seqdb::v3::operator<<(std::ostream& out, const acmacs::seqdb::v3::FastaEntry& entry)
{
    return out << entry.name << ' ' << entry.date; //  << '\n' << entry.sequence << '\n';
}

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::FastaEntry> acmacs::seqdb::v3::fasta_scan(std::string_view filename, std::string_view data)
{
    std::vector<acmacs::seqdb::v3::FastaEntry> result;
    size_t line_no = 1, name_line_no = 1;
    auto name_start = data.end(), sequence_start = data.end();
    std::string_view name_data;
    size_t errors = 0;
    bool newline = true;
    auto cp = data.begin();

    const auto emit_sequence = [&]() {
        if (sequence_start != data.end()) {
            result.emplace_back(name_data, std::string_view(&*sequence_start, static_cast<size_t>(cp - sequence_start)), source_ref_t{filename, name_line_no});
            name_data = std::string_view{};
        }
        else {
            std::cerr << "ERROR: " << filename << ':' << line_no << ": empty sequence\n";
            ++errors;
        }
    };

    for (; cp != data.end(); ++cp) {
        switch (*cp) {
            case '\n':
                ++line_no;
                newline = true;
                if (name_start != data.end()) {
                    name_data = std::string_view(&*name_start, static_cast<size_t>(cp - name_start));
                    name_start = data.end();
                    sequence_start = cp + 1;
                }
                break;
            case '\r':
                break;
            case '>':
                if (newline) {
                    if (!name_data.empty())
                        emit_sequence();
                    name_start = cp + 1;
                    name_line_no = line_no;
                }
                else {
                    std::cerr << "ERROR: " << filename << ':' << line_no << ": unexpected >\n";
                    ++errors;
                }
                break;
            default:
                newline = false;
                break;
        }
    }
    emit_sequence();

    if (errors)
        throw std::runtime_error("fasta_scan: errors encountered");
    std::for_each(std::begin(result), std::end(result), [](auto& entry) { entry.parse(); });
    return result;

} // acmacs::seqdb::v3::fasta_scan

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::FastaEntry> acmacs::seqdb::v3::fasta_scan(std::string_view filename)
{
    const std::string file_data = acmacs::file::read(filename);
    return fasta_scan(filename, file_data);

} // acmacs::seqdb::v3::fasta_scan

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
