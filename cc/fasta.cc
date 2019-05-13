#include <iostream>
#include <tuple>

#include "acmacs-base/read-file.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::FastaEntry::parse()
{

} // acmacs::seqdb::v3::FastaEntry::parse_raw_name

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::FastaEntry> acmacs::seqdb::v3::fasta_scan(std::string_view filename, std::string_view data)
{
    std::vector<acmacs::seqdb::v3::FastaEntry> result;
    size_t line_no = 1;
    auto name_start = data.end(), sequence_start = data.end();
    std::string_view name_data;
    size_t errors = 0;
    bool newline = true;
    auto cp = data.begin();

    const auto emit_sequence = [&]() {
        if (sequence_start != data.end()) {
            result.emplace_back(name_data, std::string_view(&*sequence_start, static_cast<size_t>(cp - sequence_start)));
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
