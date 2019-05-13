#include <iostream>

#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::FastaEntry::parse_raw_name()
{

} // acmacs::seqdb::v3::FastaEntry::parse_raw_name

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::FastaEntry> acmacs::seqdb::v3::fasta_scan(std::string_view data)
{
    std::vector<acmacs::seqdb::v3::FastaEntry> result;
    size_t line_no = 1;
    return result;

} // acmacs::seqdb::v3::fasta_scan

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
