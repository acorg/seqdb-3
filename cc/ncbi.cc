#include "acmacs-base/debug.hh"
#include "acmacs-base/read-file.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::fasta::scan_results_t acmacs::seqdb::v3::scan::fasta::scan_ncbi(const std::string_view directory)
{
    const std::string influenza_na_dat = acmacs::file::read(fmt::format("{}/influenza_na.dat.xz", directory));
    AD_DEBUG("influenza_na_dat: {}", influenza_na_dat.size());
    return {};

} // acmacs::seqdb::v3::scan::fasta::scan_ncbi

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
