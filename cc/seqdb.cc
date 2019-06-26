#include "acmacs-base/read-file.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/seqdb-parse.hh"

// ----------------------------------------------------------------------

seqdb::v3::Seqdb::Seqdb(const std::string& filename)
{
    json_text_ = static_cast<std::string>(acmacs::file::read(filename));
    parse(json_text_, entries_);

} // seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------

seqdb::v3::Seqdb::Seqdb(std::string&& source)
    : json_text_(std::move(source))
{
    parse(json_text_, entries_);

} // seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
