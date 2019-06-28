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

// seqdb::v3::Seqdb::Seqdb(std::string&& source)
//     : json_text_(std::move(source))
// {
//     parse(json_text_, entries_);

// } // seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------

seqdb::v3::Seqdb::refs_t seqdb::v3::Seqdb::select_by_name(std::string_view name) const
{
    refs_t refs;
    if (const auto found = std::lower_bound(std::begin(entries_), std::end(entries_), name, [](const auto& entry, std::string_view nam) { return entry.name < nam; }); found != std::end(entries_) && found->name == name) {
        for (size_t seq_no = 0; seq_no < found->seqs.size(); ++seq_no)
            refs.emplace_back(&*found, seq_no);
    }
    return refs;

} // seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
