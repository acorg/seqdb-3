#include <algorithm>

#include "acmacs-base/read-file.hh"
#include "acmacs-virus/virus-name.hh"
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

seqdb::v3::subset seqdb::v3::Seqdb::all() const
{
    subset ss;
    for (const auto& entry : entries_) {
        for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no)
            ss.refs_.emplace_back(&entry, seq_no);
    }
    return ss;

} // seqdb::v3::Seqdb::all

// ----------------------------------------------------------------------

seqdb::v3::subset seqdb::v3::Seqdb::select_by_name(std::string_view name) const
{
    subset ss;
    if (const auto found = std::lower_bound(std::begin(entries_), std::end(entries_), name, [](const auto& entry, std::string_view nam) { return entry.name < nam; }); found != std::end(entries_) && found->name == name) {
        for (size_t seq_no = 0; seq_no < found->seqs.size(); ++seq_no)
            ss.refs_.emplace_back(&*found, seq_no);
    }
    return ss;

} // seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------

std::string_view seqdb::v3::SeqdbEntry::host() const
{
    if (const auto ho = acmacs::virus::host(acmacs::virus::v2::virus_name_t{name}); !ho.empty())
        return ho;
    else
        return "HUMAN";

} // seqdb::v3::SeqdbEntry::host

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::multiple_dates(bool do_filter)
{
    if (do_filter)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.entry->dates.size() < 2; }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::multiple_dates

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::subtype(const acmacs::uppercase& virus_type)
{
    if (!virus_type.empty()) {
        std::string_view vt = virus_type;
        if (vt == "H1")
            vt = "A(H1N1)";
        else if (vt == "H3")
            vt = "A(H3N2)";
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [vt](const auto& en) { return en.entry->virus_type != vt; }), std::end(refs_));
    }
    return *this;

} // seqdb::v3::subset::subtype

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::lineage(const acmacs::uppercase& lineage)
{
    if (!lineage.empty()) {
        std::string_view lin = lineage;
        switch (lin[0]) {
          case 'V': lin = "VICTORIA"; break;
          case 'Y': lin = "YAMAGATA"; break;
        }
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lin](const auto& en) { return en.entry->lineage != lin; }), std::end(refs_));
    }
    return *this;

} // seqdb::v3::subset::lineage

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::lab(const acmacs::uppercase& lab)
{
    if (!lab.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lab=static_cast<std::string_view>(lab)](const auto& en) { return !en.has_lab(lab); }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::lab

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::host(const acmacs::uppercase& host)
{
    if (!host.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [host=static_cast<std::string_view>(host)](const auto& en) { return en.entry->host() != host; }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::host

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::continent(const acmacs::uppercase& continent)
{
    if (!continent.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [continent=static_cast<std::string_view>(continent)](const auto& en) { return en.entry->continent != continent; }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::continent

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::country(const acmacs::uppercase& country)
{
    if (!country.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [country=static_cast<std::string_view>(country)](const auto& en) { return en.entry->country != country; }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::country

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::clade(const acmacs::uppercase& clade)
{
    if (!clade.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [clade=static_cast<std::string_view>(clade)](const auto& en) { return !en.has_clade(clade); }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::clade

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::recent(size_t recent)
{
    if (recent > 0) {
        sort_by_date_recent_first();
        refs_.erase(std::next(std::begin(refs_), static_cast<ssize_t>(recent)), std::end(refs_));
    }
    return *this;

} // seqdb::v3::subset::recent

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::dates(std::string_view start, std::string_view end)
{
    if (!start.empty() || !end.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [start,end](const auto& en) { return !en.entry->date_within(start, end); }), std::end(refs_));
    return *this;

} // seqdb::v3::subset::dates

// ----------------------------------------------------------------------

seqdb::v3::subset& seqdb::v3::subset::print()
{
    for (const auto& ref : *this)
        fmt::print("{}{}{} {}\n", ref.seq_id(), ref.entry->lineage.empty() ? "" : " L:", ref.entry->lineage.empty() ? std::string_view{} : ref.entry->lineage, ref.entry->dates);
    return *this;

} // seqdb::v3::subset::print

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
