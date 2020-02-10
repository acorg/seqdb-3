#include <algorithm>
#include <random>
#include <regex>
#include <numeric>
#include <memory>
#include <cstdlib>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/in-json-parser.hh"
#include "acmacs-base/to-json.hh"
#include "acmacs-virus/virus-name.hh"
#include "acmacs-chart-2/chart-modify.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/seqdb-parse.hh"
#include "seqdb-3/hamming-distance.hh"

// ----------------------------------------------------------------------

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

static std::string sSeqdbFilename = acmacs::acmacsd_root() + "/data/seqdb.json.xz";

#pragma GCC diagnostic pop

void acmacs::seqdb::v3::setup(std::string_view filename)
{
    if (!filename.empty())
        sSeqdbFilename = filename;

} // acmacs::seqdb::v3::setup

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::Seqdb& acmacs::seqdb::v3::Seqdb::get()
{
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif
    static Seqdb sSeqdb(sSeqdbFilename);
#pragma GCC diagnostic pop

    return sSeqdb;

} // acmacs::seqdb::v3::get

// ----------------------------------------------------------------------

acmacs::seqdb::v3::Seqdb::Seqdb(std::string_view filename)
{
    try {
        json_text_ = static_cast<std::string>(acmacs::file::read(filename));
        parse(json_text_, entries_);
    }
    catch (in_json::error& err) {
        fmt::print(stderr, "{}:{}:{}: error: {}\n", filename, err.line_no, err.column_no, err.message);
        std::exit(99);
    }
    catch (std::exception& err) {
        fmt::print(stderr, "WARNING: seqdb not loaded: {}\n", err);
        json_text_.clear();
        entries_.clear();
    }

} // acmacs::seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::Seqdb::Seqdb(std::string&& source)
//     : json_text_(std::move(source))
// {
//     parse(json_text_, entries_);

// } // acmacs::seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::all() const
{
    subset ss;
    ss.refs_.reserve(entries_.size() * 2);
    for (const auto& entry : entries_) {
        for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no)
            ss.refs_.emplace_back(&entry, seq_no);
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::all

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_seq_id(std::string_view seq_id) const
{
    subset ss;
    if (const auto found = seq_id_index().find(seq_id); found != seq_id_index().end())
        ss.refs_.emplace_back(found->second);
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_seq_id(const std::vector<std::string_view>& seq_ids) const
{
    subset ss;
    for (const auto& seq_id : seq_ids) {
        if (const auto found = seq_id_index().find(seq_id); found != seq_id_index().end())
            ss.refs_.emplace_back(found->second);
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_name(std::string_view name) const
{
    subset ss;
    select_by_name(name, ss);
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_name(const std::vector<std::string_view>& names) const
{
    subset ss;
    for (const auto& name : names)
        select_by_name(name, ss);

    return ss;
} // acmacs::seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Seqdb::select_by_name(std::string_view name, subset& subs) const
{
    const auto find_name = [&subs, this](std::string_view look_for) {
        // fmt::print(stderr, "DEBUG: select_by_name \"{}\"\n", look_for);
        if (const auto found = std::lower_bound(std::begin(entries_), std::end(entries_), look_for, [](const auto& entry, std::string_view nam) { return entry.name < nam; });
            found != std::end(entries_) && found->name == look_for) {
            for (size_t seq_no = 0; seq_no < found->seqs.size(); ++seq_no)
                subs.refs_.emplace_back(&*found, seq_no);
        }
        // else if (found != std::end(entries_))
        //     fmt::print(stderr, "DEBUG: not found \"{}\"\n", found->name);
        // else
        //     fmt::print(stderr, "DEBUG: not found\n", found->name);
    };

    const auto subs_initial_size = subs.size();
    find_name(name);
    if (subs.size() == subs_initial_size) {
        if (name[0] == 'A' || name[0] == 'a' || name[0] == 'B' || name[0] == 'b') {
            const auto result = acmacs::virus::parse_name(name);
            find_name(*result.name);
            if (subs.size() == subs_initial_size && (name[0] == 'A'  || name[0] == 'a') && name[1] == '/') {
                for (const char* subtype : {"A(H1N1)/", "A(H3N2)/"}) {
                    find_name(*acmacs::virus::parse_name(std::string{subtype} + std::string{name.substr(2)}).name);
                }
            }
        }
        else {
            for (const char* subtype : {"A(H1N1)/", "A(H3N2)/", "B/"}) {
                find_name(*acmacs::virus::parse_name(std::string{subtype} + std::string{name}).name);
            }
        }
    }

} // acmacs::seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_regex(std::string_view re) const
{
    std::regex reg(std::begin(re), std::end(re), std::regex_constants::icase);
    subset ss;
    for (const auto& entry : entries_) {
        for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
            if (ref candidate{&entry, seq_no}; std::regex_search(candidate.full_name(), reg))
                ss.refs_.push_back(std::move(candidate));
        }
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_regex

// ----------------------------------------------------------------------

acmacs::seqdb::v3::ref acmacs::seqdb::v3::Seqdb::find_hi_name(std::string_view full_name) const
{
    const auto& hi_name_index = acmacs::seqdb::get().hi_name_index();
    if (const auto indp = hi_name_index.find(full_name); indp != hi_name_index.end())
        return indp->second;
    else
        return {};

} // acmacs::seqdb::v3::Seqdb::find_hi_name

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::seq_id_index_t& acmacs::seqdb::v3::Seqdb::seq_id_index() const
{
    if (seq_id_index_.empty()) {
        seq_id_index_.reserve(entries_.size() * 2);
        for (const auto& entry : entries_) {
            for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
                ref rf{&entry, seq_no};
                seq_id_index_.emplace(rf.seq_id(), std::move(rf));
            }
        }
        seq_id_index_.sort_by_key();
    }
    return seq_id_index_;

} // acmacs::seqdb::v3::Seqdb::seq_id_index

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::hi_name_index_t& acmacs::seqdb::v3::Seqdb::hi_name_index() const
{
    if (hi_name_index_.empty()) {
        hi_name_index_.reserve(entries_.size() * 2);
        for (const auto& entry : entries_) {
            for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
                for (const auto& hi_name : entry.seqs[seq_no].hi_names)
                    hi_name_index_.emplace(hi_name, ref{&entry, seq_no});
            }
        }
        hi_name_index_.sort_by_key();
    }
    return hi_name_index_;

} // acmacs::seqdb::v3::Seqdb::hi_name_index

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::Antigens& aAntigens, std::string_view /*aChartVirusType*/) const
{
    // check lineage?
    // check virus type

    subset result;

    auto find_by_hi_name = [this](const auto& antigen) -> std::optional<ref> {
        const auto& hi_name_ind = hi_name_index();
        if (auto found_ref1 = hi_name_ind.find(antigen.full_name()); found_ref1 != hi_name_ind.end())
            return found_ref1->second;
        else if (auto found_ref2 = hi_name_ind.find(antigen.full_name_for_seqdb_matching()); found_ref2 != hi_name_ind.end())
            return found_ref2->second;
        else
            return std::nullopt;
    };

    size_t matched = 0;
    for (auto antigen : aAntigens) {
        if (auto found_ref = find_by_hi_name(*antigen); found_ref.has_value()) {
            result.refs_.push_back(std::move(*found_ref));
            ++matched;
        }
        else if (antigen->passage().empty()) {
            bool found = false;
            for (const auto& selected : select_by_name(antigen->name())) {
                if (selected.seq().has_reassortant(*antigen->reassortant())) {
                    result.refs_.push_back(selected);
                    ++matched;
                    found = true;
                }
            }
            if (!found)
                result.refs_.emplace_back();
        }
        else
            result.refs_.emplace_back();
    }
    fmt::print("INFO: antigens from chart have sequences in seqdb: {}\n", matched);

    return result;

} // acmacs::seqdb::v3::Seqdb::match

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::find_by_seq_ids(const std::vector<std::string_view>& seq_ids) const
{
    const auto& index = seq_id_index();
    subset result(seq_ids.size());
    std::transform(std::begin(seq_ids), std::end(seq_ids), result.begin(), [&index](std::string_view seq_id) -> ref {
        if (const auto found = index.find(seq_id); found != index.end())
            return found->second;
        else
            return {};
    });

    return result;

} // acmacs::seqdb::v3::Seqdb::find_by_seq_ids

// ----------------------------------------------------------------------

acmacs::seqdb::v3::Seqdb::aas_indexes_t acmacs::seqdb::v3::Seqdb::aa_at_pos1_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions1) const
{
    aas_indexes_t aas_indexes;
    acmacs::enumerate(match(aAntigens), [this,&aas_indexes,&aPositions1](auto ag_no, const auto& ref) {
        if (ref) {
            std::string aa(aPositions1.size(), 'X');
            std::transform(aPositions1.begin(), aPositions1.end(), aa.begin(), [this, ref = ref](size_t pos) { return ref.aa_at_pos(*this, acmacs::seqdb::pos1_t{pos}); });
            aas_indexes[aa].push_back(ag_no);
        }
    });
    return aas_indexes;

} // acmacs::seqdb::v3::Seqdb::aa_at_pos_for_antigens

// ----------------------------------------------------------------------

acmacs::seqdb::v3::Seqdb::clades_t acmacs::seqdb::v3::Seqdb::clades_for_name(std::string_view name, clades_for_name_inclusive inclusive) const
{
    clades_t result;
    bool clades_found = false;
    for (const auto& ref : select_by_name(name)) {
        const auto& seq = ref.seq().with_sequence(*this);
        if (inclusive == clades_for_name_inclusive::yes || !clades_found) {
            std::copy(std::begin(seq.clades), std::end(seq.clades), std::back_inserter(result));
        }
        else {
            result.erase(std::remove_if(std::begin(result), std::end(result), [&seq](const auto& clade) { return !seq.has_clade_master(clade); }), std::end(result));
        }
        clades_found |= !seq.clades.empty();
    }
    return result;

} // acmacs::seqdb::v3::Seqdb::clades_for_name

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Seqdb::add_clades(acmacs::chart::ChartModify& chart) const
{
    auto antigens = chart.antigens_modify();
    acmacs::enumerate(match(*antigens, chart.info()->virus_type(acmacs::chart::Info::Compute::Yes)), [&](auto ag_no, const auto& ref) {
        if (ref) {
            const auto& seq = ref.seq().with_sequence(*this);
            auto& antigen = antigens->at(ag_no);
            if (!seq.clades.empty()) {
                for (const auto& clade : seq.clades)
                    antigen.add_clade(std::string{clade});
            }
            else {
                antigen.add_clade("SEQUENCED");
            }
        }
    });

} // acmacs::seqdb::v3::Seqdb::add_clades

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::Seqdb::sequences_of_chart_for_ace_view_1(const acmacs::chart::Chart& chart) const
{
    struct stat_per_pos_t
    {
        ssize_t shannon_index = 0; // https://en.wikipedia.org/wiki/Diversity_index
        std::map<char, size_t> aa_count;
    };

    constexpr size_t max_num_pos = 1000;
    std::vector<stat_per_pos_t> stat_per_pos(max_num_pos);
    to_json::object json_antigens;
    acmacs::enumerate(match(*chart.antigens(), chart.info()->virus_type()), [&](auto ag_no, const auto& ref) {
        if (ref) {
            const auto sequence = ref.aa_aligned(*this);
            json_antigens << to_json::key_val{std::to_string(ag_no), *sequence};
            acmacs::enumerate(*sequence, [&stat_per_pos](size_t pos, char aa) {++stat_per_pos[pos].aa_count[aa]; }, 1UL);
        }
    });
    for (auto& per_pos : stat_per_pos) {
        const auto sum = std::accumulate(per_pos.aa_count.begin(), per_pos.aa_count.end(), 0UL, [](auto accum, const auto& entry) { return accum + entry.second; });
        const auto shannon_index = -std::accumulate(per_pos.aa_count.begin(), per_pos.aa_count.end(), 0.0, [sum = double(sum)](auto accum, const auto& entry) {
            const double p = entry.second / sum;
            return accum + p * std::log(p);
        });
        per_pos.shannon_index = std::lround(shannon_index * 100);
    }
    to_json::object json_per_pos;
    for (auto [pos, entry] : acmacs::enumerate(stat_per_pos)) {
        // if (entry.size() > 1) // && (entry.find('X') == entry.end() || entry.size() > 2))
        json_per_pos << to_json::key_val{std::to_string(pos), to_json::object{to_json::key_val{"shannon", entry.shannon_index}, to_json::key_val{"aa_count", to_json::object::from(entry.aa_count)}}};
    }
    return to_json::object{to_json::key_val{"sequences", to_json::object{to_json::key_val{"antigens", json_antigens}, to_json::key_val{"per_pos", json_per_pos}}}}.compact();

} // acmacs::seqdb::v3::Seqdb::sequences_of_chart_for_ace_view_1

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::Seqdb::sequences_of_chart_as_fasta(const acmacs::chart::Chart& chart) const
{
    auto antigens = chart.antigens();
    std::string fasta;
    acmacs::enumerate(match(*antigens, chart.info()->virus_type()), [&](auto ag_no, const auto& ref) {
        if (ref)
            fasta += fmt::format(">{}\n{}\n", antigens->at(ag_no)->full_name(), ref.nuc_aligned(*this));
    });
    return fasta;

} // acmacs::seqdb::v3::Seqdb::sequences_of_chart_as_fasta

// ----------------------------------------------------------------------

std::string_view acmacs::seqdb::v3::SeqdbEntry::host() const
{
    if (const auto ho = acmacs::virus::host(acmacs::virus::v2::name_t{name}); !ho.empty())
        return ho;
    else
        return "HUMAN";

} // acmacs::seqdb::v3::SeqdbEntry::host

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::SeqdbEntry::location() const
{
    return ::virus_name::location(acmacs::virus::v2::name_t{name});

} // acmacs::seqdb::v3::SeqdbEntry::location

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::SeqdbSeq& acmacs::seqdb::v3::SeqdbSeq::find_master(const Seqdb& seqdb) const
{
    for (const auto& ref : seqdb.select_by_name(master.name)) {
        for (const auto& seq : ref.entry->seqs) {
            if (seq.is_master() && seq.matches_without_name(master))
                return seq;
        }
    }
    throw std::runtime_error{"internal in SeqdbSeq::find_master: invalid master ref"};

} // acmacs::seqdb::v3::SeqdbSeq::find_master

// ----------------------------------------------------------------------

const std::vector<acmacs::seqdb::v3::ref>& acmacs::seqdb::v3::SeqdbSeq::find_slaves(const Seqdb& seqdb, std::string_view name) const
{
    if (!slaves) {
        slaves = std::make_unique<std::vector<ref>>();
        const master_ref_t self{name, annotations, reassortants.empty() ? std::string_view{} : reassortants.front(), passages.empty() ? std::string_view{} : passages.front()};
        for (const auto& ref : seqdb.select_by_name(name)) {
            if (ref.seq().matches_without_name(self))
                slaves->push_back(ref);
        }
    }
    return *slaves;

} // acmacs::seqdb::v3::SeqdbSeq::find_slaves

// ----------------------------------------------------------------------

acmacs::seqdb::seq_id_t acmacs::seqdb::v3::ref::seq_id() const
{
    auto source = ::string::join(" ", {entry->name, seq().designation()});
    if (entry->seqs.size() > 1 && seq_index > 0) {
        // there could be multiple seqs with the same designation, but seq_id must be unique, also garli does not like name duplicates
        std::vector<std::string> designations(entry->seqs.size());
        std::transform(std::begin(entry->seqs), std::end(entry->seqs), std::begin(designations), [](const auto& en) { return en.designation(); });
        if (std::count(std::begin(designations), std::end(designations), designations[seq_index]) > 1)
            source.append(fmt::format("_d{}", seq_index));
    }
    return make_seq_id(source);

} // acmacs::seqdb::v3::ref::seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::multiple_dates(bool do_filter)
{
    if (do_filter)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.entry->dates.size() < 2; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::multiple_dates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::subtype(const acmacs::uppercase& virus_type)
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

} // acmacs::seqdb::v3::subset::subtype

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::lineage(const acmacs::uppercase& lineage)
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

} // acmacs::seqdb::v3::subset::lineage

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::lab(const acmacs::uppercase& lab)
{
    if (!lab.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lab=static_cast<std::string_view>(lab)](const auto& en) { return !en.has_lab(lab); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::lab

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::whocc_lab(bool do_filter)
{
    if (do_filter)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return ! (en.has_lab("CDC") || en.has_lab("CRICK") || en.has_lab("NIID") || en.has_lab("VIDRL")); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::whocc_lab

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::host(const acmacs::uppercase& host)
{
    if (!host.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [host=static_cast<std::string_view>(host)](const auto& en) { return en.entry->host() != host; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::host

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::continent(const acmacs::uppercase& continent)
{
    if (!continent.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [continent=static_cast<std::string_view>(continent)](const auto& en) { return en.entry->continent != continent; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::continent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::country(const acmacs::uppercase& country)
{
    if (!country.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [country=static_cast<std::string_view>(country)](const auto& en) { return en.entry->country != country; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::country

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::clade(const Seqdb& seqdb, const acmacs::uppercase& clade)
{
    if (!clade.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [&seqdb,clade=static_cast<std::string_view>(clade)](const auto& en) { return !en.has_clade(seqdb, clade); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::clade

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent(size_t recent)
{
    if (recent > 0 && refs_.size() > recent) {
        sort_by_date_recent_first();
        refs_.erase(std::next(std::begin(refs_), static_cast<ssize_t>(recent)), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent_master(size_t recent_master)
{
    if (recent_master > 0) {
        keep_master_only();
        if (refs_.size() > recent_master) {
            sort_by_date_recent_first();
            refs_.erase(std::next(std::begin(refs_), static_cast<ssize_t>(recent_master)), std::end(refs_));
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent_matched(const std::vector<size_t>& recent_matched)
{
    if (recent_matched.size() > 1 && refs_.size() > recent_matched[0]) {
        sort_by_date_recent_first();
        const auto usable_size = std::remove_if(std::next(std::begin(refs_), static_cast<ssize_t>(recent_matched[0])), std::end(refs_), [](const auto& en) { return !en.has_hi_names(); }) - std::begin(refs_);
        refs_.erase(std::next(std::begin(refs_), std::min(usable_size, static_cast<ssize_t>(recent_matched[0] + recent_matched[1]))), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent_matched

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent_matched_master(const Seqdb& seqdb, const std::vector<size_t>& recent_matched_master)
{
    if (recent_matched_master.size() > 1) {
        keep_master_only();
        sort_by_date_recent_first();
        // if ref (master) has no hi names and one of its slaves has hi name, replace ref with slave that has hi names and return false
        // if ref (master) has no hi names and none of its slaves has hi name, return true (to remove from refs_)
        const auto without_hi_names = [&seqdb](auto& ref) {
            if (ref.has_hi_names())
                return false;   // keep it
            const auto& slaves = ref.find_slaves(seqdb);
            if (const auto slave_to_use = std::find_if(std::begin(slaves), std::end(slaves), [](const auto& slave) { return slave.has_hi_names(); }); slave_to_use != std::end(slaves)) {
                ref = *slave_to_use;
                return false;   // keep it
            }
            else
                return true;    // remove
        };
        const auto usable_size = std::remove_if(std::next(std::begin(refs_), static_cast<ssize_t>(recent_matched_master[0])), std::end(refs_), without_hi_names) - std::begin(refs_);
        refs_.erase(std::next(std::begin(refs_), std::min(usable_size, static_cast<ssize_t>(recent_matched_master[0] + recent_matched_master[1]))), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent_matched

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::keep_master_only()
{
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.is_master(); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::keep_master_only

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::random(size_t random)
{
    if (random > 0 && refs_.size() > random) {
        std::mt19937 generator{std::random_device()()};
        std::uniform_int_distribution<size_t> distribution(0, refs_.size() - 1);
        set_remove_marker(true);
        for (size_t no = 0; no < random; ++no)
            refs_[distribution(generator)].to_be_removed = false;
        remove_marked();
    }
    return *this;

} // acmacs::seqdb::v3::subset::random

// ----------------------------------------------------------------------

// Eu's algortihm of subsseting 2019-07-23

// 1. Find first group master sequence. I think good starting sequence
// is the most recent one that matched against hidb. Algorithm also
// prefers matched sequences to make more antigens marked in the sig
// pages.
//
// 2. Compute hamming distance between rest sequences and the master
// sequence, sort rest sequences by hamming distance, smaller first.
//
// 3. Find group end, i.e. first sequence that has hamming distance to
// the group master bigger than dist_threshold. Assign group no to
// this group. Sort group (keep group master first) by number of hi
// names (most number of names first) and by date (most recent first).
//
// 4. Next group master is the first sequence after group end. Repeat
// 2-3-4 until all sequences are processed.
//
// 5. Select masters (first sequences) of every group. If there are
// too many groups, more than output_size, then just used first
// output_size groups. If output_size > number of groups, select the
// second sequence in each group (if group size > 1). Do it until
// output_size sequences selected.

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::group_by_hamming_distance(const Seqdb& seqdb, size_t dist_threshold, size_t output_size)
{
    if (dist_threshold > 0) {
        const auto compute_hamming_distance = [&seqdb](std::string_view master_aa, auto first, auto last) {
            std::for_each(first, last, [&seqdb,master_aa](auto& ref) { ref.hamming_distance = hamming_distance(master_aa, ref.aa_aligned(seqdb)); });
        };

        const auto sort_by_hamming_distance = [](auto first, auto last) { std::sort(first, last, [](const auto& e1, const auto& e2) { return e1.hamming_distance < e2.hamming_distance; }); };

        const auto find_group_end = [dist_threshold](auto first, auto last) { return std::find_if(first, last, [dist_threshold](const auto& en) { return en.hamming_distance >= dist_threshold; }); };

        const auto assign_group_no = [](auto first, auto last, size_t group_no) { std::for_each(first, last, [group_no](auto& en) { en.group_no = group_no; }); };

        const auto sort_by_hi_names = [](auto first, auto last) {
            std::sort(first, last, [](const auto& e1, const auto& e2) {
                return e1.seq().hi_names.size() == e2.seq().hi_names.size() ? e1.entry->date() > e2.entry->date() : e1.seq().hi_names.size() > e2.seq().hi_names.size();
            });
        };

        // ----------------------------------------------------------------------

        std::iter_swap(std::begin(refs_), most_recent_with_hi_name());
        auto group_first = std::begin(refs_);
        acmacs::Counter<ssize_t> counter_group_size;
        for (size_t group_no = 1; group_first != std::end(refs_); ++group_no) {
            const auto group_master_aa_aligned = group_first->aa_aligned(seqdb);
            const auto group_second = std::next(group_first);
            // fmt::print("DEBUG: group {} master: {} {} rest size: {}\n", group_no, group_first->seq_id(), group_first->entry->date(), std::end(refs_) - group_first);
            compute_hamming_distance(group_master_aa_aligned, group_second, std::end(refs_));
            sort_by_hamming_distance(group_second, std::end(refs_));
            const auto group_last = find_group_end(group_second, std::end(refs_));
            assign_group_no(group_first, group_last, group_no);
            sort_by_hi_names(group_no == 1 ? group_second : group_first, group_last);
            counter_group_size.count(group_last - group_first);
            group_first = group_last;
        }
        // fmt::print(stderr, "DEBUG: (num-groups:group-size): {}\n", counter_group_size.report_sorted_max_first(" {second}:{first}"));
        // fmt::print(stderr, "DEBUG: total groups: {}\n", refs_.back().group_no);
        if (refs_.back().group_no > output_size) {
            // too many groups, take one seq from each group starting with group 1, ignore groups with high numbers (furtherst from the recent strain)
            size_t prev_group = 0;
            for (auto& ref : refs_) {
                if (ref.group_no == prev_group)
                    ref.to_be_removed = true;
                else {
                    prev_group = ref.group_no;
                    if (prev_group > output_size)
                        ref.to_be_removed = true;
                }
            }
        }
        else {
            // too few groups
            set_remove_marker(true);
            size_t to_keep = 0;
            size_t prev_to_keep = output_size;
            while (to_keep < output_size && prev_to_keep != to_keep) {
                prev_to_keep = to_keep;
                size_t group_no = 1;
                for (auto& ref : refs_) {
                    if (ref.group_no >= group_no && ref.to_be_removed) {
                        ref.to_be_removed = false;
                        ++to_keep;
                        group_no = ref.group_no + 1;
                    }
                    if (to_keep >= output_size)
                        break;
                }
                // fmt::print(stderr, "DEBUG: to_keep {} group_no {}\n", to_keep, group_no);
            }
        }
        remove_marked();
    }
    return *this;

} // acmacs::seqdb::v3::subset::group_by_hamming_distance

// ----------------------------------------------------------------------

// davipatti algorithm 2019-07-23 9:58
// > 1. pick a random strain, put in selection
// > 2. pick random strain. if it has a distance < d to anything in in selection then discard it. else, add it to selection.
// > 3. repeat 3 until you have as many strains, n, as you want, or until no more strains to pick
//
// Problems: need to prioritize picking hidb matched sequences.
//
// > parameter d would have to be tuned if d=0, this is just randomly
// > sampling strains if d is very high, only very dissimilar strains will
// > make it into selection, and selection would be small ideally d would
// > be as high as possible such that the number of strains in the
// > selection is close to n
//
// Looks like we need to use a search for d, i.e. we do not stop on
// finding n strains at the step 4 and have to find all to learn how many
// redundant strains there are. And then pick d producing number of
// strains closer to n (I guess having slightly more than n is better
// than having slightly less) and cut it, if necessary.
//
// > i foresee this algorithm being run initially to make a selection when
// > new sequences come in, repeat step 3 above, but just on new strains
// > so, original members stay in selection anything novel enough gets
// > added to the selection selection slowly grows over time
//
// No. The size of selection must be the same (as close to 4k as possible).

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::subset_by_hamming_distance_random(const Seqdb& seqdb, bool do_subset, size_t output_size)
{
    if (do_subset && !refs_.empty()) {
        std::mt19937 generator{std::random_device()()};
        const auto random_from = [&generator](auto first, auto last) {
            std::uniform_int_distribution<ssize_t> distribution(0, last - first - 1);
            return std::next(first, distribution(generator));
        };

        const auto minimal_distance_less_than = [&seqdb](auto first, auto last, std::string_view picked_aa, size_t distance_threshold) -> bool {
            return std::any_of(first, last, [&seqdb, picked_aa, distance_threshold](const auto& en) { return hamming_distance(picked_aa, en.aa_aligned(seqdb)) < distance_threshold; });
        };

        decltype(refs_) best_data;
        for (size_t distance_threshold = 1; distance_threshold < 10; ++distance_threshold) {
            auto data = refs_;
            std::iter_swap(std::begin(data), random_from(std::begin(data), std::end(data)));
            auto selection_start = std::begin(data), selection_end = std::next(selection_start), discarded_start = std::end(data);
            while (discarded_start > selection_end) {
                auto picked = random_from(selection_end, discarded_start);
                if (minimal_distance_less_than(selection_start, selection_end, picked->aa_aligned(seqdb), distance_threshold)) { // discard
                    --discarded_start;
                    std::iter_swap(discarded_start, picked);
                }
                else { // put into selection
                    std::iter_swap(selection_end, picked);
                    ++selection_end;
                }
            }
            fmt::print(stderr, "DEBUG: threshold: {} selection: {}\n", distance_threshold, selection_end - selection_start);
            if (static_cast<size_t>(selection_end - selection_start) < output_size)
                break;          // use previous (best_data)
            best_data.resize(static_cast<size_t>(selection_end - selection_start));
            std::copy(selection_start, selection_end, std::begin(best_data));
        }
        if (best_data.empty())
            throw std::runtime_error(fmt::format("subset_by_hamming_distance_random: best_data is empty"));
        const auto num_seqs = std::min(output_size, best_data.size());
        refs_.resize(num_seqs);
        std::copy(std::begin(best_data), std::next(std::begin(best_data), static_cast<ssize_t>(num_seqs)), std::begin(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::subset_by_hamming_distance_random

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_nuc_duplicates(bool do_remove, bool keep_hi_matched)
{
    if (do_remove) {
        // master sequences and hi matched (if requested) in the [std::begin(refs_), to_remove_candidates_start] range
        const auto to_remove_canditates_start =
            std::partition(std::begin(refs_), std::end(refs_), [keep_hi_matched](const auto& ref) { return ref.is_master() || (keep_hi_matched && ref.is_hi_matched()); });


#if 0
        refs_.erase(to_remove_canditates_start, std::end(refs_));
#else
        // move slave seq from [to_remove_canditates_start, std::end(refs_)] that reference to
        // a sequence in [std::begin(refs_), to_remove_candidates_start]
        // to the [to_remove_start, std::end(refs_)] range
        const auto to_remove_start = std::partition(to_remove_canditates_start, std::end(refs_), [beg = std::begin(refs_), end = to_remove_canditates_start](const auto& ref1) {
            return std::find_if(beg, end, [&ref1](const auto& ref2) { return ref2.matches(ref1.seq().master); }) == end;
        });

        refs_.erase(to_remove_start, std::end(refs_));
#endif

    }
    return *this;

} // acmacs::seqdb::v3::subset::remove_nuc_duplicates

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_nuc_duplicates(const Seqdb& seqdb, bool do_remove, bool keep_hi_matched)
// {
//     if (do_remove) {
//         const acmacs::Counter counter_nuc_length(refs_, [&seqdb](const auto& en) { return en.nuc_aligned_length(seqdb); });
//         const auto nuc = [&seqdb, nuc_common_length = counter_nuc_length.max().first](const auto& en) { return en.nuc_aligned(seqdb, nuc_common_length); };
//         const auto hi_names = [](const auto& en) { return en.seq().hi_names.size(); };
//         std::sort(std::begin(refs_), std::end(refs_), [=](const auto& e1, const auto& e2) {
//             if (const auto n1 = nuc(e1), n2 = nuc(e2); n1 == n2)
//                 return hi_names(e1) > hi_names(e2);
//             else
//                 return n1 < n2;
//         });
//         if (keep_hi_matched) {
//             refs_.erase(std::unique(std::begin(refs_), std::end(refs_), [=](const auto& e1, const auto& e2) { return nuc(e1) == nuc(e2) && (hi_names(e1) == 0 || hi_names(e2) == 0); }),
//                         std::end(refs_));
//         }
//         else {
//             refs_.erase(std::unique(std::begin(refs_), std::end(refs_), [=](const auto& e1, const auto& e2) { return nuc(e1) == nuc(e2); }), std::end(refs_));
//         }
//     }
//     return *this;

// } // acmacs::seqdb::v3::subset::remove_nuc_duplicates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::refs_t::iterator acmacs::seqdb::v3::subset::most_recent_with_hi_name()
{
    auto result = std::end(refs_);
    std::string_view date;
    for (auto refp = std::begin(refs_); refp != std::end(refs_); ++refp) {
        if (refp->has_hi_names() && refp->entry->date() > date) { // refs_[no].seq().reassortants.empty() &&
            result = refp;
            date = refp->entry->date();
            // fmt::print(stderr, "DEBUG: [{}] {}\n", date, result->full_name());
        }
    }
    return result;

} // acmacs::seqdb::v3::subset::most_recent_with_hi_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::with_hi_name(bool with_hi_name)
{
    if (with_hi_name)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.has_hi_names(); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::with_hi_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::aa_at_pos(const Seqdb& seqdb, const amino_acid_at_pos1_eq_list_t& aa_at_pos)
{
    if (!aa_at_pos.empty()) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&aa_at_pos,&seqdb](const auto& en) {
                                       try {
                                           const auto& seq = en.seq().with_sequence(seqdb);
                                           return seq.amino_acids.empty() || !seq.matches(aa_at_pos); // true to remove
                                       }
                                       catch (std::exception& err) {
                                           throw std::runtime_error{fmt::format("{}, full_name: {}", err, en.full_name())};
                                       }
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::aa_at_pos

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_at_pos(const Seqdb& seqdb, const nucleotide_at_pos1_eq_list_t& nuc_at_pos)
{
    if (!nuc_at_pos.empty()) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&nuc_at_pos,&seqdb](const auto& en) {
                                       try {
                                           const auto& seq = en.seq().with_sequence(seqdb);
                                           return seq.nucs.empty() || !seq.matches(nuc_at_pos); // true to remove
                                       }
                                       catch (std::exception& err) {
                                           throw std::runtime_error{fmt::format("{}, full_name: {}", err, en.full_name())};
                                       }
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_at_pos

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::names_matching_regex(const std::vector<std::string_view>& regex_list)
{
    if (!regex_list.empty()) {
        std::vector<std::regex> re_list(regex_list.size());
        std::transform(std::begin(regex_list), std::end(regex_list), std::begin(re_list),
                       [](const auto& regex_s) { return std::regex(std::begin(regex_s), std::end(regex_s), std::regex_constants::icase); });
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&re_list](const auto& en) {
                                       return std::none_of(std::begin(re_list), std::end(re_list), [full_name = en.full_name()](const auto& re) { return std::regex_search(full_name, re); });
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::names_matching_regex

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::dates(std::string_view start, std::string_view end)
{
    if (!start.empty() || !end.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [start,end](const auto& en) { return !en.entry->date_within(start, end); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::dates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend(std::string_view seq_id, const Seqdb& seqdb)
{
    if (!seq_id.empty()) {
        auto candidates = seqdb.select_by_seq_id(seq_id);
        if (candidates.empty())
            throw std::runtime_error{fmt::format("no sequences with seq-id \"{}\" found", seq_id)};
        refs_.erase(std::remove(std::begin(refs_), std::end(refs_), candidates.front()), std::end(refs_)); // remove it, if selected earlier
        refs_.insert(std::begin(refs_), candidates.front());
    }
    return *this;

} // prepend

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend(const std::vector<std::string_view>& seq_ids, const Seqdb& seqdb)
{
    if (!seq_ids.empty()) {
        auto candidates = seqdb.select_by_seq_id(seq_ids);
        if (candidates.empty())
            throw std::runtime_error{fmt::format("no sequences by seq-ids found to prepend")};
        const auto select_to_remove = [&candidates](const auto& ref) { return std::find(std::begin(candidates), std::end(candidates), ref) != std::end(candidates); };
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), select_to_remove), std::end(refs_)); // remove it, if selected earlier
        refs_.insert(std::begin(refs_), std::begin(candidates), std::end(candidates));
    }
    return *this;

} // prepend

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend_single_matching(std::string_view re, const Seqdb& seqdb)
// {
//     if (!re.empty()) {
//         auto candidates = seqdb.select_by_regex(re);
//         if (candidates.size() != 1) {
//             fmt::print(stderr, "WARNING: Selected sequences: {}\n", candidates.size());
//             for (const auto& candidate : candidates)
//                 fmt::print(stderr, "    {}\n", candidate.seq_id());
//             throw std::runtime_error{fmt::format("regular expression must select single sequence: \"{}\", selected: {}", re, candidates.size())};
//         }
//         refs_.erase(std::remove(std::begin(refs_), std::end(refs_), candidates.front()), std::end(refs_)); // remove it, if selected earlier
//         refs_.insert(std::begin(refs_), candidates.front());
//     }
//     return *this;

// } // acmacs::seqdb::v3::subset::prepend_single_matching

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base(size_t threshold, bool do_filter)
{
    if (do_filter) {
        const auto& seqdb = acmacs::seqdb::get();
        refs_.erase(std::remove_if(std::next(std::begin(refs_)), std::end(refs_),
                                   [threshold, &seqdb, base_seq = refs_.front().aa_aligned(seqdb)](auto& en) {
                                       en.hamming_distance = hamming_distance(en.aa_aligned(seqdb), base_seq);
                                       return en.hamming_distance >= threshold;
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::sort(sorting srt)
{
    switch (srt) {
        case sorting::none:
            break;
        case sorting::name_asc:
            sort_by_name_asc();
            break;
        case sorting::name_desc:
            sort_by_name_desc();
            break;
        case sorting::date_asc:
            sort_by_date_oldest_first();
            break;
        case sorting::date_desc:
            sort_by_date_recent_first();
            break;
    }
    return *this;

} // acmacs::seqdb::v3::subset::sort

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_stat(bool do_report)
{
    if (do_report) {
        if (!refs_.empty()) {
            size_t with_hi_names = 0;
            std::string_view min_date = refs_.front().entry->date(), max_date = min_date;
            for (const auto& ref : refs_) {
                const auto date = ref.entry->date();
                if (date < min_date)
                    min_date = date;
                else if (date > max_date)
                    max_date = date;
                if (!ref.seq().hi_names.empty())
                    ++with_hi_names;
            }
            fmt::print(stderr, "Sequences: {}\nDate range: {} - {}\nHiDb matches: {}\n", refs_.size(), min_date, max_date, with_hi_names);
        }
        else {
            fmt::print(stderr, "No sequences selected\n");
        }
    }

    return *this;

} // acmacs::seqdb::v3::subset::report_stat

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::export_sequences(std::string_view filename, const Seqdb& seqdb, const export_options& options)
{
    if (!filename.empty()) {
        auto to_export = export_collect(seqdb, options);

        if (options.e_most_common_length) {
            const acmacs::Counter counter(to_export, [](const auto& en) { return en.second.size(); });
            ranges::for_each(to_export, [most_common_length = counter.max().first](auto& en) { en.second.resize(most_common_length, '-'); });
        }

        switch (options.e_format) {
            case export_options::format::fasta_aa:
            case export_options::format::fasta_nuc:
                acmacs::file::write(filename, export_fasta(to_export, options));
                break;
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::export_sequences

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_hamming_distance(bool do_report, const Seqdb& seqdb)
{
    if (do_report) {
        auto to_report = export_collect(seqdb, {export_options::format::fasta_nuc, 0, true, true, "{seq_id}"});

        struct entry_t
        {
            std::string seq_id;
            std::string nucs;
            size_t hamming_distance;
        };

        std::vector<entry_t> data(to_report.size());
        std::transform(std::begin(to_report), std::end(to_report), std::begin(data), [&to_report](const auto& source) -> entry_t {
            return {source.first, source.second, hamming_distance(to_report.front().second, source.second)};
        });
        std::sort(std::begin(data), std::end(data), [](const auto& e1, const auto& e2) { return e1.hamming_distance > e2.hamming_distance; });
        for (const auto& en : data)
            fmt::print("{:4d}  {}\n", en.hamming_distance, en.seq_id);
    }
    return *this;

} // acmacs::seqdb::v3::subset::report_hamming_distance

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::make_name(const Seqdb& seqdb, std::string_view name_format, const ref& entry) const
{
    return fmt::format(name_format,
                       fmt::arg("seq_id", entry.seq_id()),
                       fmt::arg("full_name", entry.full_name()),
                       fmt::arg("hi_name_or_full_name", entry.hi_name_or_full_name()),
                       fmt::arg("hi_names", entry.seq().hi_names),
                       fmt::arg("hi_name", !entry.seq().hi_names.empty() ? entry.seq().hi_names.front() : std::string_view{}),
                       fmt::arg("lineage", entry.entry->lineage),
                       fmt::arg("name", entry.entry->name),
                       fmt::arg("date", entry.entry->date()),
                       fmt::arg("dates", entry.entry->dates),
                       fmt::arg("lab_id", entry.seq().lab_id()),
                       fmt::arg("passage", entry.seq().passage()),
                       fmt::arg("clades", entry.seq().with_sequence(seqdb).clades),
                       fmt::arg("lab", entry.seq().lab()),
                       fmt::arg("country", entry.entry->country),
                       fmt::arg("continent", entry.entry->continent),
                       fmt::arg("group_no", entry.group_no ? fmt::format("group:{}", entry.group_no) : std::string{}),
                       fmt::arg("hamming_distance", entry.hamming_distance)
                       );

} // acmacs::seqdb::v3::subset::make_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::collected_t acmacs::seqdb::v3::subset::export_collect(const Seqdb& seqdb, const export_options& options) const
{
    const auto get_seq = [&options,&seqdb](const auto& entry) -> std::string_view {
        const auto& seq = entry.seq().with_sequence(seqdb);
        if (options.e_format == export_options::format::fasta_aa) {
            if (options.e_aligned)
                return *seq.aa_aligned_master();
            else
                return std::get<std::string_view>(seq.amino_acids);
        }
        else {
            if (options.e_aligned)
                return *seq.nuc_aligned_master();
            else
                return std::get<std::string_view>(seq.nucs);
        }
    };

    collected_t result(refs_.size()); // {seq_id, sequence}
    std::transform(std::begin(refs_), std::end(refs_), std::begin(result),
                   [this, &options, &get_seq, &seqdb](const auto& en) -> collected_entry_t { return std::pair(make_name(seqdb, options.e_name_format, en), std::string{get_seq(en)}); });
    return result;

} // acmacs::seqdb::v3::subset::export_collect

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::export_fasta(const collected_t& entries, const export_options& options)
{
    std::string output;
    const auto output_size =
        std::accumulate(std::begin(entries), std::end(entries), 0UL, [](size_t size, const auto& en) { return size + en.first.size() + en.second.size() + 2 + en.second.size() / 40; });
    output.reserve(output_size);
    for (const auto& en : entries) {
        output.append(1, '>');
        output.append(en.first);
        output.append(1, '\n');
        if (options.e_wrap_at == 0 || options.e_wrap_at >= en.second.size()) {
            output.append(en.second);
            output.append(1, '\n');
        }
        else {
            for (const auto chunk : en.second | ranges::views::chunk(options.e_wrap_at)) {
                output.append(ranges::to<std::string>(chunk));
                output.append(1, '\n');
            }
        }
    }
    return output;

} // acmacs::seqdb::v3::subset::export_fasta

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::append(const subset& another)
{
    std::copy(std::begin(another), std::end(another), std::back_inserter(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::append

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::subset::filter_by_indexes(const acmacs::chart::PointIndexList& indexes, enum matched_only matched_only) const
{
    subset result;
    for (auto index : indexes) {
        if (index < refs_.size() && (matched_only == matched_only::no || refs_[index]))
            result.refs_.push_back(refs_[index]);
    }
    return result;

} // acmacs::seqdb::v3::subset::filter_by_indexes

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::print(const Seqdb& seqdb, std::string_view name_format, bool do_print)
{
    if (do_print) {
        for (const auto& ref : *this)
            fmt::print("{}\n", make_name(seqdb, name_format, ref));
    }
    return *this;

} // acmacs::seqdb::v3::subset::print

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
