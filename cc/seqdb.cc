#include <algorithm>
#include <random>
#include <regex>
#include <numeric>
#include <memory>
#include <cstdlib>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/in-json-parser.hh"
#include "acmacs-base/to-json.hh"
#include "acmacs-base/hash.hh"
#include "acmacs-base/date.hh"
#include "acmacs-base/string-matcher.hh"
#include "acmacs-virus/virus-name-normalize.hh"
#include "acmacs-virus/virus-name-v1.hh"
#include "acmacs-chart-2/chart-modify.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/seqdb-parse.hh"
#include "seqdb-3/hamming-distance.hh"
#include "seqdb-3/log.hh"

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
        find_slaves();
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

std::pair<acmacs::seqdb::v3::Seqdb::seq_id_iter, acmacs::seqdb::v3::Seqdb::seq_id_iter> acmacs::seqdb::v3::Seqdb::find_seq_id(std::string_view seq_id) const
{
    auto bounds = seq_id_index().find(seq_id);
    if (bounds.first == bounds.second && seq_id.size() > 3 && seq_id[seq_id.size() - 3] == '_' && seq_id[seq_id.size() - 2] == 'd')
        bounds = seq_id_index().find(seq_id.substr(0, seq_id.size() - 3));
    return bounds;

} // acmacs::seqdb::v3::Seqdb::find_seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_seq_id(std::string_view seq_id) const
{
    subset ss;
    if (const auto [first, last] = find_seq_id(seq_id); first != last)
        ss.refs_.emplace_back(first->second);
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_seq_id(const std::vector<std::string_view>& seq_ids) const
{
    subset ss;
    for (const auto& seq_id : seq_ids) {
        if (const auto [first, last] = find_seq_id(seq_id); first != last)
            ss.refs_.emplace_back(first->second);
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::find_by_seq_ids(const std::vector<std::string_view>& seq_ids) const
{
    subset result(seq_ids.size());
    std::transform(std::begin(seq_ids), std::end(seq_ids), result.begin(), [this](std::string_view seq_id) -> ref {
        if (const auto [first, last] = find_seq_id(seq_id); first != last)
            return first->second;
        else
            return {};
    });

    return result;

} // acmacs::seqdb::v3::Seqdb::find_by_seq_ids

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
    using namespace std::string_view_literals;

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
    if (subs.size() == subs_initial_size && (name[0] == 'A' || name[0] == 'a' || name[0] == 'B' || name[0] == 'b')) {
        const auto result = acmacs::virus::name::parse(name);
        find_name(result.name());
        if (subs.size() == subs_initial_size && (name[0] == 'A' || name[0] == 'a') && name[1] == '/') {
            for (const auto subtype : {"A(H1N1)/"sv, "A(H3N2)/"sv, "A(H1)/"sv, "A(H3)/"sv}) {
                find_name(acmacs::virus::name::parse(fmt::format("{}{}", subtype, name.substr(2))).name());
            }
        }
    }
    if (subs.size() == subs_initial_size) {
        for (const auto subtype : {"A(H1N1)/"sv, "A(H3N2)/"sv, "B/"sv, "A(H1)/"sv, "A(H3)/"sv}) {
            find_name(acmacs::virus::name::parse(fmt::format("{}{}", subtype, name)).name());
        }
    }

} // acmacs::seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_accession_number(const std::vector<std::string_view>& accession_numbers) const
{
    const auto intersect = [&accession_numbers](const std::vector<std::string_view> ids) {
        for (const auto& id : ids) {
            if (const auto found = std::find(std::begin(accession_numbers), std::end(accession_numbers), id); found != std::end(accession_numbers)) {
                // AD_DEBUG("select_by_accession_number {}", ids);
                return true;
            }
        }
        return false;
    };

    subset ss;
    for (const auto& entry : entries_) {
        for (auto [seq_no, seq] : acmacs::enumerate(entry.seqs)) {
            if (intersect(seq.gisaid.isolate_ids) || intersect(seq.gisaid.sample_ids_by_sample_provider))
                ss.refs_.emplace_back(&entry, seq_no);
        }
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_accession_number

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_name_hash(std::string_view name, std::string_view hash) const
{
    subset ss;
    if (auto [found_first, found_last] = hash_index().find(hash); found_first != found_last) {
        bool ref_found{false};
        for (; found_first != found_last; ++found_first) {
            if (found_first->second.entry->name == name) {
                ss.refs_.push_back(found_first->second);
                ref_found = true;
                // fmt::print(stderr, "DEBUG: select_by_name_hash {} {} -> {}\n", name, hash, found_first->second.full_name());
            }
        }
        if (!ref_found)
            fmt::print(stderr, "WARNING: Seqdb::select_by_name_hash: name difference for hash {}, no \"{}\"\n", hash, name);
    }
    // else
    //     fmt::print(stderr, "DEBUG: select_by_name_hash {} {} -> NOT FOUND\n", name, hash);
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_name_hash

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

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_slaves() const
{
    subset ss;
    for (const auto& entry : entries_) {
        for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
            if (ref candidate{&entry, seq_no}; !candidate.is_master())
                ss.refs_.push_back(std::move(candidate));
        }
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_slaves

// ----------------------------------------------------------------------

acmacs::seqdb::v3::ref acmacs::seqdb::v3::Seqdb::find_hi_name(std::string_view full_name) const
{
    if (const auto* ref = acmacs::seqdb::get().hi_name_index().find(full_name); ref)
        return *ref;
    else
        return {};

} // acmacs::seqdb::v3::Seqdb::find_hi_name

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::seq_id_index_t& acmacs::seqdb::v3::Seqdb::seq_id_index() const
{
    std::lock_guard<std::mutex> index_guard(index_access_);
    if (seq_id_index_.empty()) {
        for (const auto& entry : entries_) {
            for (auto [seq_no, seq] : acmacs::enumerate(entry.seqs)) {
                for (const auto& designation : seq.designations())
                    seq_id_index_.emplace(make_seq_id(acmacs::string::join(acmacs::string::join_space, entry.name, designation)), ref{entry, seq_no});
            }
        }
        seq_id_index_.sort();     // force sorting to avoid future raise condition during access from different threads
    }
    return seq_id_index_;

} // acmacs::seqdb::v3::Seqdb::seq_id_index

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::hi_name_index_t& acmacs::seqdb::v3::Seqdb::hi_name_index() const
{
    std::lock_guard<std::mutex> index_guard(index_access_);
    if (hi_name_index_.empty()) {
        for (const auto& entry : entries_) {
            for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
                for (const auto& hi_name : entry.seqs[seq_no].hi_names)
                    hi_name_index_.emplace(hi_name, ref{&entry, seq_no});
            }
        }
        hi_name_index_.sort();     // force sorting to avoid future raise condition during access from different threads
    }
    return hi_name_index_;

} // acmacs::seqdb::v3::Seqdb::hi_name_index

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::hash_index_t& acmacs::seqdb::v3::Seqdb::hash_index() const
{
    std::lock_guard<std::mutex> index_guard(index_access_);
    if (hash_index_.empty()) {
        using namespace ranges::views;
        hash_index_.collect(
            entries_
            | for_each([](const auto& entry) {
                return ranges::yield_from(
                    ints(0UL, entry.seqs.size())
                    | transform([&entry](auto seq_no) { return std::pair{entry.seqs[seq_no].hash, seq_no}; })
                    | filter([](const auto& hash_seq_no) { return !hash_seq_no.first.empty(); })
                    | transform([&entry](const auto& hash_seq_no) -> std::pair<std::string_view, ref> { return {hash_seq_no.first, ref{&entry, hash_seq_no.second}}; }));
            }));
        hash_index_.sort();     // force sorting to avoid future raise condition during access from different threads
    }
    return hash_index_;

} // acmacs::seqdb::v3::Seqdb::hash_index

// ----------------------------------------------------------------------

inline std::optional<acmacs::seqdb::v3::ref> match(const acmacs::seqdb::v3::subset& sequences, const acmacs::virus::Reassortant& ag_reassortant, const acmacs::virus::Passage& ag_passage)
{
    if (sequences.empty())
        return std::nullopt;
    std::vector<string_match::score_t> score_per_seq(sequences.size(), string_match::score_t{-1});
    for (size_t seq_no{0}; seq_no < sequences.size(); ++seq_no) {
        const auto& seq = sequences[seq_no].seq();

        AD_LOG(acmacs::log::hi_name_matching, "      {} -- reassortants:{} passages:{}", sequences[seq_no].seq_id(), seq.reassortants, seq.passages);
        if ((seq.reassortants.empty() && ag_reassortant.empty()) || std::any_of(std::begin(seq.reassortants), std::end(seq.reassortants), [&ag_reassortant](std::string_view reass) { return ag_reassortant == reass; })) {
            for (const auto& s_passage : seq.passages) {
                if (acmacs::virus::passages_match(ag_passage, acmacs::virus::Passage{s_passage})) {
                    if (const auto score = string_match::match(s_passage, *ag_passage); score > score_per_seq[seq_no]) {
                        AD_LOG(acmacs::log::hi_name_matching, "      score: {}", score);
                        score_per_seq[seq_no] = score;
                    }
                }
            }
        }
    }
    if (const auto best_seq = std::max_element(std::begin(score_per_seq), std::end(score_per_seq)); *best_seq >= 0)
        return sequences[static_cast<size_t>(best_seq - std::begin(score_per_seq))];
    else
        return std::nullopt;
}

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::Antigens& aAntigens, std::string_view /*aChartVirusType*/) const
{
    // check lineage?
    // check virus type

    subset result;

    auto find_by_hi_name = [this](const auto& antigen) -> std::optional<ref> {
        const auto& hi_name_ind = hi_name_index();
        if (const auto* found_ref1 = hi_name_ind.find(antigen.full_name()); found_ref1)
            return *found_ref1;
        else if (const auto* found_ref2 = hi_name_ind.find(antigen.full_name_for_seqdb_matching()); found_ref2)
            return *found_ref2;
        else
            return std::nullopt;
    };

    size_t num_matched = 0;
    for (auto antigen : aAntigens) {
        if (auto found_ref = find_by_hi_name(*antigen); found_ref.has_value()) {
            result.refs_.push_back(std::move(*found_ref));
            ++num_matched;
        }
        else {
            const auto name_fields = acmacs::virus::name::parse(antigen->name());
            if (name_fields.mutations.empty()) {
                const acmacs::virus::Reassortant ag_reassortant{antigen->reassortant().empty() ? name_fields.reassortant : antigen->reassortant()};
                const acmacs::virus::Passage ag_passage{antigen->passage().empty() ? name_fields.passage : antigen->passage()};
                const auto sequences{select_by_name(name_fields.name())};
                AD_LOG(acmacs::log::hi_name_matching, "select_by_name \"{}\" ({}) \"{}\" sequences:{}", antigen->name(), name_fields.name(), antigen->full_name(), sequences.size());
                if (const auto matched = ::match(sequences, ag_reassortant, ag_passage); matched.has_value()) {
                    AD_LOG(acmacs::log::hi_name_matching, "  --> {}", matched->seq_id());
                    result.refs_.push_back(*matched);
                    ++num_matched;
                }
                else
                    result.refs_.emplace_back();
            }
            else // contains mutation data, nothing to look in seqdb
                result.refs_.emplace_back();
        }
    }
    fmt::print("INFO: antigens from chart have sequences in seqdb: {}\n", num_matched);

    return result;

} // acmacs::seqdb::v3::Seqdb::match

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

void acmacs::seqdb::v3::Seqdb::add_clades(acmacs::chart::ChartModify& chart, verbose verb) const
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
            if (verb == verbose::yes)
                fmt::print(stderr, "DEBUG: Seqdb::add_clades AG {:4d} {} -- {}\n", ag_no, antigen.full_name(), antigen.clades());
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
            const double p = static_cast<double>(entry.second) / sum;
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

// cannot return std::string_view because ho below would point to temporary
std::string acmacs::seqdb::v3::SeqdbEntry::host() const
{
    if (const std::string ho{acmacs::virus::host(acmacs::virus::v2::name_t{name})}; !ho.empty())
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

std::string_view acmacs::seqdb::v3::SeqdbEntry::date() const
{
    if (!dates.empty())
        return dates.front();
    if (name.size() > 5 && name[name.size() - 5] == '/') {
        if (const auto year = static_cast<size_t>(static_cast<int>(date::year_from_string(name.substr(name.size() - 4)))); year > 1900 && year <= date::current_year())
            return name.substr(name.size() - 4);
    }
    return std::string_view{};

} // acmacs::seqdb::v3::SeqdbEntry::date

// ----------------------------------------------------------------------


void acmacs::seqdb::v3::Seqdb::find_slaves() const
{
    if (!slaves_found_) {
        for (const auto& slave : select_slaves())
            slave.seq().find_master(*this).add_slave(slave);
        slaves_found_ = true;
    }

} // acmacs::seqdb::v3::Seqdb::find_slaves

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::SeqdbSeq& acmacs::seqdb::v3::SeqdbSeq::find_master(const Seqdb& seqdb) const
{
    if (master.name.empty())
        throw std::runtime_error{fmt::format("internal in SeqdbSeq::find_master: not a slave (name empty): {} {}", master.name, master.hash)}; // master.annotations, master.reassortant, master.passage)};

    for (const auto& ref : seqdb.select_by_name_hash(master.name, master.hash)) {
        if (ref)
            return ref.seq();
    }

    // for (const auto& ref : seqdb.select_by_name(master.name)) {
    //     for (const auto& seq : ref.entry->seqs) {
    //         if (seq.is_master() && seq.matches_without_name(master))
    //             return seq;
    //     }
    // }

    throw std::runtime_error{fmt::format("internal in SeqdbSeq::find_master: invalid master ref: {} {}", master.name, master.hash)}; //, master.annotations, master.reassortant, master.passage)};

} // acmacs::seqdb::v3::SeqdbSeq::find_master

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::SeqdbSeq::add_slave(const ref& slave) const
{
    if (!slaves_)
        slaves_ = std::make_unique<std::vector<ref>>();
    slaves_->push_back(slave);

} // acmacs::seqdb::v3::SeqdbSeq::add_slave

// ----------------------------------------------------------------------

const std::vector<acmacs::seqdb::v3::ref>& acmacs::seqdb::v3::SeqdbSeq::slaves() const
{
    if (!slaves_)
        slaves_ = std::make_unique<std::vector<ref>>();
    return *slaves_;

} // acmacs::seqdb::v3::SeqdbSeq::slaves

// ----------------------------------------------------------------------

// returns designations with and without hash
std::vector<std::string> acmacs::seqdb::v3::SeqdbSeq::designations(bool just_first) const
{
    const auto prefix = acmacs::string::join(acmacs::string::join_space, annotations, string::join(acmacs::string::join_space, reassortants));
    std::string my_hash{hash};
    if (my_hash.empty() && !is_master())
        my_hash = master.hash;
    const auto prefixed_hash = fmt::format("h{}", my_hash);
    if (passages.empty()) {
        return {acmacs::string::join(acmacs::string::join_space, prefix, prefixed_hash), prefix}; // not seq-id with hash must be first to support just_first
    }
    else if (just_first) {
        return {acmacs::string::join(acmacs::string::join_space, prefix, passages.front(), prefixed_hash)};
    }
    else {
        using namespace ranges::views;
        using hash_t = decltype(prefixed_hash);
        const std::array hashes{prefixed_hash, hash_t{}};
        return passages
                | for_each([&prefix,&hashes](std::string_view psg) {
                    return ranges::yield_from(
                        hashes
                        | transform([&prefix,psg](std::string_view a_hash) -> std::string { return acmacs::string::join(acmacs::string::join_space, prefix, psg, a_hash); }));
                })
                | ranges::to<std::vector>
                | ranges::actions::sort
                | ranges::actions::unique;
    }

} // acmacs::seqdb::v3::SeqdbSeq::designations

// ----------------------------------------------------------------------

acmacs::seqdb::seq_id_t acmacs::seqdb::v3::ref::seq_id() const
{
    auto source = acmacs::string::join(acmacs::string::join_space, entry->name, seq().designation());
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

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent(size_t recent, master_only master)
{
    if (recent > 0) {
        if (master == master_only::yes)
            keep_master_only();
        if (refs_.size() > recent) {
            sort_by_date_recent_first();
            refs_.erase(std::next(std::begin(refs_), static_cast<ssize_t>(recent)), std::end(refs_));
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent_matched(const std::vector<size_t>& recent_matched, master_only master)
{
    if (recent_matched.size() > 1 && refs_.size() > recent_matched[0]) {
        if (recent_matched.size() != 2)
            throw std::runtime_error{fmt::format("invalid recent-matched specification: {} {}", recent_matched, recent_matched.size())};
        if (master == master_only::yes)
            keep_master_only();
        if ((recent_matched[0] + recent_matched[1]) < refs_.size()) {
            sort_by_date_recent_first();
            if (master == master_only::yes) {
                // if ref (master) has no hi names and one of its slaves has hi name, replace ref with slave that has hi names and return false
                // if ref (master) has no hi names and none of its slaves has hi name, return true (to remove from refs_)
                size_t number_to_keep = recent_matched[1];
                const auto without_hi_names = [&number_to_keep, remove = true, keep = false](auto& ref) {
                    if (number_to_keep == 0)
                        return remove;
                    if (ref.has_hi_names()) {
                        --number_to_keep;
                        return keep;
                    }
                    const auto& slaves = ref.seq().slaves();
                    if (const auto slave_to_use = std::find_if(std::begin(slaves), std::end(slaves), [](const auto& slave) { return slave.has_hi_names(); }); slave_to_use != std::end(slaves)) {
                        //     ref = *slave_to_use;
                        --number_to_keep;
                        return keep;
                    }
                    else
                        return remove;
                };

                const auto end = std::remove_if(std::next(std::begin(refs_), static_cast<ssize_t>(recent_matched[0])), std::end(refs_), without_hi_names);
                refs_.erase(end, std::end(refs_));
            }
            else {
                const auto usable_size =
                    std::remove_if(std::next(std::begin(refs_), static_cast<ssize_t>(recent_matched[0])), std::end(refs_), [](const auto& en) { return !en.has_hi_names(); }) - std::begin(refs_);
                refs_.erase(std::next(std::begin(refs_), std::min(usable_size, static_cast<ssize_t>(recent_matched[0] + recent_matched[1]))), std::end(refs_));
            }
        }
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

void acmacs::seqdb::v3::subset::remove(ref_indexes& to_remove)
{
    std::sort(std::begin(to_remove), std::end(to_remove));
    const auto rm_end = std::unique(std::begin(to_remove), std::end(to_remove));
    auto rm_iter = std::begin(to_remove);
    size_t current_index{0};
    const auto remove_predicate = [&current_index,&rm_iter,rm_end](const auto&) {
        if (rm_iter != rm_end && *rm_iter == current_index++) {
            ++rm_iter;
            return true;
        }
        else
            return false;
    };
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), remove_predicate), std::end(refs_));

} // acmacs::seqdb::v3::subset::remove

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::subset::keep(ref_indexes& to_keep)
{
    std::sort(std::begin(to_keep), std::end(to_keep));
    const auto keep_end = std::unique(std::begin(to_keep), std::end(to_keep));
    auto keep_iter = std::begin(to_keep);
    size_t current_index{0};
    const auto remove_predicate = [&current_index,&keep_iter,keep_end](const auto&) {
        if (keep_iter != keep_end && *keep_iter == current_index++) {
            ++keep_iter;
            return false;
        }
        else
            return true;
    };
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), remove_predicate), std::end(refs_));

} // acmacs::seqdb::v3::subset::keep

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::random(size_t random)
{
    if (random > 0 && refs_.size() > random) {
        std::mt19937 generator{std::random_device()()};
        std::uniform_int_distribution<size_t> distribution(0, refs_.size() - 1);
        ref_indexes to_keep(random);
        std::generate_n(to_keep.begin(), random, [&distribution,&generator]() { return distribution(generator); });
        keep(to_keep);
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
            ref_indexes to_remove;
            size_t prev_group = 0;
            for (auto [index, ref] : acmacs::enumerate(refs_)) {
                if (ref.group_no == prev_group)
                    to_remove.push_back(index);
                else {
                    prev_group = ref.group_no;
                    if (prev_group > output_size)
                        to_remove.push_back(index);
                }
            }
            remove(to_remove);
        }
        else {
            // too few groups
            ref_indexes to_keep_indexes;
            size_t to_keep = 0;
            size_t prev_to_keep = output_size;
            while (to_keep < output_size && prev_to_keep != to_keep) {
                prev_to_keep = to_keep;
                size_t group_no = 1;
                for (auto [index, ref] : acmacs::enumerate(refs_)) {
                    if (ref.group_no >= group_no) {
                        to_keep_indexes.push_back(index);
                        ++to_keep;
                        group_no = ref.group_no + 1;
                    }
                    if (to_keep >= output_size)
                        break;
                }
                // fmt::print(stderr, "DEBUG: to_keep {} group_no {}\n", to_keep, group_no);
            }
            keep(to_keep_indexes);
        }
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

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_empty(const Seqdb& seqdb, bool nuc)
{
    const auto is_empty = [&seqdb, nuc](const auto& ref) {
        const auto& seq = ref.seq_with_sequence(seqdb);
        // AD_LOG(acmacs::log::sequences, "      master aa:{} nuc:{} orig:{}", seq.aa_aligned_length_master(), seq.nuc_aligned_length_master(), ref.seq_id());
        return nuc ? seq.nuc_aligned_length_master() == 0 : seq.aa_aligned_length_master() == 0;
    };

    AD_LOG(acmacs::log::sequences, "removing empty ({}) from {} sequences", nuc ? "nuc" : "aa", refs_.size());
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), is_empty), std::end(refs_));
    AD_LOG(acmacs::log::sequences, "    {} sequences left", refs_.size());
    return *this;

} // acmacs::seqdb::v3::subset::remove_empty

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

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::min_aa_length(const Seqdb& seqdb, size_t length)
{
    if (length) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [length, &seqdb](const auto& en) { return en.aa_aligned_length(seqdb) < length; }), std::end(refs_));
    }
    return *this;


} // acmacs::seqdb::v3::subset::min_aa_length

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::min_nuc_length(const Seqdb& seqdb, size_t length)
{
    if (length) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [length, &seqdb](const auto& en) { return en.nuc_aligned_length(seqdb) < length; }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::min_nuc_length

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_with_front_back_deletions(const Seqdb& seqdb, bool remove, size_t length)
{
    if (remove) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [length, &seqdb](const auto& en) {
            const auto nucs = en.nuc_aligned(seqdb);
            if (nucs.at(pos1_t{1}) == '-')
                return true;
            if (length > 0 && (nucs.size() < pos0_t{length} || nucs.at(pos1_t{length}) == '-'))
                return true;    // too short or has deletion in the last nuc
            return false;
        }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::remove_with_front_back_deletions

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

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_to(size_t threshold, std::string_view seq_id)
{
    if (!seq_id.empty()) {
        const auto& seqdb = acmacs::seqdb::get();
        const auto compare_to = seqdb.select_by_seq_id(seq_id);
        if (compare_to.empty())
            throw std::runtime_error{fmt::format("no sequences with seq-id \"{}\" found", seq_id)};
        const auto before{refs_.size()};
        refs_.erase(std::remove_if(std::next(std::begin(refs_)), std::end(refs_),
                                   [threshold, &seqdb, comapre_to_seq = compare_to.front().nuc_aligned(seqdb)](auto& en) {
                                       en.hamming_distance = hamming_distance(en.nuc_aligned(seqdb), comapre_to_seq, hamming_distance_by_shortest::no);
                                       return en.hamming_distance >= threshold;
                                   }),
                    std::end(refs_));
        const auto after{refs_.size()};
        AD_LOG(acmacs::log::sequences, "{} sequences removed ({} left) which are too far from {}, threshold: {}", before - after, after, seq_id, threshold);
        if ((before - after) > (before / 4))
            AD_WARNING("too many sequences removed ({} or {:.1f}%) that are too far from {}, hamming distance threshold: {}", before - after, static_cast<double>(before - after) / static_cast<double>(before) * 100.0, seq_id, threshold);
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_hamming_distance_to

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base(size_t threshold, bool do_filter)
{
    if (do_filter) {
        const auto& seqdb = acmacs::seqdb::get();
        const auto before{refs_.size()};
        refs_.erase(std::remove_if(std::next(std::begin(refs_)), std::end(refs_),
                                   [threshold, &seqdb, base_seq = refs_.front().nuc_aligned(seqdb)](auto& en) {
                                       en.hamming_distance = hamming_distance(en.nuc_aligned(seqdb), base_seq, hamming_distance_by_shortest::no);
                                       return en.hamming_distance >= threshold;
                                   }),
                    std::end(refs_));
        const auto after{refs_.size()};
        AD_LOG(acmacs::log::sequences, "{} sequences removed ({} left) which are too far from the base seq, threshold: {}", before - after, after, threshold);
        if ((before - after) > (before / 4))
            AD_WARNING("too many sequences removed ({} or {:.1f}%) that are too far from the base sequence, hamming distance threshold: {}", before - after, static_cast<double>(before - after) / static_cast<double>(before) * 100.0, threshold);
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

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_stat(const Seqdb& seqdb, bool do_report)
{
    if (do_report) {
        if (!refs_.empty()) {
            size_t with_hi_names = 0;
            std::string_view min_date = refs_.front().entry->date(), max_date = min_date;
            Counter<std::string> by_year;
            Counter<size_t> aa_length, nuc_length;
            for (const auto& ref : refs_) {
                const auto date = ref.entry->date();
                if (date < min_date)
                    min_date = date;
                else if (date > max_date)
                    max_date = date;
                if (date.size() >= 4)
                    by_year.count(date.substr(0, 4));
                if (!ref.seq().hi_names.empty())
                    ++with_hi_names;
                aa_length.count(ref.seq_with_sequence(seqdb).aa_aligned_length_master());
                nuc_length.count(ref.seq_with_sequence(seqdb).nuc_aligned_length_master());
            }
            fmt::print(stderr, "Selected sequences: {:6d}\n      HiDb matches: {:6d}\n        Date range: {} - {}\n", refs_.size(), with_hi_names, min_date, max_date);
            fmt::print(stderr, "         AA length:{}\nNucleotide lengths:{}\n", aa_length.report_sorted_max_first(" {first}:{second}"), nuc_length.report_sorted_max_first(" {first}:{second}"));
            fmt::print(stderr, "           by Year:{}\n", by_year.report_sorted_max_first(" {first}:{second}"));
        }
        else {
            fmt::print(stderr, "No sequences selected\n");
        }
    }

    return *this;

} // acmacs::seqdb::v3::subset::report_stat

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_aa_at(const Seqdb& seqdb, const pos1_list_t& pos1_list)
{
    if (!pos1_list.empty() && !refs_.empty()) {
        std::vector<CounterChar> counters(pos1_list.size());
        for (const auto& ref : refs_) {
            for (auto index : acmacs::range(pos1_list.size()))
                counters[index].count(ref.aa_at_pos(seqdb, pos1_list[index]));
        }
        fmt::print(stderr, "AA at pos stat:\n");
        for (auto index : acmacs::range(pos1_list.size()))
            fmt::print(stderr, "  {}\n{}", pos1_list[index], counters[index].report_sorted_max_first(fmt::format("    {:3d}{{first}}  {{second:5d}}\n", pos1_list[index])));
    }
    return *this;

} // acmacs::seqdb::v3::subset::report_aa_at

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::export_sequences(std::string_view filename, const Seqdb& seqdb, const export_options& options)
{
    if (!filename.empty()) {
        auto to_export = export_collect(seqdb, options);

        if (options.e_most_common_length == export_options::most_common_length::yes) {
            const acmacs::Counter counter(to_export, [](const auto& en) { return en.sequence.size(); });
            const auto most_common_length = counter.max().first;
            AD_LOG(acmacs::log::fasta, "most common length: {}", most_common_length);
            ranges::for_each(to_export, [most_common_length](auto& en) { en.sequence.resize(most_common_length, '-'); });
        }
        else if (options.e_length > 0) {
            AD_LOG(acmacs::log::fasta, "sequence length for exporting: {}", options.e_length);
            ranges::for_each(to_export, [length = options.e_length](auto& en) { en.sequence.resize(length, '-'); });
        }

        ranges::for_each(to_export, [deletion_report_threshold=options.e_deletion_report_threshold](auto& en) {
            const auto dels = static_cast<size_t>(ranges::count_if(en.sequence, [](char nuc_aa) { return nuc_aa == '-' || nuc_aa == 'X'; }));
            const auto dels_at_the_end = en.sequence.back() == '-' || en.sequence.back() == 'X';
            if (dels_at_the_end || dels > deletion_report_threshold)
                AD_WARNING("{}: {} deletions or unknown AAs or deletions at the end", en.seq_id, dels);
        });

        AD_LOG(acmacs::log::fasta, "writing {} sequences to {}", to_export.size(), filename);
        acmacs::file::write(filename, export_fasta(to_export, options));
    }
    return *this;

} // acmacs::seqdb::v3::subset::export_sequences

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_hamming_distance(bool do_report)
{
    if (do_report) {
        std::vector<const ref*> refs(refs_.size());
        std::transform(std::begin(refs_), std::end(refs_), std::begin(refs), [](const auto& rr) { return &rr; });
        std::sort(std::begin(refs), std::end(refs), [](const ref* r1, const ref* r2) { return r1->hamming_distance > r2->hamming_distance; });
        for (const auto& en : refs)
            fmt::print("{:4d}  {}\n", en->hamming_distance, en->seq_id());
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
                       fmt::arg("hamming_distance", entry.hamming_distance),
                       fmt::arg("nuc_length", entry.seq().nuc_aligned_length_master()),
                       fmt::arg("aa_length", entry.seq().aa_aligned_length_master()),
                       fmt::arg("gisaid_accession_numbers", acmacs::string::join(acmacs::string::join_sep_t{"|"}, entry.seq().gisaid.isolate_ids)),
                       fmt::arg("ncbi_accession_numbers", acmacs::string::join(acmacs::string::join_sep_t{"|"}, entry.seq().gisaid.sample_ids_by_sample_provider))
                       );

} // acmacs::seqdb::v3::subset::make_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::collected_t acmacs::seqdb::v3::subset::export_collect(const Seqdb& seqdb, const export_options& options) const
{
    const auto get_seq = [&options,&seqdb](const auto& entry) -> std::string_view {
        const auto& seq = entry.seq().with_sequence(seqdb);
        AD_LOG(acmacs::log::fasta, "{} has-seq:{}", entry.seq_id(), entry.is_master());
        if (!entry.is_master())
            AD_LOG(acmacs::log::fasta, "    ref:({} {})", entry.seq().master.name, entry.seq().master.hash);
        AD_LOG(acmacs::log::fasta, "    aa:{} nuc:{}", seq.aa_aligned_length_master(), seq.nuc_aligned_length_master());
        if (options.e_format == export_options::format::fasta_aa) {
            if (options.e_aligned == export_options::aligned::yes)
                return *seq.aa_aligned_master();
            else
                return std::get<std::string_view>(seq.amino_acids);
        }
        else {
            if (options.e_aligned == export_options::aligned::yes)
                return *seq.nuc_aligned_master();
            else
                return std::get<std::string_view>(seq.nucs);
        }
    };

    collected_t result(refs_.size()); // {seq_id, sequence}
    std::transform(std::begin(refs_), std::end(refs_), std::begin(result),
                   [this, &options, &get_seq, &seqdb](const auto& en) -> collected_entry_t { return {make_name(seqdb, options.e_name_format, en), std::string{get_seq(en)}}; });
    // remove entries with empty sequences
    result.erase(std::remove_if(std::begin(result), std::end(result), [](const auto& en) { return en.sequence.empty(); }), std::end(result));
    AD_LOG(acmacs::log::fasta, "collected for exporting: {}", result.size());
    return result;

} // acmacs::seqdb::v3::subset::export_collect

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::export_fasta(const collected_t& entries, const export_options& options)
{
    std::string output;
    const auto output_size =
        std::accumulate(std::begin(entries), std::end(entries), 0UL, [](size_t size, const auto& en) { return size + en.seq_id.size() + en.sequence.size() + 2 + en.sequence.size() / 40; });
    output.reserve(output_size);
    for (const auto& en : entries) {
        output.append(1, '>');
        output.append(en.seq_id);
        output.append(1, '\n');
        if (options.e_wrap_at == 0 || options.e_wrap_at >= en.sequence.size()) {
            output.append(en.sequence);
            output.append(1, '\n');
        }
        else {
            for (const auto chunk : en.sequence | ranges::views::chunk(options.e_wrap_at)) {
                output.append(ranges::to<std::string>(chunk));
                output.append(1, '\n');
            }
        }
    }
    fmt::print("INFO: exported to fasta: {}\n", entries.size());
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
