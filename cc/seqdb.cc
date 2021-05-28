#include <algorithm>
#include <random>
#include <regex>
#include <numeric>
#include <memory>
#include <cstdlib>

#include "acmacs-base/read-file.hh"
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
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/seqdb-parse.hh"
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

const acmacs::seqdb::v3::lab_id_index_t& acmacs::seqdb::v3::Seqdb::lab_id_index() const
{
    std::lock_guard<std::mutex> index_guard(index_access_);
    if (lab_id_index_.empty()) {
        for (const auto& entry : entries_) {
            for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
                for (const auto& [lab, lab_ids] : entry.seqs[seq_no].lab_ids) {
                    for (const auto& lab_id : lab_ids) {
                        const auto lab_and_id = fmt::format("{}#{}", lab, lab_id);
                        lab_id_index_.emplace(lab_and_id, ref{&entry, seq_no});
                    }
                }
            }
        }
        lab_id_index_.sort();     // force sorting to avoid future raise condition during access from different threads
        // duplicates are possible!
    }
    return lab_id_index_;

} // acmacs::seqdb::v3::Seqdb::lab_id_index

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

        AD_LOG(acmacs::log::hi_name_matching, "{} R:{} P:{}", sequences[seq_no].seq_id(), seq.reassortants, seq.passages);
        AD_LOG_INDENT;
        if ((seq.reassortants.empty() && ag_reassortant.empty()) ||
            std::any_of(std::begin(seq.reassortants), std::end(seq.reassortants), [&ag_reassortant](std::string_view reass) { return ag_reassortant == reass; })) {
            if (!seq.passages.empty()) {
                for (const auto& s_passage : seq.passages) {
                    if (acmacs::virus::passages_match(ag_passage, acmacs::virus::Passage{s_passage})) {
                        const auto score = string_match::match(s_passage, *ag_passage);
                        score_per_seq[seq_no] = std::max(score_per_seq[seq_no], score);
                        AD_LOG(acmacs::log::hi_name_matching, "score: {} P:{}", score, s_passage);
                    }
                }
            }
            else {
                score_per_seq[seq_no] = ag_passage.empty() ? 2 : 1;
                AD_LOG(acmacs::log::hi_name_matching, "score: {} seq has no passage", score_per_seq[seq_no]);
            }
        }
        else
            AD_LOG(acmacs::log::hi_name_matching, "reassortant mismatch");
    }
    if (const auto best_seq = std::max_element(std::begin(score_per_seq), std::end(score_per_seq)); *best_seq >= 0)
        return sequences[static_cast<size_t>(best_seq - std::begin(score_per_seq))];
    else
        return std::nullopt;
}

// ----------------------------------------------------------------------

template <typename AgSr> acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const AgSr& antigens_sera, std::string_view /*aChartVirusType*/) const
{
    // check lineage?
    // check virus type

    subset result;

    auto find_by_hi_name = [this](const auto& antigen) -> std::optional<ref> {
        const auto& hi_name_ind = hi_name_index();
        if (const auto* found_ref1 = hi_name_ind.find(antigen.format("{name_full}")); found_ref1)
            return *found_ref1;
        else if (const auto* found_ref2 = hi_name_ind.find(antigen.format("{name}{ }{reassortant}{ }{passage}{ }{annotations}")); found_ref2)
            return *found_ref2;
        else
            return std::nullopt;
    };

    auto find_by_parsed_name = [&](const auto& antigen) -> std::optional<ref> {
        if (const auto name_fields = acmacs::virus::name::parse(antigen.name()); name_fields.mutations.empty()) {
            const acmacs::virus::Reassortant ag_reassortant{antigen.reassortant().empty() ? name_fields.reassortant : antigen.reassortant()};
            const acmacs::virus::Passage ag_passage{antigen.passage().empty() ? name_fields.passage : antigen.passage()};
            const auto sequences{select_by_name(name_fields.name())};
            AD_LOG(acmacs::log::hi_name_matching, "match find_by_parsed_name \"{}\" ({}) \"{}\" sequences:{}", antigen.name(), name_fields.name(), antigen.format("{name_full}"), sequences.size());
            AD_LOG_INDENT;
            if (const auto matched = ::match(sequences, ag_reassortant, ag_passage); matched.has_value()) {
                AD_LOG(acmacs::log::hi_name_matching, "--> {}", matched->seq_id());
                return *matched;
            }
        }
        return std::nullopt;
    };

    auto find_by_lab_id = [&](const auto& lab_id, const auto& antigen) -> std::optional<ref> {
        if (const auto [first, last] = lab_id_index().find(lab_id); first != last) {
            if (std::distance(first, last) == 1)
                return first->second;
            // multiple refs for the same lab_id
            subset sequences;
            for (auto r1 = first; r1 != last; ++r1) {
                // AD_DEBUG("{} {} -> {}", antigen.name_full(), r1->first, r1->second.seq_id());
                sequences.append(r1->second);
            }
            if (const auto matched = ::match(sequences, antigen.reassortant(), antigen.passage()); matched.has_value()) {
                // AD_DEBUG("--> {}", matched->seq_id());
                return *matched;
            }
            AD_WARNING("multiple refs for {} {} (first is selected): ({}) {}", antigen.name_full(), lab_id, sequences.size(), sequences);
            return first->second;
        }
        else
            return std::nullopt;
    };

    size_t num_matched = 0;
    for (auto antigen : antigens_sera) {
        std::optional<ref> found_ref{std::nullopt};
        if constexpr (std::is_same_v<AgSr, chart::Antigens> || std::is_same_v<AgSr, chart::AntigensModify>) {
            for (const auto& lab_id : antigen->lab_ids()) {
                found_ref = find_by_lab_id(lab_id, *antigen);
                if (found_ref.has_value())
                    break;
            }
        }
        if (!found_ref.has_value())
            found_ref = find_by_hi_name(*antigen);
        if (!found_ref.has_value())
            found_ref = find_by_parsed_name(*antigen);
        if (found_ref.has_value()) {
            result.refs_.push_back(std::move(*found_ref));
            ++num_matched;
        }
        else
            result.refs_.emplace_back();
    }
    if constexpr (std::is_same_v<AgSr, acmacs::chart::Antigens> || std::is_same_v<AgSr, acmacs::chart::AntigensModify>)
        AD_INFO("antigens from chart have sequences in seqdb: {}", num_matched);
    else
        AD_INFO("sera from chart have sequences in seqdb: {}", num_matched);

    return result;

} // acmacs::seqdb::v3::Seqdb::match

template acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::Antigens&, std::string_view) const;
template acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::AntigensModify&, std::string_view) const;
template acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::Sera&, std::string_view) const;
template acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::SeraModify&, std::string_view) const;

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

void acmacs::seqdb::v3::Seqdb::populate(acmacs::chart::ChartModify& chart) const
{
    const auto populate_ag_sr = [this, &chart]<typename AgSr>(AgSr& antigens_sera) {
        acmacs::enumerate(match(antigens_sera, chart.info()->virus_type(acmacs::chart::Info::Compute::Yes)), [&](auto no, const auto& ref) {
            if (ref) {
                const auto& seq = ref.seq().with_sequence(*this);
                auto& antigen_serum = antigens_sera.at(no);
                antigen_serum.sequence_aa(seq.aa_aligned_master());
                antigen_serum.sequence_nuc(seq.nuc_aligned_master());
                if (!seq.clades.empty()) {
                    for (const auto& clade : seq.clades)
                        antigen_serum.add_clade(std::string{clade});
                }
                else {
                    antigen_serum.add_clade("SEQUENCED");
                }
                if (const auto& lineage = ref.entry->lineage; !lineage.empty()) {
                    if (const auto ag_lineage = antigen_serum.lineage(); ag_lineage == acmacs::chart::BLineage::Unknown)
                        antigen_serum.lineage(lineage);
                    else if (ag_lineage != lineage) {
                        AD_WARNING("{} lineage difference, seqdb: {}, antigen_serum lineage in chart updated",
                                   acmacs::chart::format_antigen_serum<AgSr>("{ag_sr} {no0:{num_digits}d} {full_name} {lineage}", chart, no, acmacs::chart::collapse_spaces_t::yes), lineage);
                        antigen_serum.lineage(lineage);
                    }
                }
                AD_LOG(acmacs::log::hi_name_matching, "Seqdb::populate {} <-- {}",
                       acmacs::chart::format_antigen_serum<AgSr>("{ag_sr} {no0:{num_digits}d} {full_name}{ }{lineage}{ }{clades}", chart, no, acmacs::chart::collapse_spaces_t::yes), ref.seq_id());
            }
        });
    };

    populate_ag_sr(chart.antigens_modify());
    populate_ag_sr(chart.sera_modify());

} // acmacs::seqdb::v3::Seqdb::populate

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
            fasta += fmt::format(">{}\n{}\n", antigens->at(ag_no)->format("{name_full}"), ref.nuc_aligned(*this));
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
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
