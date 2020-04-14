#include "acmacs-base/string-matcher.hh"
#include "hidb-5/hidb.hh"
#include "seqdb-3/scan-match-hidb.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

using lab_id_index_entry_t = std::pair<std::string_view, const hidb::bin::Antigen*>;

struct hidb_ref_t
{
    std::shared_ptr<hidb::Antigens> antigens;
    std::vector<lab_id_index_entry_t> lab_id_index;
};

struct score_size_t
{
    string_match::score_t score;
    size_t len;

    constexpr bool operator<(const score_size_t& a) const { return score < a.score; }
};

struct score_seq_found_t : public score_size_t
{
    size_t seq_no;
    size_t found_no;

    score_seq_found_t(const score_size_t& ss, size_t sn, size_t fn) : score_size_t{ss.score, ss.len}, seq_no{sn}, found_no{fn} {}
    constexpr bool operator<(const score_seq_found_t& a) const { return score > a.score; }
};

using Matching = std::vector<std::vector<score_seq_found_t>>;

using seq_iter_t = std::vector<acmacs::seqdb::scan::fasta::scan_result_t>::iterator;

static bool match(const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last, std::string_view subtype);
static void find_by_lab_id(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last);
static void find_by_name(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last);
static Matching make_matching(seq_iter_t first, seq_iter_t last, const hidb::AntigenPList& found);
static bool match_greedy(seq_iter_t first, const hidb::AntigenPList& found, const Matching& matching);
// static bool match_normal(seq_iter_t first, const hidb::AntigenPList& found, const Matching& matching);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::match_hidb(std::vector<fasta::scan_result_t>& sequences)
{
    AD_INFO("INFO: matching against hidb");
        // sequences msut be sorted by name!
    std::map<std::string, hidb_ref_t, std::less<>> hidbs;
    for (const std::string_view subtype : {"B", "H1", "H3"}) {
        auto hidb_antigens = hidb::get(acmacs::virus::type_subtype_t{subtype}, report_time::no).antigens();
        hidbs.emplace(std::string{subtype}, hidb_ref_t{hidb_antigens, hidb_antigens->sorted_by_labid()});
    }

    size_t matched = 0;
    for (auto en_first = sequences.begin(); en_first != sequences.end(); ) {
        const auto en_last = std::find_if(std::next(en_first), sequences.end(), [name=en_first->sequence.name()](const auto& en) { return en.sequence.name() != name; });
        if (const auto hb = en_first->sequence.type_subtype().h_or_b(); hb == "B" || hb == "H1" || hb == "H3") {
            if (match(hidbs.find(hb)->second, en_first, en_last, hb))
                ++matched;
        }
        en_first = en_last;
    }
    AD_INFO("INFO: matched against hidb: {}", matched);

} // acmacs::seqdb::v3::scan::match_hidb

// ----------------------------------------------------------------------

bool match(const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last, std::string_view subtype)
{
    hidb::AntigenPList found_hidb_antigens;
    find_by_lab_id(found_hidb_antigens, hidb_ref, first, last);
    find_by_name(found_hidb_antigens, hidb_ref, first, last);
    std::sort(found_hidb_antigens.begin(), found_hidb_antigens.end());
    found_hidb_antigens.erase(std::unique(found_hidb_antigens.begin(), found_hidb_antigens.end()), found_hidb_antigens.end());

    bool matched = false;
    if (!found_hidb_antigens.empty()) {
        const auto& seq = first->sequence;
        if (subtype == "B") {
            if (const auto hidb_lineage = found_hidb_antigens.front()->lineage(); hidb_lineage != acmacs::chart::BLineage::Unknown && hidb_lineage != acmacs::chart::BLineage{seq.lineage()})
                AD_WARNING("lineage mismatch seq: {} vs. hidb: {} {}", seq.full_name(), found_hidb_antigens.front()->name(), found_hidb_antigens.front()->lineage());
        }
        for (const auto& en: found_hidb_antigens)
            first->sequence.add_date(*en->date());

        const auto matching = make_matching(first, last, found_hidb_antigens); // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
        matched = match_greedy(first, found_hidb_antigens, matching);
        // matched = match_normal(first, found, matching);
    }
    return matched;

} // match

// ----------------------------------------------------------------------

Matching make_matching(seq_iter_t first, seq_iter_t last, const hidb::AntigenPList& found)
{
    Matching matching;
    size_t seq_no = 0;
    for (; first != last; ++first) {
        const auto& seq = first->sequence;
        // fmt::print(stderr, "DEBUG: seq passages: {}\n", seq.passages());
        std::vector<score_seq_found_t> matching_for_seq;
        size_t found_no = 0;
        for (auto f : found) {
            if (seq.reassortant() == f->reassortant()) {
                const auto f_passage = f->passage();
                for (const auto& s_passage : seq.passages()) {
                    if (acmacs::virus::passages_match(f_passage, s_passage))
                        matching_for_seq.emplace_back(score_size_t{string_match::match(*s_passage, *f_passage), std::min(s_passage.size(), f_passage.size())}, seq_no, found_no);
                }
            }
            ++found_no;
        }
        std::sort(matching_for_seq.begin(), matching_for_seq.end());
        matching.push_back(std::move(matching_for_seq));
        ++seq_no;
    }
    std::sort(matching.begin(), matching.end(), [](const auto& a, const auto& b) -> bool { return a.empty() ? false : (b.empty() || a[0] < b[0]); });
    return matching;

} // make_matching

// ----------------------------------------------------------------------

// greedy matching: add all hi-names having matching reassortant and passage type (egg/cell) regardless of score
// if antigen is in multiple matching entries, use the one with the highest score
// returns if at least one seq matched
bool match_greedy(seq_iter_t first, const hidb::AntigenPList& found, const Matching& matching)
{
    std::map<size_t, score_seq_found_t> antigen_to_matching; // antigen index in found to (matching index and score)
    for (const auto& mp : matching) {
        for (const auto& sf: mp) {
            // fmt::print(stderr, "DEBUG: hidb:{} score:{} seqdb:{}\n", found[sf.found_no]->full_name(), sf.score, sf.seq_no);
            const auto ampi = antigen_to_matching.emplace(sf.found_no, sf);
            if (!ampi.second && ampi.first->second.score < sf.score) {        // already present, replace if sf has hi higher score
                ampi.first->second = sf;
            }
        }
    }

    bool matched = false;
    for (const auto& e: antigen_to_matching) {
        const auto name = found[e.first]->full_name();
        auto& sequence = std::next(first, static_cast<ssize_t>(e.second.seq_no))->sequence;
        sequence.add_hi_name(name);
        // fmt::print(stderr, "DEBUG: add {} -> {}\n", name, sequence.full_name()); // std::next(first, static_cast<ssize_t>(e.second.seq_no))->fasta.name);
        if (const size_t subtype_size = name.find('/'); subtype_size > 1 && subtype_size <= 8)
            sequence.update_subtype(acmacs::virus::type_subtype_t{name.substr(0, subtype_size)});
        if (const auto& date = found[e.first]->date(); !date.empty())
            first->sequence.add_date(*date);
        matched = true;
    }
    return matched;

} // match_greedy

// ----------------------------------------------------------------------

// bool match_normal(seq_iter_t first, const hidb::AntigenPList& found, const Matching& matching)
// {
//     bool matched = false;
//     if (matching.size() == 1) {
//         for (const auto& ms: matching[0]) {
//             if (ms.score == matching[0][0].score) {
//                 const auto name = found[ms.found_no]->full_name();
//                 first->sequence.add_hi_name(name);
//                 matched = true;
//             }
//         }
//     }
//     else {
//         std::set<size_t> found_assigned;
//         for (const auto& m: matching) {
//             for (const auto& sf: m) {
//                 if (sf.score == m[0].score && found_assigned.count(sf.found_no) == 0) {
//                     const auto name = found[sf.found_no]->full_name();
//                     std::next(first, static_cast<ssize_t>(sf.seq_no))->sequence.add_hi_name(name);
//                     matched = true;
//                     found_assigned.insert(sf.found_no);
//                 }
//             }
//         }
//     }
//     return matched;

// } // match_normal

// ----------------------------------------------------------------------

void find_by_lab_id(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last)
{
    for (; first != last; ++first) {
        const auto& sequence = first->sequence;
        if (const auto cdcids = sequence.lab_ids().find(acmacs::uppercase{"CDC"}); cdcids != sequence.lab_ids().end() && !cdcids->second.empty()) {
            for (const auto& cdcid_raw : cdcids->second) {
                const auto cdcid = fmt::format("CDC#{}", cdcid_raw);
                auto [hidb_first, hidb_last] = std::equal_range(std::begin(hidb_ref.lab_id_index), std::end(hidb_ref.lab_id_index), lab_id_index_entry_t{std::string_view(cdcid), nullptr},
                                                                [](const auto& e1, const auto& e2) { return e1.first < e2.first; });
                // fmt::print(stderr, "find_by_lab_id {} {}\n", cdcid, last - first);
                for (; hidb_first != hidb_last; ++hidb_first)
                    found.push_back(hidb_ref.antigens->make(hidb_first->second));
            }
        }
    }

} // find_by_lab_id

// ----------------------------------------------------------------------

void find_by_name(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last)
{
    for (; first != last; ++first) {
        const auto& sequence = first->sequence;
        // AD_DEBUG("find antigen in hidb: \"{}\"", *sequence.name());
        const auto antigen_index_list = hidb_ref.antigens->find(*sequence.name(), hidb::fix_location::no, hidb::find_fuzzy::no);
        std::transform(antigen_index_list.begin(), antigen_index_list.end(), std::back_inserter(found), [](const hidb::AntigenPIndex& antigen_index) -> hidb::AntigenP { return antigen_index.first; });
    }

} // find_by_name

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
