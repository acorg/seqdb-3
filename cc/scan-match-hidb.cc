#include "acmacs-base/string-matcher.hh"
#include "hidb-5/hidb.hh"
#include "seqdb-3/scan-match-hidb.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
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

} // namespace acmacs::seqdb::inline v3

using lab_id_index_entry_t = std::pair<std::string_view, const hidb::bin::Antigen*>;

struct hidb_ref_t
{
    const hidb::HiDb& hidb;
    std::shared_ptr<hidb::Antigens> antigens;
    std::vector<lab_id_index_entry_t> lab_id_index;
};

using Matching = std::vector<std::vector<acmacs::seqdb::v3::score_seq_found_t>>;

using seq_iter_t = std::vector<acmacs::seqdb::scan::fasta::scan_result_t>::iterator;

using seq_ptr_list_t = std::vector<acmacs::seqdb::scan::sequence_t*>;
using hi_to_seq_t = std::map<std::pair<const hidb_ref_t*, hidb::AntigenIndex>, seq_ptr_list_t>;

static void match(const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last, std::string_view subtype, hi_to_seq_t& hi_to_seq);
static void find_by_lab_id(hidb::AntigenIndexList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last);
static void find_by_name(hidb::AntigenIndexList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last);
static Matching make_matching(seq_iter_t first, seq_iter_t last, const hidb::AntigenPList& found);
static void match_greedy(seq_iter_t first, const hidb::AntigenIndexList& found, const hidb::AntigenPList& antigens, const Matching& matching, const hidb_ref_t& hidb_ref, hi_to_seq_t& hi_to_seq);
// static bool match_normal(seq_iter_t first, const hidb::AntigenPList& found, const Matching& matching);
static void update_seqdb(hi_to_seq_t& hi_to_seq);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::match_hidb(std::vector<fasta::scan_result_t>& sequences)
{
    AD_INFO("INFO: matching against hidb");
    // sequences must be sorted by name!
    std::map<std::string, hidb_ref_t, std::less<>> hidbs;
    for (const std::string_view subtype : {"B", "H1", "H3"}) {
        const auto& hidb = hidb::get(acmacs::virus::type_subtype_t{subtype}, report_time::no);
        auto hidb_antigens = hidb.antigens();
        hidbs.emplace(std::string{subtype}, hidb_ref_t{hidb, hidb_antigens, hidb_antigens->sorted_by_labid()});
    }

    hi_to_seq_t hi_to_seq;
    for (auto en_first = sequences.begin(); en_first != sequences.end();) {
        const auto en_last = std::find_if(std::next(en_first), sequences.end(), [name = en_first->sequence.name()](const auto& en) { return en.sequence.name() != name; });
        if (const auto hb = en_first->sequence.type_subtype().h_or_b(); hb == "B" || hb == "H1" || hb == "H3")
            match(hidbs.find(hb)->second, en_first, en_last, hb, hi_to_seq);
        en_first = en_last;
    }

    update_seqdb(hi_to_seq);

    AD_INFO("INFO: matched against hidb: {}", hi_to_seq.size());

} // acmacs::seqdb::v3::scan::match_hidb

// ----------------------------------------------------------------------

void update_seqdb(hi_to_seq_t& hi_to_seq)
{
    const auto update = [](const auto& hi, acmacs::seqdb::v3::scan::sequence_t& seq) {
        auto antigen = hi.first->antigens->at(hi.second);
        const auto name = antigen->full_name();
        seq.add_hi_name(name);
        // // AD_DEBUG("    add {} -> {}", name, sequence.full_name()); // std::next(first, static_cast<ssize_t>(en.second.seq_no))->fasta.name);
        if (const size_t subtype_size = name.find('/'); subtype_size > 1 && subtype_size <= 8)
            seq.update_subtype(acmacs::virus::type_subtype_t{name.substr(0, subtype_size)});
        if (const auto& date = antigen->date(); !date.empty())
            seq.add_date(*date);
    };

    // if a hi-name can be used with multiple sequences, choose just
    // one sequence, the one with closest passage and preferably with
    // the same lab.
    for (auto& [ag, sequences] : hi_to_seq) {
        if (sequences.size() == 1) {
            update(ag, *sequences[0]);
        }
        else {
            auto antigen = ag.first->antigens->at(ag.second);
            const acmacs::uppercase hi_lab{ag.first->hidb.lab(*antigen)};
            // AD_DEBUG("    hidb {} lab: {} passage: \"{}\" ({})", antigen->full_name(), hi_lab, antigen->passage(), sequences.size());
            std::sort(std::begin(sequences), std::end(sequences), [&antigen, &hi_lab](auto* s1, auto* s2) {
                auto score1 = acmacs::virus::passage_compare(s1->passage(), antigen->passage()) * 10 + (s1->lab_in({hi_lab}) ? 0 : 1);
                auto score2 = acmacs::virus::passage_compare(s2->passage(), antigen->passage()) * 10 + (s2->lab_in({hi_lab}) ? 0 : 1);
                return score1 < score2;
            });
            update(ag, *sequences.front());
            // for (auto* seq : sequences) {
            //     update(ag, *seq);
            //     std::string seq_lab;
            //     if (!seq->lab_ids().empty())
            //         seq_lab = seq->lab_ids().begin()->first;
            //     AD_DEBUG("        seqdb {} lab: {} passage:\"{}\"", seq->full_name(), seq_lab, seq->passage());
            // }
        }
    }
}

// ----------------------------------------------------------------------

void match(const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last, std::string_view subtype, hi_to_seq_t& hi_to_seq)
{
    hidb::AntigenIndexList found_hidb_antigens;
    find_by_lab_id(found_hidb_antigens, hidb_ref, first, last);
    find_by_name(found_hidb_antigens, hidb_ref, first, last);
    std::sort(found_hidb_antigens.begin(), found_hidb_antigens.end());
    found_hidb_antigens.erase(std::unique(found_hidb_antigens.begin(), found_hidb_antigens.end()), found_hidb_antigens.end());

    if (!found_hidb_antigens.empty()) {
        const hidb::AntigenPList antigens = hidb_ref.antigens->list(found_hidb_antigens);
        const auto& seq = first->sequence;
        if (subtype == "B") {
            if (const auto hidb_lineage = antigens[0]->lineage(); hidb_lineage != acmacs::chart::BLineage::Unknown && hidb_lineage != acmacs::chart::BLineage{seq.lineage()})
                AD_WARNING("lineage mismatch seq: {} vs. hidb: {} {}", seq.full_name(), antigens[0]->name(), antigens[0]->lineage());
        }
        for (const auto& antigen: antigens)
            first->sequence.add_date(*antigen->date());

        const auto matching = make_matching(first, last, antigens); // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
        match_greedy(first, found_hidb_antigens, antigens, matching, hidb_ref, hi_to_seq);
        // matched = match_normal(first, found, matching);
    }

} // match

// ----------------------------------------------------------------------

Matching make_matching(seq_iter_t first, seq_iter_t last, const hidb::AntigenPList& found)
{
    Matching matching;
    size_t seq_no = 0;
    for (; first != last; ++first) {
        const auto& seq = first->sequence;
        // fmt::print(stderr, "DEBUG: seq passages: {}\n", seq.passages());
        std::vector<acmacs::seqdb::v3::score_seq_found_t> matching_for_seq;
        size_t found_no = 0;
        for (auto f : found) {
            if (seq.reassortant() == f->reassortant()) {
                const auto f_passage = f->passage();
                for (const auto& s_passage : seq.passages()) {
                    if (acmacs::virus::passages_match(f_passage, s_passage))
                        matching_for_seq.emplace_back(acmacs::seqdb::v3::score_size_t{string_match::match(*s_passage, *f_passage), std::min(s_passage.size(), f_passage.size())}, seq_no, found_no);
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
void match_greedy(seq_iter_t first, const hidb::AntigenIndexList& found, const hidb::AntigenPList& antigens, const Matching& matching, const hidb_ref_t& hidb_ref, hi_to_seq_t& hi_to_seq)
{
    std::map<size_t, acmacs::seqdb::v3::score_seq_found_t> antigen_to_matching; // antigen index in found to (matching index and score)
    for (const auto& mp : matching) {
        for (const auto& sf: mp) {
            // fmt::print(stderr, "DEBUG: hidb:{} score:{} seqdb:{}\n", found[sf.found_no]->full_name(), sf.score, sf.seq_no);
            const auto ampi = antigen_to_matching.emplace(sf.found_no, sf);
            if (!ampi.second && ampi.first->second.score < sf.score) {        // already present, replace if sf has hi higher score
                ampi.first->second = sf;
            }
        }
    }

    // AD_DEBUG("match_greedy");
    for (const auto& en : antigen_to_matching) {
        auto& sequence = std::next(first, static_cast<ssize_t>(en.second.seq_no))->sequence;
        hi_to_seq.try_emplace(std::make_pair(&hidb_ref, found[en.first]), seq_ptr_list_t{}).first->second.push_back(&sequence);
    }

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

void find_by_lab_id(hidb::AntigenIndexList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last)
{
    for (; first != last; ++first) {
        const auto& sequence = first->sequence;
        if (const auto cdcids = sequence.lab_ids().find(acmacs::uppercase{"CDC"}); cdcids != sequence.lab_ids().end() && !cdcids->second.empty()) {
            for (const auto& cdcid_raw : cdcids->second) {
                const auto cdcid = fmt::format("CDC#{}", cdcid_raw);
                auto [hidb_first, hidb_last] = std::equal_range(std::begin(hidb_ref.lab_id_index), std::end(hidb_ref.lab_id_index), lab_id_index_entry_t{std::string_view(cdcid), nullptr},
                                                                [](const auto& e1, const auto& e2) { return e1.first < e2.first; });
                for (; hidb_first != hidb_last; ++hidb_first)
                    found.push_back(hidb_ref.antigens->index(hidb_first->second));
            }
        }
    }

} // find_by_lab_id

// ----------------------------------------------------------------------

void find_by_name(hidb::AntigenIndexList& found, const hidb_ref_t& hidb_ref, seq_iter_t first, seq_iter_t last)
{
    for (; first != last; ++first) {
        const auto& sequence = first->sequence;
        // AD_DEBUG("find antigen in hidb: \"{}\"", *sequence.name());
        const auto antigen_index_list = hidb_ref.antigens->find(*sequence.name(), hidb::fix_location::no, hidb::find_fuzzy::no);
        std::copy(antigen_index_list.begin(), antigen_index_list.end(), std::back_inserter(found));
    }

} // find_by_name

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
