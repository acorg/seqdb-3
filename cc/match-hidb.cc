#include "acmacs-base/string-matcher.hh"
#include "hidb-5/hidb.hh"
#include "seqdb-3/match-hidb.hh"
#include "seqdb-3/fasta.hh"

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

static bool match(const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& seq, std::string_view subtype);
static void find_by_lab_id(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& sequence);
static void find_by_name(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& sequence);
static Matching make_matching(const acmacs::seqdb::v3::sequence_t& seq, const hidb::AntigenPList& found);
static bool match_greedy(const acmacs::seqdb::v3::sequence_t& seq, const hidb::AntigenPList& found, const Matching& matching);
static bool match_normal(const acmacs::seqdb::v3::sequence_t& seq, const hidb::AntigenPList& found, const Matching& matching);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::match_hidb(std::vector<fasta::scan_result_t>& sequences)
{
    std::map<std::string, hidb_ref_t, std::less<>> hidbs;
    for (const std::string_view subtype : {"B", "H1", "H3"}) {
        auto hidb_antigens = hidb::get(subtype).antigens();
        hidbs.emplace(std::string{subtype}, hidb_ref_t{hidb_antigens, hidb_antigens->sorted_by_labid()});
    }

    size_t matched = 0;
    for (auto& en : sequences) {
        auto& seq = en.sequence;
        if (const auto hb = seq.type_subtype().h_or_b(); hb == "B" || hb == "H1" || hb == "H3") {
            if (match(hidbs.find(hb)->second, seq, hb))
                ++matched;
        }
    }
    fmt::print("INFO: matched against hidb: {}\n", matched);

} // acmacs::seqdb::v3::match_hidb

// ----------------------------------------------------------------------

bool match(const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& seq, std::string_view subtype)
{
    hidb::AntigenPList found_hidb_antigens;
    find_by_lab_id(found_hidb_antigens, hidb_ref, seq);
    find_by_name(found_hidb_antigens, hidb_ref, seq);
    std::sort(found_hidb_antigens.begin(), found_hidb_antigens.end());
    found_hidb_antigens.erase(std::unique(found_hidb_antigens.begin(), found_hidb_antigens.end()), found_hidb_antigens.end());

    bool matched = false;
    if (!found_hidb_antigens.empty()) {
        // update country, continent?
        if (subtype == "B") {
            if (const auto hidb_lineage = found_hidb_antigens.front()->lineage(); hidb_lineage != acmacs::chart::BLineage::Unknown && hidb_lineage != *seq.lineage())
                fmt::print(stderr, "WARNING: lineage mismatch seq: {} vs. hidb: {} {}\n", seq.full_name(), found_hidb_antigens.front()->name(),
                           static_cast<std::string>(found_hidb_antigens.front()->lineage()));
        }
        // update date
        // for (const auto& en: found_hidb_antigens)
        //     seq.add_date(en->date());

        const auto matching = make_matching(seq, found_hidb_antigens); // for each seq list of matching [[score, min passage len], found_no] - sorted by score desc
        matched = match_greedy(seq, found_hidb_antigens, matching);
        // matched = match_normal(seq, found, matching);
    }
    return matched;

} // match

// ----------------------------------------------------------------------

Matching make_matching(const acmacs::seqdb::v3::sequence_t& seq, const hidb::AntigenPList& found)
{
    // enum class CellOrEgg {Unknown, Cell, Egg, CellAndEgg};

    // const auto cell_or_egg = [](std::string aPassage) -> CellOrEgg
    // {
    //     acmacs::virus::Passage passage(aPassage);
    //     return passage.is_egg() ? CellOrEgg::Egg : (passage.is_cell() ? CellOrEgg::Cell : CellOrEgg::CellAndEgg); // OR is CellAndEgg
    // };

    // const auto cell_or_egg_v = [cell_or_egg](const std::vector<std::string>& variants) -> CellOrEgg
    // {
    //     CellOrEgg r = CellOrEgg::Unknown;
    //     for (const auto& passage: variants) {
    //         CellOrEgg p = cell_or_egg(passage);
    //         if (r == CellOrEgg::Unknown)
    //             r = p;
    //         else if (p != r) {
    //             r = CellOrEgg::CellAndEgg;
    //             break;
    //         }
    //     }
    //     return r == CellOrEgg::Unknown ? CellOrEgg::CellAndEgg : r; // unknown (or no passage) gives CellOrEgg::CellAndEgg to match anything
    // };

    // const auto match_cell_egg = [](CellOrEgg a, CellOrEgg b) -> bool
    // {
    //     return a == b || a == CellOrEgg::CellAndEgg || b == CellOrEgg::CellAndEgg;

    // };

    Matching matching;
    // size_t seq_no = 0;
    // for (auto& seq: entry.seqs()) {
    //     std::vector<score_seq_found_t> matching_for_seq;
    //     size_t found_no = 0;
    //     const auto seq_cell_or_egg = cell_or_egg_v(seq.passages());
    //     for (auto f: found) {
    //         const auto f_passage = f->passage();
    //           // std::cerr << "match_cell_egg: " << match_cell_egg(cell_or_egg(f_passage), seq_cell_or_egg) << " -- " << f_passage << ':' << static_cast<int>(passage::cell_or_egg(f_passage)) << " " << seq.passages() << ':' << static_cast<int>(seq_cell_or_egg) << '\n';
    //         if (seq.reassortant_match(f->reassortant()) && match_cell_egg(cell_or_egg(f_passage), seq_cell_or_egg)) {
    //             std::vector<score_size_t> scores; // score and min passage length (to avoid too incomplete matches)
    //             if (!seq.passages().empty())
    //                 std::transform(seq.passages().begin(), seq.passages().end(), std::back_inserter(scores),
    //                                [&f_passage](const auto& passage) -> score_size_t { return {string_match::match(passage, f_passage), std::min(passage.size(), f_passage.size())}; });
    //             else
    //                 scores.emplace_back(string_match::match(std::string{}, f_passage), 0);
    //             matching_for_seq.emplace_back(*std::max_element(scores.begin(), scores.end() /*$, [](const auto& a, const auto& b) { return a.first < b.first; }*/), seq_no, found_no);
    //               // report_stream << "  @" << seq.passages() << " @ " << f_passage << " " << score_size->first << " " << score_size->second << '\n';
    //         }
    //         ++found_no;
    //     }
    //     std::sort(matching_for_seq.begin(), matching_for_seq.end());
    //     matching.push_back(std::move(matching_for_seq));
    //     ++seq_no;
    // }
    // std::sort(matching.begin(), matching.end(), [](const auto& a, const auto& b) -> bool { return a.empty() ? false : (b.empty() || a[0] < b[0]); });
    return matching;

} // make_matching

// ----------------------------------------------------------------------

bool match_greedy(const acmacs::seqdb::v3::sequence_t& seq, const hidb::AntigenPList& found, const Matching& matching)
{
    return false;
} // match_greedy

// ----------------------------------------------------------------------

bool match_normal(const acmacs::seqdb::v3::sequence_t& seq, const hidb::AntigenPList& found, const Matching& matching)
{
    return false;

} // match_normal

// ----------------------------------------------------------------------

void find_by_lab_id(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& sequence)
{
    if (sequence.lab() == "CDC" && !sequence.lab_id().empty()) {
        const auto cdcid = fmt::format("CDC#{}", sequence.lab_id());
        auto [first, last] = std::equal_range(std::begin(hidb_ref.lab_id_index), std::end(hidb_ref.lab_id_index), lab_id_index_entry_t{std::string_view(cdcid), nullptr},
                                              [](const auto& e1, const auto& e2) { return e1.first < e2.first; });
        // fmt::print(stderr, "find_by_lab_id {} {}\n", cdcid, last - first);
        for (; first != last; ++first)
            found.push_back(hidb_ref.antigens->make(first->second));
    }

} // find_by_lab_id

// ----------------------------------------------------------------------

void find_by_name(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& sequence)
{
    const auto antigen_index_list = hidb_ref.antigens->find(*sequence.name(), hidb::fix_location::no, hidb::find_fuzzy::no);
    std::transform(antigen_index_list.begin(), antigen_index_list.end(), std::back_inserter(found), [](const hidb::AntigenPIndex& antigen_index) -> hidb::AntigenP { return antigen_index.first; });

} // find_by_name

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
