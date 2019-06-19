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

static void find_by_lab_id(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& sequence);
static void find_by_name(hidb::AntigenPList& found, const hidb_ref_t& hidb_ref, const acmacs::seqdb::v3::sequence_t& sequence);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::match_hidb(std::vector<fasta::scan_result_t>& sequences)
{
    std::map<std::string, hidb_ref_t, std::less<>> hidbs;
    for (const std::string_view subtype : {"B", "H1", "H3"}) {
        auto hidb_antigens = hidb::get(subtype).antigens();
        hidbs.emplace(std::string{subtype}, hidb_ref_t{hidb_antigens, hidb_antigens->sorted_by_labid()});
    }

    for (auto& en : sequences) {
        auto& seq = en.sequence;
        if (const auto hb = seq.type_subtype().h_or_b(); hb == "B" || hb == "H1" || hb == "H3") {
            const auto& hidb_ref = hidbs.find(hb)->second;
            hidb::AntigenPList found_hidb_antigens;
            find_by_lab_id(found_hidb_antigens, hidb_ref, seq);
            // find_by_name(found_hidb_antigens, hidb_ref, seq);
            std::sort(found_hidb_antigens.begin(), found_hidb_antigens.end());
            found_hidb_antigens.erase(std::unique(found_hidb_antigens.begin(), found_hidb_antigens.end()), found_hidb_antigens.end());
            if (!found_hidb_antigens.empty()) {
                // update country, continent?
                // check lineage
                if (hb == "B") {
                    if (const auto hidb_lineage = found_hidb_antigens.front()->lineage(); hidb_lineage != acmacs::chart::BLineage::Unknown && hidb_lineage != *seq.lineage())
                        fmt::print(stderr, "WARNING: lineage mismatch seq: {} vs. hidb: {} {}\n", seq.full_name(), found_hidb_antigens.front()->name(), static_cast<std::string>(found_hidb_antigens.front()->lineage()));
                }
                // update date
                // for (const auto& en: found_hidb_antigens)
                //     seq.add_date(en->date());
            }
        }
    }

} // acmacs::seqdb::v3::match_hidb

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
