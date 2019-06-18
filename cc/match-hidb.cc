#include "hidb-5/hidb.hh"
#include "seqdb-3/match-hidb.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

static void find_by_lab_id(hidb::AntigenPList& found, const hidb::Antigens& hidb_antigens, const acmacs::seqdb::v3::sequence_t& sequence);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::match_hidb(std::vector<fasta::scan_result_t>& sequences)
{
    for (auto& en : sequences) {
        auto& seq = en.sequence;
        if (const auto hb = seq.type_subtype().h_or_b(); hb == "B" || hb == "H1" || hb == "H3") {
            auto hidb_antigens = hidb::get(hb).antigens();
            hidb::AntigenPList found_hidb_antigens;
            find_by_lab_id(found_hidb_antigens, *hidb_antigens, seq);
        }
    }

} // acmacs::seqdb::v3::match_hidb

// ----------------------------------------------------------------------

void find_by_lab_id(hidb::AntigenPList& found, const hidb::Antigens& hidb_antigens, const acmacs::seqdb::v3::sequence_t& sequence)
{
    if (sequence.lab() == "CDC" && !sequence.lab_id().empty()) {
        // const auto f_cdcid = hidb_antigens.find_labid(sequence.lab_id());
        // std::copy(f_cdcid.begin(), f_cdcid.end(), std::back_inserter(found));
    }

} // find_by_lab_id

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
