#include "seqdb-3/clades.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

namespace local
{
    void detect_B_lineage(acmacs::seqdb::v3::sequence_t& sequence);
}

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::detect_lineages_clades(std::vector<fasta::scan_result_t>& sequences)
{
#pragma omp parallel for default(shared) schedule(static, 256)
    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        if (auto& entry = sequences[e_no]; fasta::is_aligned(entry)) {
            const auto& type_subtype = entry.sequence.type_subtype();
            if (*type_subtype == "B") {
                local::detect_B_lineage(entry.sequence);
            }
        }
    }
        // detect B lineage and VIC deletion mutants, adjust deletions
        // detect clades

} // acmacs::seqdb::v3::detect_lineages_clades

// ----------------------------------------------------------------------

namespace local
{
    inline bool is_victoria_del2017(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 161 && deletions.deletions.front().num == 2 && deletions.insertions.empty();
    }

    inline bool is_victoria_tripledel2017(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 161 && deletions.deletions.front().num == 3 && deletions.insertions.empty();
    }

    inline bool is_victoria_tripledel2017_pos_shifted_164(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 163 && deletions.deletions.front().num == 3 && deletions.insertions.empty();
    }

    inline bool is_yamagata_shifted(acmacs::seqdb::v3::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return ((deletions.size() == 1 && deletions.front().pos == 158 && deletions.front().num == 1 && sequence.aa_aligned_substr(155, 6) == "MAWVIP")
                || (deletions.size() == 1 && deletions.front().pos == 161 && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 2) == "VP")
                || (deletions.size() == 1 && deletions.front().pos == 163 && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 5) == "VPKXN")
                )
                && sequence.deletions().insertions.empty();
    }

    inline bool is_yamagata(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return !deletions.deletions.empty() && deletions.deletions.front().pos == 162 && deletions.deletions.front().num == 1
                && (deletions.deletions.size() == 1 || deletions.deletions[1].pos > 500)
                && deletions.insertions.empty()
                ;
    }

}

// B/Yamagata/16/88
// B/Victoria/2/87

// VICTORIA del2017: 162, 163
// VICTORIA tripledel2017: 162, 163, 164 by convention

// YAMAGATA: deletion must be at 163
// David Burke 2017-08-17: deletions ( and insertions) of amino acids usually occur in regions of the protein structure where it changes direction ( loops ).
// In the case of HA, this is after VPK and before NKTAT/YKNAT.

void local::detect_B_lineage(acmacs::seqdb::v3::sequence_t& sequence)
{
    const auto warn = [&](const char* infix, const char* prefix = "WARNING") {
        fmt::print(stderr, "{}: {} lineage {} and {} deletions {} {}\n{}\n", prefix, sequence.year(), *sequence.lineage(), infix, sequence.full_name(),
                   acmacs::seqdb::format(sequence.deletions()), acmacs::seqdb::format(sequence.deletions().deletions, sequence.aa_aligned()));
    };

    auto& deletions = sequence.deletions();
    if (deletions.empty()) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("no");
    }
    else if (is_victoria_del2017(sequence.deletions())) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria del2017");
        // sequence.add_clade("DEL2017");
    }
    else if (is_victoria_tripledel2017(sequence.deletions())) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria tripledel2017");
        // sequence.add_clade("TRIPLEDEL2017");
    }
    else if (is_victoria_tripledel2017_pos_shifted_164(sequence.deletions())) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria tripledel2017 (pos shifted)");
        sequence.deletions().deletions.front().pos = 161;
        // sequence.add_clade("TRIPLEDEL2017");
    }
    else if (is_yamagata_shifted(sequence)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
            warn("yamagata-shifted");
        sequence.deletions().deletions = std::vector<acmacs::seqdb::deletions_insertions_t::pos_num_t>{{162, 1}};
    }
    else if (is_yamagata(sequence.deletions())) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
            warn("yamagata");
    }
    else {
        warn("unknown", "ERROR");
    }

} // local::detect_B_lineage

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
