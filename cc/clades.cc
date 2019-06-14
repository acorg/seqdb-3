#include "seqdb-3/clades.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

namespace local
{
    void detect_B_lineage(acmacs::seqdb::v3::sequence_t& sequence);
    void detect_B_clade(acmacs::seqdb::v3::sequence_t& sequence);
    void detect_H1_clade(acmacs::seqdb::v3::sequence_t& sequence);
    void detect_H3_clade(acmacs::seqdb::v3::sequence_t& sequence);
}

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::detect_lineages_clades(std::vector<fasta::scan_result_t>& sequences)
{
#pragma omp parallel for default(shared) schedule(static, 256)
    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        if (auto& entry = sequences[e_no]; fasta::is_aligned(entry)) {
            const auto subtype = entry.sequence.type_subtype().h_or_b();
            if (subtype == "B") {
                local::detect_B_lineage(entry.sequence);
                local::detect_B_clade(entry.sequence);
            }
            else if (subtype == "H1") {
                local::detect_H1_clade(entry.sequence);
            }
            else if (subtype == "H3") {
                local::detect_H3_clade(entry.sequence);
            }
        }
    }
        // detect B lineage and VIC deletion mutants, adjust deletions
        // detect clades

} // acmacs::seqdb::v3::detect_lineages_clades

// ----------------------------------------------------------------------

namespace local
{
    inline bool is_victoria(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.empty();
    }

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

    inline bool is_victoria_sixdel2019(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 163 && deletions.deletions.front().num == 6 && deletions.insertions.empty();
    }

    inline bool is_victoria_deletions_at_the_end(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos > 500 && deletions.insertions.empty();
    }

    inline bool is_yamagata_shifted(acmacs::seqdb::v3::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return ((deletions.size() == 1 && deletions.front().pos == 158 && deletions.front().num == 1 && sequence.aa_aligned_substr(155, 6) == "MAWVIP")
                || (deletions.size() == 1 && deletions.front().pos == 161 && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 2) == "VP")
                || (deletions.size() == 1 && deletions.front().pos == 160 && deletions.front().num == 1 && sequence.aa_aligned_substr(157, 3) == "WAV")
                || (deletions.size() == 1 && deletions.front().pos == 163 && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 3) == "VPK")
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

    inline bool is_yamagata_doubledel(acmacs::seqdb::v3::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return deletions.size() == 1 && deletions.front().pos == 162 && deletions.front().num == 2
                && sequence.year() <= 2013
                && sequence.deletions().insertions.empty()
                ;
    }

    // 12 sequences from TAIWAN 2010 have deletions 169:2
    inline bool is_taiwan_169_2(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 168 && deletions.deletions.front().num == 2 && deletions.insertions.empty();
    }

    inline bool is_semi_ignored(acmacs::seqdb::v3::sequence_t& sequence)
    {
        return *sequence.name() == "B/MIE/1/2019"         // DEL[1](162:4)<pos-1-based>  NIID:20190314   -- B/MIE/1/2019 |  2019-01-22 | MDCK 1 +1 |  18/19-498 | National Institute of Infectious Diseases (NIID) | B / H0N0 |  Victoria
                || *sequence.name() == "B/INDONESIA/NIHRDSB183950/2018" // DEL[1](164:2)<pos-1-based> VIDRL:20180913 -- B/Indonesia/NIHRDSB183950/2018 |  2018-04-01 | X, MDCK1 |  10004643 VW10005052 | WHO Collaborating Centre for Reference and Research on Influenza | B / H0N0 |  Victoria
                ;
    }

    inline bool is_ignored(acmacs::seqdb::v3::sequence_t& sequence)
    {
        return *sequence.name() == "B/ONTARIO/RV1769/2019"         // DEL[1](163:3)<pos-1-based>  B/Ontario/RV1769/2019 |  2019-04-11 | P1 |  RV1769/19 | Public Health Agency of Canada (PHAC) | B / H0N0 |  Victoria
                || *sequence.name() == "B/KENYA/4/2018"            // DEL[1](160:1)<pos-1-based>  B/Kenya/004/2018 |  2018-01-05 |  |   | Other Database Import | B / H0N0 |  unknown
                || *sequence.name() == "B/KENYA/11/2018"           // DEL[1](160:1)<pos-1-based>  B/Kenya/011/2018 |  2018-01-15 |  |   | Other Database Import | B / H0N0 |  unknown
                || *sequence.name() == "B/ORENBURG/CRIE-100/2018"  // DEL[1](160:1)<pos-1-based>  B/Orenburg/CRIE/100/2018 |  2018-02-08 |  |   | Central Research Institute of Epidemiology | B / H0N0 |  Yamagata
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
    if (is_victoria(deletions) || is_victoria_deletions_at_the_end(deletions)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("no");
    }
    else if (is_victoria_del2017(deletions)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria del2017");
        sequence.add_clade(acmacs::seqdb::v3::clade_t{"DEL2017"});
    }
    else if (is_victoria_tripledel2017(deletions)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria tripledel2017");
        sequence.add_clade(acmacs::seqdb::v3::clade_t{"TRIPLEDEL2017"});
    }
    else if (is_victoria_tripledel2017_pos_shifted_164(deletions)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria tripledel2017 (pos shifted)");
        deletions.deletions.front().pos = 161;
        sequence.add_clade(acmacs::seqdb::v3::clade_t{"TRIPLEDEL2017"});
    }
    else if (is_victoria_sixdel2019(deletions)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
            warn("victoria sixdel2019 (pos shifted)");
        sequence.add_clade(acmacs::seqdb::v3::clade_t{"SIXDEL2019"});
    }
    else if (is_yamagata_shifted(sequence)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
            warn("yamagata-shifted");
        deletions.deletions = std::vector<acmacs::seqdb::deletions_insertions_t::pos_num_t>{{162, 1}};
    }
    else if (is_yamagata(deletions)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
            warn("yamagata");
    }
    else if (is_yamagata_doubledel(sequence)) {
        if (sequence.lineage().empty())
            sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
        else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
            warn("yamagata");
    }
    else if (is_taiwan_169_2(deletions)) {
        // 12 sequences from TAIWAN 2010 have deletions 169:2
        sequence.lineage(acmacs::virus::lineage_t{});
        sequence.add_clade(acmacs::seqdb::v3::clade_t{"TAIWAN2010"});
    }
    else if (is_semi_ignored(sequence)) {
        fmt::print(stderr, "INFO: {} {}\n", sequence.full_name(), acmacs::seqdb::format(deletions));
    }
    else if (is_ignored(sequence)) {
        // do not issue warning
    }
    else {
        warn("unknown", "ERROR");
    }

} // local::detect_B_lineage

// ----------------------------------------------------------------------

void local::detect_B_clade(acmacs::seqdb::v3::sequence_t& sequence)
{
    if (sequence.lineage() == acmacs::virus::lineage_t{"VICTORIA"}) {
        // 2018-09-03, Sarah: clades should (technically) be defined by a phylogenetic tree rather than a set of amino acids
        if (sequence.aa_at_pos1(75) == 'K' && sequence.aa_at_pos1(172) == 'P' && sequence.aa_at_pos1(58) != 'P')
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"1A"});
        else if (sequence.aa_at_pos1(58) == 'P')
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"1B"});
        else
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"1"});
    }
    else if (sequence.lineage() == acmacs::virus::lineage_t{"YAMAGATA"}) {
        // 165N -> Y2, 165Y -> Y3 (yamagata numeration, 163 is not -)
        // 166N -> Y2, 166Y -> Y3 (victoria numeration, 163 is -)
        switch (sequence.aa_at_pos1(166)) {
            case 'N':
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"Y2"});
                break;
            case 'Y':
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"Y3"});
                break;
            default:
                break;
        }
    }

} // local::detect_B_clade

// ----------------------------------------------------------------------

void local::detect_H1_clade(acmacs::seqdb::v3::sequence_t& sequence)
{
    // // ----------------------------------------------------------------------
    // // 2018-09-19 clade definitions changed by Sarah before SSM
    // // ----------------------------------------------------------------------
    // // 6B: 163Q
    // // 6B1: 162N, 163Q
    // // 6B2: 152T, 163Q
    // auto r = std::vector<std::string>();
    // const auto pos152 = static_cast<size_t>(151 - aShift),
    //         pos162 = static_cast<size_t>(161 - aShift),
    //         pos163 = static_cast<size_t>(162 - aShift);
    // if (pos163 > 0 && aSequence.size() > pos163) {
    //     if (aSequence[pos163] == 'Q') {
    //         r.push_back("6B");
    //         if (aSequence[pos162] == 'N')
    //             r.push_back("6B1");
    //         if (aSequence[pos152] == 'T')
    //             r.push_back("6B2");
    //     }
    // }
    // return r;

} // local::detect_H1_clade

// ----------------------------------------------------------------------

void local::detect_H3_clade(acmacs::seqdb::v3::sequence_t& sequence)
{

} // local::detect_H3_clade

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
