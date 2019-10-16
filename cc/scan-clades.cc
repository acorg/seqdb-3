#include <array>

#include "seqdb-3/scan-clades.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

namespace local
{
    namespace B
    {
        static void lineage(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
        static void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
    } // namespace B

    namespace H1
    {
        static void deletions(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
        static void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
    } // namespace H1

    namespace H3
    {
        static void deletions(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
        static void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
    } // namespace H3

} // namespace local

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::detect_lineages_clades(std::vector<fasta::scan_result_t>& sequences)
{
#pragma omp parallel for default(shared) schedule(static, 256)
    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        if (auto& entry = sequences[e_no]; fasta::is_aligned(entry)) {
            const auto subtype = entry.sequence.type_subtype().h_or_b();
            const auto fasta_ref = fmt::format("{}:{}: note:  {}", entry.fasta.filename, entry.fasta.line_no, entry.fasta.entry_name);
            if (subtype == "B") {
                local::B::lineage(entry.sequence, fasta_ref);
                local::B::clade(entry.sequence, fasta_ref);
            }
            else if (subtype == "H1") {
                local::H1::deletions(entry.sequence, fasta_ref);
                local::H1::clade(entry.sequence, fasta_ref);
            }
            else if (subtype == "H3") {
                local::H3::deletions(entry.sequence, fasta_ref);
                local::H3::clade(entry.sequence, fasta_ref);
            }
        }
    }
        // detect B lineage and VIC deletion mutants, adjust deletions
        // detect clades

} // acmacs::seqdb::v3::scan::detect_lineages_clades

// ****************************************************************************************************
// B
// ****************************************************************************************************

namespace local::B
{
    using deletions_insertions_t = acmacs::seqdb::v3::scan::deletions_insertions_t;

    inline bool no_deletions_after_before(const deletions_insertions_t& deletions, acmacs::seqdb::pos1_t on_or_after_pos, acmacs::seqdb::pos1_t on_or_before_pos)
    {
        for (const auto& del : deletions.deletions) {
            if (del.pos >= on_or_after_pos && del.pos <= on_or_before_pos)
                return false;
        }
        return true;
    }

    inline bool N_deletions_at(const deletions_insertions_t& deletions, size_t num_deletions, acmacs::seqdb::pos1_t pos)
    {
        return deletions.deletions.front().pos == pos && deletions.deletions.front().num == num_deletions && deletions.insertions.empty();
    }

    // ----------------------------------------------------------------------

    inline bool is_yamagata_shifted(acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return ((deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{158} && deletions.front().num == 1 && sequence.aa_aligned_substr(155, 6) == "MAWVIP") ||
                (deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{161} && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 2) == "VP") ||
                (deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{160} && deletions.front().num == 1 && sequence.aa_aligned_substr(157, 3) == "WAV") ||
                (deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{163} && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 3) == "VPK")) &&
               sequence.deletions().insertions.empty();
    }

    inline bool is_semi_ignored(acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        return *sequence.name() == "B/MIE/1/2019" // DEL[1](162:4)<pos-1-based>  NIID:20190314   -- B/MIE/1/2019 |  2019-01-22 | MDCK 1 +1 |  18/19-498 | National Institute of Infectious Diseases
                                                  // (NIID) | B / H0N0 |  Victoria
               || *sequence.name() == "B/INDONESIA/NIHRDSB183950/2018" // DEL[1](164:2)<pos-1-based> VIDRL:20180913 -- B/Indonesia/NIHRDSB183950/2018 |  2018-04-01 | X, MDCK1 |  10004643 VW10005052 |
                                                                       // WHO Collaborating Centre for Reference and Research on Influenza | B / H0N0 |  Victoria
            ;
    }

    inline bool is_ignored(acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        return *sequence.name() ==
                   "B/ONTARIO/RV1769/2019" // DEL[1](163:3)<pos-1-based>  B/Ontario/RV1769/2019 |  2019-04-11 | P1 |  RV1769/19 | Public Health Agency of Canada (PHAC) | B / H0N0 |  Victoria
               || *sequence.name() == "B/KENYA/4/2018"  // DEL[1](160:1)<pos-1-based>  B/Kenya/004/2018 |  2018-01-05 |  |   | Other Database Import | B / H0N0 |  unknown
               || *sequence.name() == "B/KENYA/11/2018" // DEL[1](160:1)<pos-1-based>  B/Kenya/011/2018 |  2018-01-15 |  |   | Other Database Import | B / H0N0 |  unknown
               || *sequence.name() ==
                      "B/ORENBURG/CRIE-100/2018" // DEL[1](160:1)<pos-1-based>  B/Orenburg/CRIE/100/2018 |  2018-02-08 |  |   | Central Research Institute of Epidemiology | B / H0N0 |  Yamagata
            ;
    }

    // B/Yamagata/16/88
    // B/Victoria/2/87

    // YAMAGATA: deletion must be at 163
    // David Burke 2017-08-17: deletions ( and insertions) of amino acids usually occur in regions of the protein structure where it changes direction ( loops ).
    // In the case of HA, this is after VPK and before NKTAT/YKNAT.

    // DISABLED:
    //     Sarah 2018-08, David Burke disagrees 2019-07-16
    //     VICTORIA del2017: 162, 163
    //     VICTORIA tripledel2017: 162, 163, 164 by convention

    void lineage(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref)
    {
        const auto warn = [&](const char* infix, const char* prefix = "WARNING") {
            fmt::print(stderr, "{}: {} lineage {} and {} deletions {} {}\n{}\n{}\n", prefix, sequence.year(), *sequence.lineage(), infix, sequence.full_name(),
                       acmacs::seqdb::scan::format(sequence.deletions()), fasta_ref, sequence.aa_format());
        };

        // const auto rep = [&]() {
        //     fmt::print(stderr, "{}\n{}\n{}\n", sequence.full_name(), sequence.aa_format(), sequence.nuc_format());
        // };

        auto& deletions = sequence.deletions();
        constexpr acmacs::seqdb::pos1_t b_vic_del_mutants_pos{162}; // Must be 162 according to Sarah and CDC

        //---------- VICTORIA ----------

        if (no_deletions_after_before(deletions, acmacs::seqdb::pos1_t{1}, acmacs::seqdb::pos1_t{500})) { // may have deletions after 500
            // VICTORIA
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
                warn("no");
        }
        else if (N_deletions_at(deletions, 2, acmacs::seqdb::pos1_t{162}) || N_deletions_at(deletions, 2, acmacs::seqdb::pos1_t{163})) {
            // VICTORIA (double) del 2017
            // according to David Burke 2019-07-16 14:27, also see https://jvi.asm.org/content/jvi/73/9/7343.full.pdf
            // B/GUATEMALA/581/2017      VPN--KNKTAT
            // B/COLORADO/6/2017_MDCK1   VPD--KNKTAT
            deletions.deletions.front().pos = b_vic_del_mutants_pos;
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
                warn("victoria del2017");
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"DEL2017"});
        }
        else if (N_deletions_at(deletions, 3, acmacs::seqdb::pos1_t{162}) || N_deletions_at(deletions, 3, acmacs::seqdb::pos1_t{163}) || N_deletions_at(deletions, 3, acmacs::seqdb::pos1_t{164})) {
            // VICTORIA triple del 2017
            // according to David Burke 2019-07-16 14:27
            // VPK---NKTAT
            deletions.deletions.front().pos = b_vic_del_mutants_pos;
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
                warn("victoria tripledel2017");
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"TRIPLEDEL2017"});
        }
        else if (N_deletions_at(deletions, 6, acmacs::seqdb::pos1_t{164})) {
            // VICTORIA sixdel2019 (only from Japan as of 2019-07-19)
            // David Burke 2019-07-19 15:40: These look really
            // unusual. Based on the geometry of the loop, I would
            // tend to align the N with C-terminal side: B/KANAGAWA/AC1867/2019 VPK------NTNP
            deletions.deletions.front().pos = b_vic_del_mutants_pos;
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
                warn("victoria sixdel2019 (pos shifted)");
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"SIXDEL2019"});
            // rep();
        }

            //---------- YAMAGATA ----------

        else if (N_deletions_at(deletions, 1, acmacs::seqdb::pos1_t{163}) && no_deletions_after_before(deletions, acmacs::seqdb::pos1_t{164}, acmacs::seqdb::pos1_t{500})) {
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
                warn("yamagata");
        }
        else if (is_yamagata_shifted(sequence)) {
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
                warn("yamagata-shifted");
            deletions.deletions = std::vector<acmacs::seqdb::scan::deletions_insertions_t::pos_num_t>{{acmacs::seqdb::pos1_t{162}, 1}};
        }
        else if (N_deletions_at(deletions, 2, acmacs::seqdb::pos1_t{163}) && sequence.year() <= 2013) {
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
                warn("yamagata");
        }
        else if (N_deletions_at(deletions, 2, acmacs::seqdb::pos1_t{169})) {
            // 12 sequences from TAIWAN 2010 have deletions 169:2
            sequence.lineage(acmacs::virus::lineage_t{});
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"TAIWAN2010"});
        }
        else if (N_deletions_at(deletions, 1, acmacs::seqdb::pos1_t{160}) && no_deletions_after_before(deletions, acmacs::seqdb::pos1_t{161}, acmacs::seqdb::pos1_t{500}) && sequence.aa_at_pos(acmacs::seqdb::pos1_t{161}) == 'E' && sequence.aa_at_pos(acmacs::seqdb::pos1_t{163}) == 'K') {
            // deletion detection was invalid, most probably due to 162X. B/ALICANTE/19_0649/20171219
            // fmt::print(stderr, "DEBUG: {} 160{} 161{} 162{} 163{}\n", acmacs::seqdb::scan::format(deletions), sequence.aa_at_pos(acmacs::seqdb::pos1_t{160}), sequence.aa_at_pos(acmacs::seqdb::pos1_t{161}), sequence.aa_at_pos(acmacs::seqdb::pos1_t{162}), sequence.aa_at_pos(acmacs::seqdb::pos1_t{163}));
            sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
            deletions.deletions = std::vector<acmacs::seqdb::scan::deletions_insertions_t::pos_num_t>{{acmacs::seqdb::pos1_t{162}, 1}};
        }
        else if (is_semi_ignored(sequence)) {
            fmt::print(stderr, "INFO: {} {}\n", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
        }
        else if (is_ignored(sequence)) {
            // do not issue warning
        }
        else {
            fmt::print(stderr, "DEBUG: 1-at-163:{} no-between-164-500:{}\n", N_deletions_at(deletions, 1, acmacs::seqdb::pos1_t{163}), no_deletions_after_before(deletions, acmacs::seqdb::pos1_t{164}, acmacs::seqdb::pos1_t{500}));
            warn("unknown", "ERROR");
        }

    } // lineage

    // ----------------------------------------------------------------------

    // V1A: !58P 75K 172P
    // V1A.1 = del2017
    // V1A.2 = tripledel2017
    // V1B  58P

    void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view /*fasta_ref*/)
    {
        if (sequence.lineage() == acmacs::virus::lineage_t{"VICTORIA"}) {
            // 2018-09-03, Sarah: clades should (technically) be defined by a phylogenetic tree rather than a set of amino acids
            if (sequence.aa_at_pos(acmacs::seqdb::pos1_t{75}) == 'K' && sequence.aa_at_pos(acmacs::seqdb::pos1_t{172}) == 'P' && sequence.aa_at_pos(acmacs::seqdb::pos1_t{58}) != 'P')
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"V1A"});
            else if (sequence.aa_at_pos(acmacs::seqdb::pos1_t{58}) == 'P')
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"V1B"});
            else
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"V1"});
        }
        else if (sequence.lineage() == acmacs::virus::lineage_t{"YAMAGATA"}) {
            // 165N -> Y2, 165Y -> Y3 (yamagata numeration, 163 is not -)
            // 166N -> Y2, 166Y -> Y3 (victoria numeration, 163 is -)
            switch (sequence.aa_at_pos(acmacs::seqdb::pos1_t{166})) {
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

    } // clade

} // namespace local::B

// ****************************************************************************************************
// H1
// ****************************************************************************************************

namespace local::H1
{
    // inline bool is_127_1(const acmacs::seqdb::v3::deletions_insertions_t& deletions)
    // {
    //     return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 126 && deletions.deletions.front().num == 1 && deletions.insertions.empty();
    // }

    void deletions(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref)
    {
        const auto warn = [&](const char* prefix = "WARNING") {
            fmt::print(stderr, "{}: {} {} {} {} :: {}\n{}\n", prefix, sequence.year(), sequence.date_simulated(), sequence.full_name(), acmacs::seqdb::scan::format(sequence.deletions()), fasta_ref,
                       sequence.aa_format());
        };

        const auto& deletions = sequence.deletions();
        const auto host = acmacs::virus::host(sequence.name());
        const auto year = sequence.year();
        if (deletions.deletions.size() == 1) {
            const auto& del1 = deletions.deletions[0];
            if (sequence.type_subtype() == acmacs::virus::type_subtype_t{"A(H1N2)"} || !host.empty() || year < 2010)
                sequence.add_clade(acmacs::seqdb::clade_t{"*DEL"});
            else if (del1.pos == acmacs::seqdb::pos1_t{127} && del1.num == 1 && (year < 2018 || fasta_ref.find("seasonal") != std::string_view::npos))
                sequence.add_clade(acmacs::seqdb::clade_t{"*DEL-127:1"});
            else if (del1.pos == acmacs::seqdb::pos1_t{160} && del1.num == 4 && sequence.name() == acmacs::virus::virus_name_t{"A(H1N1)/NEWPORT/323/2019"})
                fmt::print(stderr, "INFO: {} {}\n", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
            else if (del1.pos > acmacs::seqdb::pos1_t{400})
                ; // ignore
            else
                warn();
        }
        else if (deletions.deletions.size() > 1) {
            if (!host.empty() || year < 2010)
                sequence.add_clade(acmacs::seqdb::clade_t{"*DEL"});
            else
                warn();
        }
        else if (!deletions.insertions.empty())
            sequence.add_clade(acmacs::seqdb::clade_t{"*INS"});
        else if (!deletions.empty()) {
            warn();
        }

    } // deletions

    // ----------------------------------------------------------------------
    // Before 2018-09-19
    // ----------------------------------------------------------------------
    //   // 84N+162N+216T - 6B.1, 152T+173I+501E - 6B.2
    //   // ? 156 (see A/PUERTO RICO/15/2018 of CDC:20180511)

    // ----------------------------------------------------------------------
    // 2018-09-19 clade definitions changed by Sarah before SSM
    // ----------------------------------------------------------------------
    // 6B: 163Q
    // 6B1: 162N, 163Q
    // 6B2: 152T, 163Q

    void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view /*fasta_ref*/)
    {
        if (sequence.aa_at_pos(acmacs::seqdb::pos1_t{163}) == 'Q') {
            sequence.add_clade(acmacs::seqdb::clade_t{"6B"});
            if (sequence.aa_at_pos(acmacs::seqdb::pos1_t{162}) == 'N')
                sequence.add_clade(acmacs::seqdb::clade_t{"6B1"});
            if (sequence.aa_at_pos(acmacs::seqdb::pos1_t{74}) == 'R' && sequence.aa_at_pos(acmacs::seqdb::pos1_t{164}) == 'T' && sequence.aa_at_pos(acmacs::seqdb::pos1_t{295}) == 'V')
                sequence.add_clade(acmacs::seqdb::clade_t{"6B1.A"});
            if (sequence.aa_at_pos(acmacs::seqdb::pos1_t{152}) == 'T')
                sequence.add_clade(acmacs::seqdb::clade_t{"6B2"});
        }

    } // clade

} // namespace local::H1

// ****************************************************************************************************
// H3
// ****************************************************************************************************

namespace local::H3
{

    void deletions(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref)
    {
        const auto warn = [&](const char* prefix = "WARNING") {
            fmt::print(stderr, "{}: {} <{}> {} {} :: {}\n{}\n", prefix, sequence.year(), sequence.aa_aligned_length(), sequence.full_name(), acmacs::seqdb::scan::format(sequence.deletions()), fasta_ref,
                       sequence.aa_format());
        };

        const auto& deletions = sequence.deletions();
        if (!deletions.insertions.empty())
            sequence.add_clade(acmacs::seqdb::clade_t{"*INS"});
        else if (!deletions.empty()) {
            if (sequence.aa_aligned_length() < 500)
                ; // ignore short
            else if (!acmacs::virus::host(sequence.name()).empty())
                ; // ignore
            else if (sequence.year() < 2018)
                sequence.add_clade(acmacs::seqdb::clade_t{"*DEL"});
            else
                warn();
        }

    } // deletions

    // ----------------------------------------------------------------------

    using p1 = acmacs::seqdb::pos1_t;

    struct clade_desc_t
    {
        struct pos1_aa_t
        {
            p1 pos1;
            char aa;
        };

        acmacs::seqdb::clade_t clade;
        std::vector<pos1_aa_t> pos1_aa;
    };

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

    static const std::array sClades{
        clade_desc_t{acmacs::seqdb::clade_t{"3C.3"}, {{p1{158}, 'N'}, {p1{159}, 'F'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"3A"},   {{p1{138}, 'S'}, {p1{159}, 'S'}, {p1{225}, 'D'}}}, // R326K causes split in the tree for 2019-0814-tc1, removed on 2019-08-21 {p1{326}, 'R'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"3B"},   {{ p1{62}, 'K'}, { p1{83}, 'R'}, {p1{261}, 'Q'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A"},   {{p1{158}, 'N'}, {p1{159}, 'Y'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1"},  {{p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{406}, 'V'}, {p1{484}, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1A"}, {{p1{121}, 'K'}, {p1{135}, 'K'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{406}, 'V'}, {p1{479}, 'E'}, {p1{484}, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1B"}, {{ p1{92}, 'R'}, {p1{121}, 'K'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{311}, 'Q'}, {p1{406}, 'V'}, {p1{484}, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1B-135K"},              {{ p1{92}, 'R'}, {p1{121}, 'K'}, {p1{135}, 'K'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{311}, 'Q'}, {p1{406}, 'V'}, {p1{484}, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1B-135K-137F-193S"},    {{ p1{92}, 'R'}, {p1{121}, 'K'}, {p1{135}, 'K'}, {p1{137}, 'F'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{193}, 'S'}, {p1{311}, 'Q'}, {p1{406}, 'V'}, {p1{484}, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1B-131K"},              {{ p1{92}, 'R'}, {p1{121}, 'K'}, {p1{131}, 'K'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{311}, 'Q'}, {p1{406}, 'V'}, {p1{484}, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A2"},  {{p1{131}, 'K'}, {p1{142}, 'K'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{261}, 'Q'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A3"},  {{p1{121}, 'K'}, {p1{135}, 'K'}, {p1{144}, 'K'}, {p1{150}, 'K'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{261}, 'Q'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A4"},  {{ p1{31}, 'S'}, { p1{53}, 'N'}, {p1{142}, 'G'}, {p1{144}, 'R'}, {p1{158}, 'N'}, {p1{159}, 'Y'}, {p1{171}, 'K'}, {p1{192}, 'T'}, {p1{197}, 'H'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"159S"}, {{p1{159}, 'S'}}}, // explicit Derek's request on 2019-04-18
        clade_desc_t{acmacs::seqdb::clade_t{"159F"}, {{p1{159}, 'F'}}}, // explicit Derek's request on 2019-04-18
        clade_desc_t{acmacs::seqdb::clade_t{"159Y"}, {{p1{159}, 'Y'}}}, // explicit Derek's request on 2019-04-18
    };

    // Removed because it makes no sense, GLY cannot be difined this way, search email for sequon
        // clade_desc_t{acmacs::seqdb::clade_t{"GLY"},  {{p1{160}, 'S'}}},
        // clade_desc_t{acmacs::seqdb::clade_t{"GLY"},  {{p1{160}, 'T'}}},

#pragma GCC diagnostic pop

    void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view /*fasta_ref*/)
    {
        const auto has_aa = [&](const auto& pos1_aa) -> bool { return sequence.aa_at_pos(pos1_aa.pos1) == pos1_aa.aa; };

        for (const auto& clade_desc : sClades) {
            if (std::all_of(clade_desc.pos1_aa.begin(), clade_desc.pos1_aa.end(), has_aa))
                sequence.add_clade(clade_desc.clade);
        }

    } // clade

} // namespace local::H3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
