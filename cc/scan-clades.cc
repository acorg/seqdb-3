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
            const auto fasta_ref = fmt::format("{}:{}  {}", entry.fasta.filename, entry.fasta.line_no, entry.fasta.entry_name);
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
    inline bool is_victoria(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions) { return deletions.empty(); }

    inline bool is_victoria_del2017(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 161 && deletions.deletions.front().num == 2 && deletions.insertions.empty();
    }

    inline bool is_victoria_tripledel2017(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 161 && deletions.deletions.front().num == 3 && deletions.insertions.empty();
    }

    inline bool is_victoria_tripledel2017_pos_shifted_163_164(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && (deletions.deletions.front().pos == 162 || deletions.deletions.front().pos == 163) && deletions.deletions.front().num == 3 && deletions.insertions.empty();
    }

    inline bool is_victoria_sixdel2019(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 163 && deletions.deletions.front().num == 6 && deletions.insertions.empty();
    }

    inline bool is_victoria_deletions_at_the_end(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos > 500 && deletions.insertions.empty();
    }

    inline bool is_yamagata_shifted(acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return ((deletions.size() == 1 && deletions.front().pos == 158 && deletions.front().num == 1 && sequence.aa_aligned_substr(155, 6) == "MAWVIP") ||
                (deletions.size() == 1 && deletions.front().pos == 161 && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 2) == "VP") ||
                (deletions.size() == 1 && deletions.front().pos == 160 && deletions.front().num == 1 && sequence.aa_aligned_substr(157, 3) == "WAV") ||
                (deletions.size() == 1 && deletions.front().pos == 163 && deletions.front().num == 1 && sequence.aa_aligned_substr(159, 3) == "VPK")) &&
               sequence.deletions().insertions.empty();
    }

    inline bool is_yamagata(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return !deletions.deletions.empty() && deletions.deletions.front().pos == 162 && deletions.deletions.front().num == 1 &&
               (deletions.deletions.size() == 1 || deletions.deletions[1].pos > 500) && deletions.insertions.empty();
    }

    inline bool is_yamagata_doubledel(acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return deletions.size() == 1 && deletions.front().pos == 162 && deletions.front().num == 2 && sequence.year() <= 2013 && sequence.deletions().insertions.empty();
    }

    // 12 sequences from TAIWAN 2010 have deletions 169:2
    inline bool is_taiwan_169_2(const acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return deletions.deletions.size() == 1 && deletions.deletions.front().pos == 168 && deletions.deletions.front().num == 2 && deletions.insertions.empty();
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

    // VICTORIA del2017: 162, 163
    // VICTORIA tripledel2017: 162, 163, 164 by convention

    // YAMAGATA: deletion must be at 163
    // David Burke 2017-08-17: deletions ( and insertions) of amino acids usually occur in regions of the protein structure where it changes direction ( loops ).
    // In the case of HA, this is after VPK and before NKTAT/YKNAT.

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
        else if (is_victoria_tripledel2017_pos_shifted_163_164(deletions)) {
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
                warn("victoria tripledel2017 (pos shifted)");
            deletions.deletions.front().pos = 161;
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"TRIPLEDEL2017"});
            // rep();
        }
        else if (is_victoria_sixdel2019(deletions)) {
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"VICTORIA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"VICTORIA"})
                warn("victoria sixdel2019 (pos shifted)");
            sequence.add_clade(acmacs::seqdb::v3::clade_t{"SIXDEL2019"});
            // rep();
        }
        else if (is_yamagata_shifted(sequence)) {
            if (sequence.lineage().empty())
                sequence.lineage(acmacs::virus::lineage_t{"YAMAGATA"});
            else if (sequence.lineage() != acmacs::virus::lineage_t{"YAMAGATA"})
                warn("yamagata-shifted");
            deletions.deletions = std::vector<acmacs::seqdb::scan::deletions_insertions_t::pos_num_t>{{162, 1}};
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
            fmt::print(stderr, "INFO: {} {}\n", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
        }
        else if (is_ignored(sequence)) {
            // do not issue warning
        }
        else {
            warn("unknown", "ERROR");
        }

    } // lineage

    // ----------------------------------------------------------------------

    void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view /*fasta_ref*/)
    {
        if (sequence.lineage() == acmacs::virus::lineage_t{"VICTORIA"}) {
            // 2018-09-03, Sarah: clades should (technically) be defined by a phylogenetic tree rather than a set of amino acids
            if (sequence.aa_at_pos1(75) == 'K' && sequence.aa_at_pos1(172) == 'P' && sequence.aa_at_pos1(58) != 'P')
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"V1A"});
            else if (sequence.aa_at_pos1(58) == 'P')
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"V1B"});
            else
                sequence.add_clade(acmacs::seqdb::v3::clade_t{"V1"});
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
            else if (del1.pos == 126 && del1.num == 1 && (year < 2018 || fasta_ref.find("seasonal") != std::string_view::npos))
                sequence.add_clade(acmacs::seqdb::clade_t{"*DEL-127:1"});
            else if (del1.pos == 159 && del1.num == 4 && sequence.name() == acmacs::virus::virus_name_t{"A(H1N1)/NEWPORT/323/2019"})
                fmt::print(stderr, "INFO: {} {}\n", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
            else if (del1.pos > 400)
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
        if (sequence.aa_at_pos1(163) == 'Q') {
            sequence.add_clade(acmacs::seqdb::clade_t{"6B"});
            if (sequence.aa_at_pos1(162) == 'N')
                sequence.add_clade(acmacs::seqdb::clade_t{"6B1"});
            if (sequence.aa_at_pos1(74) == 'R' && sequence.aa_at_pos1(164) == 'T' && sequence.aa_at_pos1(295) == 'V')
                sequence.add_clade(acmacs::seqdb::clade_t{"6B1.A"});
            if (sequence.aa_at_pos1(152) == 'T')
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

    struct clade_desc_t
    {
        struct pos1_aa_t
        {
            size_t pos1;
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
        clade_desc_t{acmacs::seqdb::clade_t{"3C.3"}, {{158, 'N'}, {159, 'F'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"3A"},   {{138, 'S'}, {159, 'S'}, {225, 'D'}, {326, 'R'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"3B"},   {{ 62, 'K'}, { 83, 'R'}, {261, 'Q'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A"},   {{158, 'N'}, {159, 'Y'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1"},  {{158, 'N'}, {159, 'Y'}, {171, 'K'}, {406, 'V'}, {484, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1A"}, {{121, 'K'}, {135, 'K'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {406, 'V'}, {479, 'E'}, {484, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A1B"}, {{ 92, 'R'}, {121, 'K'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {311, 'Q'}, {406, 'V'}, {484, 'E'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A2"},  {{131, 'K'}, {142, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A3"},  {{121, 'K'}, {135, 'K'}, {144, 'K'}, {150, 'K'}, {158, 'N'}, {159, 'Y'}, {261, 'Q'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"2A4"},  {{ 31, 'S'}, { 53, 'N'}, {142, 'G'}, {144, 'R'}, {158, 'N'}, {159, 'Y'}, {171, 'K'}, {192, 'T'}, {197, 'H'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"GLY"},  {{160, 'S'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"GLY"},  {{160, 'T'}}},
        clade_desc_t{acmacs::seqdb::clade_t{"159S"}, {{159, 'S'}}}, // explicit Derek's request on 2019-04-18
        clade_desc_t{acmacs::seqdb::clade_t{"159F"}, {{159, 'F'}}}, // explicit Derek's request on 2019-04-18
        clade_desc_t{acmacs::seqdb::clade_t{"159Y"}, {{159, 'Y'}}}, // explicit Derek's request on 2019-04-18
    };

#pragma GCC diagnostic pop

    void clade(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view /*fasta_ref*/)
    {
        const auto has_aa = [&](const auto& pos1_aa) -> bool { return sequence.aa_at_pos1(pos1_aa.pos1) == pos1_aa.aa; };

        for (const auto& clade_desc : sClades) {
            if (std::all_of(clade_desc.pos1_aa.begin(), clade_desc.pos1_aa.end(), has_aa))
                sequence.add_clade(clade_desc.clade);
        }

    } // clade

} // namespace local::H3

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
