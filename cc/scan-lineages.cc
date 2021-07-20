#include <map>

#include "acmacs-base/filesystem.hh"
#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/settings-v3.hh"
#include "acmacs-base/rjson-v3.hh"
#include "acmacs-base/counter.hh"
#include "seqdb-3/scan-lineages.hh"
#include "seqdb-3/scan-fasta.hh"
#include "seqdb-3/aa-at-pos.hh"

// ----------------------------------------------------------------------

class CladeDefinitions : public acmacs::settings::v3::Data
{
  public:
    CladeDefinitions()
    {
        using namespace std::string_literals;
        using namespace std::string_view_literals;

        if (const auto filename = fmt::format("{}/share/conf/clades.json", acmacs::acmacsd_root()); fs::exists(filename))
            acmacs::settings::v3::Data::load(filename);
        else
            throw std::runtime_error{fmt::format("WARNING: cannot load \"{}\": file not found\n", filename)};

        using pp = std::pair<std::string, std::string_view>;
        for (const auto& [virus_type, tag] : {pp{"H1"s, "clades-A(H1N1)2009pdm"sv}, pp{"H3"s, "clades-A(H3N2)"sv}, pp{"BVICTORIA"s, "clades-B/Vic"sv}, pp{"BYAMAGATA"s, "clades-B/Yama"sv}}) {
            current_virus_type_ = virus_type;
            apply(tag);
        }
    }

    bool apply_built_in(std::string_view name) override // returns true if built-in command with that name found and applied
    {
        using namespace std::string_view_literals;
        if (name == "clade"sv) {
            const auto& aa_field = getenv("aa"sv);
            auto aa = aa_field.is_null() ? acmacs::seqdb::amino_acid_at_pos1_eq_list_t{} : acmacs::seqdb::extract_aa_at_pos1_eq_list(aa_field);
            const auto& nuc_field = getenv("nuc"sv);
            auto nuc = nuc_field.is_null() ? acmacs::seqdb::nucleotide_at_pos1_eq_list_t{} : acmacs::seqdb::extract_nuc_at_pos1_eq_list(nuc_field);
            add(current_virus_type_, getenv_or("name"sv, ""sv), std::move(aa), std::move(nuc));
        }
        else
            return acmacs::settings::v3::Data::apply_built_in(name);
        return true;
    }

    bool matches(const acmacs::seqdb::v3::scan::sequence_t& seq, const acmacs::seqdb::amino_acid_at_pos1_eq_list_t& aa_at_pos) const
    {
        return std::all_of(std::begin(aa_at_pos), std::end(aa_at_pos), [&seq](const acmacs::seqdb::amino_acid_at_pos1_eq_t& pos1_aa) {
            return (seq.aa_at_pos(std::get<acmacs::seqdb::pos1_t>(pos1_aa)) == std::get<char>(pos1_aa)) == std::get<bool>(pos1_aa);
        });
    }

    bool matches(const acmacs::seqdb::v3::scan::sequence_t& seq, const acmacs::seqdb::nucleotide_at_pos1_eq_list_t& nuc_at_pos) const
    {
        return std::all_of(std::begin(nuc_at_pos), std::end(nuc_at_pos), [&seq](const acmacs::seqdb::nucleotide_at_pos1_eq_t& pos1_nuc) {
            return (seq.nuc_at_pos(std::get<acmacs::seqdb::pos1_t>(pos1_nuc)) == std::get<char>(pos1_nuc)) == std::get<bool>(pos1_nuc);
        });
    }

    void add_clades(acmacs::seqdb::v3::scan::sequence_t& sequence, const std::string& virus_type)
    {
        if (auto found = data_.find(virus_type); found != std::end(data_)) {
            for (const auto& [clade_name, aa_at_pos, nuc_at_pos] : found->second) {
                if ((aa_at_pos.empty() || matches(sequence, aa_at_pos)) && (nuc_at_pos.empty() || matches(sequence, nuc_at_pos)))
                    sequence.add_clade(clade_name);
            }
        }
        else
            AD_WARNING("no clade definitions for {} seq: {}", virus_type, sequence.name());
    }

  private:
    std::string current_virus_type_;
    std::map<std::string, std::vector<std::tuple<acmacs::seqdb::v3::clade_t, acmacs::seqdb::amino_acid_at_pos1_eq_list_t, acmacs::seqdb::nucleotide_at_pos1_eq_list_t>>, std::less<>> data_;

    void add(const std::string& virus_type, std::string_view clade_name, acmacs::seqdb::amino_acid_at_pos1_eq_list_t&& aa_at_pos, acmacs::seqdb::nucleotide_at_pos1_eq_list_t&& nuc_at_pos)
    {
        data_[virus_type].emplace_back(acmacs::seqdb::v3::clade_t{clade_name}, std::move(aa_at_pos), std::move(nuc_at_pos));
    }
};

// ----------------------------------------------------------------------

namespace local
{
    namespace B
    {
        static void lineage(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref, const acmacs::virus::lineage_t& fasta_lineage);
    } // namespace B

    namespace H1
    {
        static void deletions(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
    } // namespace H1

    namespace H3
    {
        static void deletions(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref);
    } // namespace H3

} // namespace local

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::detect_lineages_clades(std::vector<fasta::scan_result_t>& sequences)
{

    CladeDefinitions clade_definitions;

#pragma omp parallel for default(shared) schedule(static, 256)
    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        if (auto& entry = sequences[e_no]; !entry.reference && fasta::is_aligned(entry)) {
            const auto subtype = entry.sequence.type_subtype().h_or_b();
            const auto fasta_ref = fmt::format("{}:{}: note:  {}", entry.fasta.filename, entry.fasta.line_no, entry.fasta.entry_name);
            if (subtype == "B") {
                local::B::lineage(entry.sequence, fasta_ref, entry.fasta.lineage);
                if (!entry.sequence.lineage().empty()) // no clade definitions without lineage
                    clade_definitions.add_clades(entry.sequence, fmt::format("{}{}", subtype, entry.sequence.lineage()));
            }
            else if (subtype == "H1") {
                local::H1::deletions(entry.sequence, fasta_ref);
                clade_definitions.add_clades(entry.sequence, std::string{subtype});
            }
            else if (subtype == "H3") {
                local::H3::deletions(entry.sequence, fasta_ref);
                clade_definitions.add_clades(entry.sequence, std::string{subtype});
            }
        }
    }

    // populate lineage for references
    std::map<acmacs::virus::name_t, std::vector<fasta::scan_result_t*>> referenced;
    for (auto& seq : sequences) {
        if (seq.reference)
            referenced[seq.reference->name].push_back(&seq);
    }
    for (auto& seq : sequences) {
        if (const auto found = referenced.find(seq.sequence.name()); found != referenced.end() && !seq.sequence.lineage().empty()) {
            for (auto* ref : found->second)
                ref->sequence.lineage(seq.sequence.lineage());
        }
    }

} // acmacs::seqdb::v3::scan::detect_lineages_clades

// ****************************************************************************************************
// B
// ****************************************************************************************************

namespace local::B
{
    using pos0_t = acmacs::seqdb::pos0_t;
    using pos1_t = acmacs::seqdb::pos1_t;
    using deletions_insertions_t = acmacs::seqdb::v3::scan::deletions_insertions_t;
    using pos_num_t = acmacs::seqdb::v3::scan::deletions_insertions_t::pos_num_t;
    using v_pos_num_t = std::vector<acmacs::seqdb::v3::scan::deletions_insertions_t::pos_num_t>;

    constexpr pos1_t b_vic_del_mutants_pos{162}; // Must be 162 according to Sarah and CDC

    inline bool no_deletions_after_before(const deletions_insertions_t& deletions, pos1_t on_or_after_pos, pos1_t on_or_before_pos)
    {
        for (const auto& del : deletions.deletions) {
            if (del.pos >= on_or_after_pos && del.pos <= on_or_before_pos)
                return false;
        }
        return true;
    }

    inline bool N_deletions_at(const deletions_insertions_t& deletions, size_t num_deletions, pos1_t pos)
    {
        return !deletions.deletions.empty() && deletions.deletions.front().pos == pos && deletions.deletions.front().num == num_deletions && deletions.insertions.empty();
    }

    inline bool N_deletions_at(const deletions_insertions_t& deletions, size_t num_deletions, pos1_t pos_min, pos1_t pos_max)
    {
        return !deletions.deletions.empty() && deletions.deletions.front().pos >= pos_min && deletions.deletions.front().pos <= pos_max && deletions.deletions.front().num == num_deletions && deletions.insertions.empty();
    }

    inline void replace_front_deletions(deletions_insertions_t& deletions, const v_pos_num_t& replace_with)
    {
        deletions.deletions.erase(deletions.deletions.begin());
        deletions.deletions.insert(deletions.deletions.begin(), std::begin(replace_with), std::end(replace_with));
    }

    // ----------------------------------------------------------------------

    inline bool is_yamagata_shifted(const acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        const auto& deletions = sequence.deletions().deletions;
        return ((deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{158} && deletions.front().num == 1 && sequence.aa_aligned_without_deletions_substr(155, 6) == "MAWVIP") ||
                (deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{161} && deletions.front().num == 1 && sequence.aa_aligned_without_deletions_substr(159, 2) == "VP") ||
                (deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{160} && deletions.front().num == 1 && sequence.aa_aligned_without_deletions_substr(157, 3) == "WAV") ||
                (deletions.size() == 1 && deletions.front().pos == acmacs::seqdb::pos0_t{163} && deletions.front().num == 1 && sequence.aa_aligned_without_deletions_substr(159, 3) == "VPK")) &&
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

    inline bool is_ignored(const acmacs::seqdb::v3::scan::sequence_t& sequence)
    {
        return *sequence.name() ==
                   "B/ONTARIO/RV1769/2019" // DEL[1](163:3)<pos-1-based>  B/Ontario/RV1769/2019 |  2019-04-11 | P1 |  RV1769/19 | Public Health Agency of Canada (PHAC) | B / H0N0 |  Victoria
               || *sequence.name() == "B/KENYA/4/2018"  // DEL[1](160:1)<pos-1-based>  B/Kenya/004/2018 |  2018-01-05 |  |   | Other Database Import | B / H0N0 |  unknown
               || *sequence.name() == "B/KENYA/11/2018" // DEL[1](160:1)<pos-1-based>  B/Kenya/011/2018 |  2018-01-15 |  |   | Other Database Import | B / H0N0 |  unknown
               || *sequence.name() ==
                      "B/ORENBURG/CRIE-100/2018" // DEL[1](160:1)<pos-1-based>  B/Orenburg/CRIE/100/2018 |  2018-02-08 |  |   | Central Research Institute of Epidemiology | B / H0N0 |  Yamagata
            ;
    }

    // ----------------------------------------------------------------------

    inline void set_lineage(acmacs::seqdb::v3::scan::sequence_t& sequence, const acmacs::virus::lineage_t& lineage, std::string_view fasta_ref, const char* warn_infix)
    {
        if (sequence.lineage().empty())
            sequence.lineage(lineage);
        else if (sequence.lineage() != lineage)
            AD_WARNING("{} lineage {} and {} deletions {} {}\n{}\n{}", sequence.year(), *sequence.lineage(), warn_infix, sequence.full_name(), acmacs::seqdb::scan::format(sequence.deletions()),
                       fasta_ref, sequence.aa_format());
    }

    // ----------------------------------------------------------------------

    // B/Yamagata/16/88
    // B/Victoria/2/87

    // YAMAGATA: deletion must be at 163
    // David Burke 2017-08-17: deletions ( and insertions) of amino acids usually occur in regions of the protein structure where it changes direction ( loops ).
    // In the case of HA, this is after VPK and before NKTAT/YKNAT.

    // DISABLED:
    //     Sarah 2018-08, David Burke disagrees 2019-07-16
    //     VICTORIA del2017: 162, 163
    //     VICTORIA tripledel2017: 162, 163, 164 by convention

    void lineage(acmacs::seqdb::v3::scan::sequence_t& sequence, std::string_view fasta_ref, const acmacs::virus::lineage_t& fasta_lineage)
    {
        auto& deletions = sequence.deletions();
        const auto date{sequence.date_simulated()};

        // AD_DEBUG("{} {}", sequence.full_name(), acmacs::seqdb::scan::format(deletions));

        // if (sequence.name() == acmacs::virus::name_t{"B/BEIJING/15/1984"})
        //     AD_DEBUG("{} {} {}", sequence.name(), sequence.aa_aligned_substr(159, 14), deletions);

        //---------- Special deletions ----------

        bool special_deletions = false;
        if (N_deletions_at(deletions, 2, pos1_t{162}, pos1_t{164})) {
            const auto at162{sequence.aa_at_pos_without_deletions(pos1_t{162})}, at163{sequence.aa_at_pos_without_deletions(pos1_t{163})};
            if ((at162 == 'R' || at162 == 'K') && at163 == 'D') { // https://jvi.asm.org/content/jvi/73/9/7343.full.pdf page 7346
                // VPR-D-NKTATN  VPK-D-NKTATN
                replace_front_deletions(deletions, v_pos_num_t{{pos1_t{163}, 1}, {pos1_t{164}, 1}});
                special_deletions = true;
            }
            else if (at162 == 'K' && at163 == 'N') {
                // VPK--NNKTATN
                replace_front_deletions(deletions, v_pos_num_t{{pos1_t{163}, 2}});
                special_deletions = true;
            }
            else if (at162 == 'R' && (at163 == 'E' || at163 == 'N')) {
                // VPRE--NNKTATN VPRN--
                replace_front_deletions(deletions, v_pos_num_t{{pos1_t{164}, 2}});
                special_deletions = true;
            }
            else if ((at162 == 'N' || at162 == 'X') && at163 == 'K') {
                // VP-N-KNKTATNPLTI
                // VP-X-KNKTATNPLTI
                replace_front_deletions(deletions, v_pos_num_t{{pos1_t{162}, 1}, {pos1_t{163}, 1}});
                special_deletions = true;
            }
            else if (at162 == 'D' || at162 == 'G') { // VPDK -> VP--DK, VPDD -> VP--DD vic 2del
                // 2del
            }
            else {
                AD_DEBUG("del-spec {} {} {}", sequence.aa_aligned_without_deletions_substr(159, 14), deletions, sequence.name());
            }
        }
        else if (N_deletions_at(deletions, 1, pos1_t{159}, pos1_t{162}) && sequence.aa_aligned_without_deletions_substr(159, 4) == "VPRD") {
            // B/BEIJING/258/1993, B/NEW YORK/1044/2001, mistake in deletion detection, it's YAMAGATA
            replace_front_deletions(deletions, v_pos_num_t{{pos1_t{163}, 1}});
        }

        //---------- VICTORIA ----------

        if (no_deletions_after_before(deletions, pos1_t{10}, pos1_t{500})) { // may have deletions before 10 (e.g at the beginning due to trancation) and after 500
            set_lineage(sequence, acmacs::virus::lineage::VICTORIA, fasta_ref, "no");
        }
        else if (N_deletions_at(deletions, 2, pos1_t{162})) {
            // VICTORIA 2del 2017
            set_lineage(sequence, acmacs::virus::lineage::VICTORIA, fasta_ref, "victoria del2017");
        }
        else if (N_deletions_at(deletions, 3, pos1_t{162}, pos1_t{164})) {
            // VICTORIA triple del 2017
            // according to David Burke 2019-07-16 14:27
            // VPK---NKTAT
            deletions.deletions.front().pos = b_vic_del_mutants_pos;
            set_lineage(sequence, acmacs::virus::lineage::VICTORIA, fasta_ref, "victoria tripledel2017");
        }
        else if (N_deletions_at(deletions, 6, pos1_t{164})) {
            // VICTORIA sixdel2019 (only from Japan as of 2019-07-19)
            // David Burke 2019-07-19 15:40: These look really
            // unusual. Based on the geometry of the loop, I would
            // tend to align the N with C-terminal side: B/KANAGAWA/AC1867/2019 VPK------NTNP
            deletions.deletions.front().pos = b_vic_del_mutants_pos;
            set_lineage(sequence, acmacs::virus::lineage::VICTORIA, fasta_ref, "victoria sixdel2019 (pos shifted)");
        }

        //---------- YAMAGATA ----------

        else if (N_deletions_at(deletions, 1, pos1_t{163})) { // B/YAMAGATA/16/1988 has also del-164-1: VPR-D-NKTA
            set_lineage(sequence, acmacs::virus::lineage::YAMAGATA, fasta_ref, "yamagata");
        }
        else if (is_yamagata_shifted(sequence)) {
            // AD_DEBUG("{} {}", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
            set_lineage(sequence, acmacs::virus::lineage::YAMAGATA, fasta_ref, "yamagata-shifted");
            deletions.deletions = std::vector<acmacs::seqdb::scan::deletions_insertions_t::pos_num_t>{{pos1_t{163}, 1}};
        }
        else if (N_deletions_at(deletions, 2, pos1_t{163}) && sequence.year() <= 2013) {
            set_lineage(sequence, acmacs::virus::lineage::YAMAGATA, fasta_ref, "yamagata");
        }
        else if (N_deletions_at(deletions, 2, pos1_t{169})) {
            // 12 sequences from TAIWAN 2010 have deletions 169:2
            // sequence.lineage(acmacs::virus::lineage_t{});
        }
        else if (N_deletions_at(deletions, 1, pos1_t{160}) && no_deletions_after_before(deletions, pos1_t{161}, pos1_t{500}) && sequence.aa_at_pos(pos1_t{161}) == 'E' &&
                 sequence.aa_at_pos(pos1_t{163}) == 'K') {
            // deletion detection was invalid, most probably due to 162X. B/ALICANTE/19_0649/20171219
            // fmt::print(stderr, "DEBUG: {} 160{} 161{} 162{} 163{}\n", acmacs::seqdb::scan::format(deletions), sequence.aa_at_pos(pos1_t{160}),
            // sequence.aa_at_pos(pos1_t{161}), sequence.aa_at_pos(pos1_t{162}), sequence.aa_at_pos(pos1_t{163}));
            set_lineage(sequence, acmacs::virus::lineage::YAMAGATA, fasta_ref, "yamagata");
            deletions.deletions = std::vector<acmacs::seqdb::scan::deletions_insertions_t::pos_num_t>{{pos1_t{163}, 1}};
        }
        else if (is_semi_ignored(sequence)) {
            AD_INFO("{} {}", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
        }
        else if (is_ignored(sequence)) {
            // do not issue warning
        }
        else if (!special_deletions) {
            // AD_DEBUG("1-at-163:{} no-between-164-500:{}", N_deletions_at(deletions, 1, pos1_t{163}), no_deletions_after_before(deletions, pos1_t{164}, pos1_t{500}));
            AD_WARNING("not-special {} {} {} {}", sequence.aa_aligned_without_deletions_substr(159, 14), deletions, sequence.name(), fasta_ref); // sequence.aa_format());
        }

        if (sequence.lineage() != fasta_lineage && fasta_lineage != acmacs::virus::lineage::UNKNOWN && !fasta_lineage.empty())
            AD_WARNING("B-lineage fas:{} seq:{} {} {} {}", fasta_lineage[0], sequence.lineage().empty() ? ' ' : sequence.lineage()[0], sequence.aa_aligned_without_deletions_substr(159, 14), deletions, sequence.name());

        // if (sequence.name() == acmacs::virus::name_t{"B/BEIJING/15/1984"})
        //     AD_DEBUG("{} {} {}", sequence.name(), sequence.aa_aligned_substr(159, 14), deletions);
        // if (sequence.aa_aligned_substr(159, 4) == "VPRD" && sequence.lineage() != acmacs::virus::lineage::YAMAGATA)
        //     AD_WARNING("VPRD: {} {} fas:{:8s} seq:{:8s} {}", sequence.aa_aligned_substr(159, 14), deletions, fasta_lineage, sequence.lineage(), sequence.name());

    } // lineage

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
                ; // sequence.add_clade(acmacs::seqdb::clade_t{"*DEL"});
            else if (del1.pos == acmacs::seqdb::pos1_t{127} && del1.num == 1 && (year < 2018 || fasta_ref.find("seasonal") != std::string_view::npos))
                ; // sequence.add_clade(acmacs::seqdb::clade_t{"*DEL-127:1"});
            else if (del1.pos == acmacs::seqdb::pos1_t{160} && del1.num == 4 && sequence.name() == acmacs::virus::name_t{"A(H1N1)/NEWPORT/323/2019"})
                fmt::print(stderr, "INFO: {} {}\n", sequence.full_name(), acmacs::seqdb::scan::format(deletions));
            else if (del1.pos > acmacs::seqdb::pos1_t{400})
                ; // ignore
            else
                warn();
        }
        else if (deletions.deletions.size() > 1) {
            if (!host.empty() || year < 2010)
                ; // sequence.add_clade(acmacs::seqdb::clade_t{"*DEL"});
            else
                warn();
        }
        else if (!deletions.insertions.empty())
            ; // sequence.add_clade(acmacs::seqdb::clade_t{"*INS"});
        else if (!deletions.empty()) {
            warn();
        }

    } // deletions

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
            ; // sequence.add_clade(acmacs::seqdb::clade_t{"*INS"});
        else if (!deletions.empty()) {
            if (sequence.aa_aligned_length() < 500)
                ; // ignore short
            else if (!acmacs::virus::host(sequence.name()).empty())
                ; // ignore
            else if (sequence.year() < 2018)
                ; // sequence.add_clade(acmacs::seqdb::clade_t{"*DEL"});
            else
                warn();
        }

    } // deletions

} // namespace local::H3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
