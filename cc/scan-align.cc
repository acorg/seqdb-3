// #include <optional>
// #include <array>
// #include <map>
// #include <tuple>

#include "acmacs-base/named-type.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/fmt.hh"
#include "seqdb-3/scan-align.hh"
#include "seqdb-3/scan-fasta.hh"
#include "seqdb-3/scan-deletions.hh"
#include "seqdb-3/hamming-distance.hh"

// http://signalpeptide.com

namespace local
{
    // returns if it was aligned
    static bool align(acmacs::seqdb::v3::scan::sequence_t& sequence, const acmacs::virus::type_subtype_t& type_subtype_hint);
    static std::optional<std::tuple<int, acmacs::virus::type_subtype_t>> align(std::string_view amino_acids, const acmacs::virus::type_subtype_t& type_subtype_hint);

    // ----------------------------------------------------------------------

    // inline std::string_view type_subtype_hint(std::string_view type_subtype)
    // {
    //     return type_subtype.size() > 5 && type_subtype[0] == 'A' ? type_subtype.substr(0, 5) : type_subtype;
    // }

    inline bool has_infix(std::string_view source, size_t pos, std::string_view match)
    {
        return source.substr(pos, match.size()) == match;
    }

    inline std::string::size_type find_in_sequence(std::string_view sequence, size_t limit, std::initializer_list<const char*> look_for)
    {
        const auto source = sequence.substr(0, limit);
        for (const char* str : look_for) {
            if (const auto pos = source.find(str); pos != std::string::npos)
                return pos;
        }
        return std::string::npos;
    }

    struct start_aa_t
    {
        const char* type_subtype_h_or_b;
        char start_aa;
    };

    static constexpr const std::array start_aa_table{
        start_aa_t{"H1",  'D'}, // DTIC, DTLC
        start_aa_t{"H2",  'D'}, // DQIC
        start_aa_t{"H3",  'Q'},
        start_aa_t{"H4",  'Q'},
        start_aa_t{"H5",  'D'}, // DQIC
        start_aa_t{"H6",  'D'},
        start_aa_t{"H7",  'D'}, // DKIC
        start_aa_t{"H8",  'D'}, // DRIC
        start_aa_t{"H9",  'D'}, // DKIC
        start_aa_t{"H10", 'D'}, // DKIC
        start_aa_t{"H11", 'D'}, // DEIC
        start_aa_t{"H12", 'D'}, // DKIC
        start_aa_t{"H13", 'D'}, // DRIC
        start_aa_t{"H14", 'Q'}, // QITN
        start_aa_t{"H15", 'D'}, // DKIC
        start_aa_t{"H16", 'D'}, // DKIC
        start_aa_t{"H17", 'D'}, // DRIC
        start_aa_t{"B",     'D'}, // DRIC
    };

    inline char start_aa(const acmacs::virus::type_subtype_t& hint)
    {
        if (const auto found = ranges::find_if(start_aa_table, [hint](const auto& entry) { return entry.type_subtype_h_or_b == hint.h_or_b(); }); found != ranges::end(start_aa_table))
            return found->start_aa;
        throw std::runtime_error(fmt::format("local::start_aa: unsupported type_subtype: {}", hint));
    }

    // ----------------------------------------------------------------------

    class Aligner
    {
      public:
        Aligner() = default;

        void update(std::string_view amino_acids, const acmacs::virus::type_subtype_t& type_subtype)
        {
            tables_.try_emplace(std::string(type_subtype.h_or_b())).first->second.update(amino_acids, 0);
        }

        std::optional<std::tuple<int, acmacs::virus::type_subtype_t>> align(std::string_view amino_acids, const acmacs::virus::type_subtype_t& type_subtype_hint) const
        {
            if (const auto found = tables_.find(type_subtype_hint.h_or_b()); found != tables_.end()) {
                if (const auto res = found->second.align(start_aa(type_subtype_hint), amino_acids); res.has_value())
                    return std::tuple(*res, type_subtype_hint);
            }
            return std::nullopt;
        }

        void report() const
        {
            fmt::print(stderr, "Aligner {}\n", tables_.size());
            for (const auto& [type_subtype, table] : tables_)
                table.report(fmt::format(" {:<8s} ", type_subtype));
            fmt::print(stderr, "\n");
        }

      private:
        constexpr static size_t max_sequence_length{1000};
        constexpr static size_t number_of_symbols{128};
        constexpr static size_t table_size{number_of_symbols * max_sequence_length};

        struct table_t
        {
            using contribution_t = int;

            table_t()
            {
                ranges::fill(data, 1);
                for (auto pos : ranges::views::iota(0UL, max_sequence_length)) { // X and - do not contribute at any position
                    data[number_of_symbols * pos + static_cast<size_t>('X')] = 0;
                    data[number_of_symbols * pos + static_cast<size_t>('-')] = 0;
                }
            }

            void update(std::string_view amino_acids, size_t prefix_size)
            {
                for (char aa : amino_acids) {
                    data[number_of_symbols * prefix_size + static_cast<size_t>(aa)] = 0;
                    ++prefix_size;
                }
            }

            std::optional<int> align(char start_aa, std::string_view amino_acids) const;
            void report(std::string prefix) const;

            std::array<contribution_t, table_size> data;
        };

        std::map<std::string, table_t, std::less<>> tables_;
    };

} // namespace local

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::translate_align(std::vector<fasta::scan_result_t>& sequences)
{
#pragma omp parallel for default(shared) schedule(static, 256)
    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        auto& entry = sequences[e_no];
        entry.sequence.translate();
        local::align(entry.sequence, entry.fasta.type_subtype);
    }

    // remove not translated
    sequences.erase(std::remove_if(std::begin(sequences), std::end(sequences), [](const auto& entry) { return entry.sequence.aa().empty(); }), std::end(sequences));

    // fmt::print(stderr, "translate_align aligned 1: {} :: {}\n", ranges::count_if(sequences, fasta::is_aligned), date::current_date_time());

    local::Aligner aligner;
    for (const auto& entry : sequences | ranges::views::filter(fasta::is_aligned)) {
        aligner.update(entry.sequence.aa_aligned(), entry.sequence.type_subtype());
    }
    // aligner.report();

    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        auto& entry = sequences[e_no];
        if (!entry.sequence.aligned()) {
            if (const auto align_data = aligner.align(entry.sequence.aa(), entry.fasta.type_subtype); align_data.has_value()) {
                const auto [shift, type_subtype] = *align_data;
                entry.sequence.set_shift(shift, type_subtype);
            }
        }
    }

    // fmt::print(stderr, "translate_align aligned 2: {} :: {}\n", ranges::count_if(sequences, fasta::is_aligned), date::current_date_time());

} // acmacs::seqdb::v3::scan::translate_align_clade

// ----------------------------------------------------------------------

bool local::align(acmacs::seqdb::v3::scan::sequence_t& sequence, const acmacs::virus::type_subtype_t& type_subtype_hint)
{
    if (const auto shift_type = align(sequence.aa(), type_subtype_hint); shift_type.has_value()) {
        const auto [shift, type_subtype] = *shift_type;
        sequence.set_shift(shift, type_subtype);
        return true;
    }
    else
        return false;

} // acmacs::seqdb::v3::scan::align

// ----------------------------------------------------------------------

std::optional<std::tuple<int, acmacs::virus::type_subtype_t>> local::align(std::string_view amino_acids, const acmacs::virus::type_subtype_t& type_subtype_hint)
{
    const auto make_type_subtype = [&type_subtype_hint](const char* detected_type_subtype) -> acmacs::virus::type_subtype_t {
        const auto dts = acmacs::virus::type_subtype_t{detected_type_subtype};
        if (type_subtype_hint.h_or_b() == dts.h_or_b())
            return type_subtype_hint;
        else
            return dts;
    };

    // --------------------------------------------------
    // first stage

    // H3
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MKTII"});
        pos != std::string::npos &&
        (amino_acids[pos + 16] == 'Q' || amino_acids[pos + 15] == 'A')) // amino_acids.substr(pos + 15, 2) != "DR") { // DR[ISV]C - start of the B sequence (signal peptide is 15 aas!)
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H3)")};

    // H1
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MKV", "MKA", "MEA", "MEV"});
        pos != std::string::npos && (has_infix(amino_acids, pos + 17, "DTLC") || has_infix(amino_acids, pos + 17, "DTIC")))
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H1)")};

    // B
    {
        // Only B has CTDL at first 100 AAs
        if (const auto pos = find_in_sequence(amino_acids, 100, {"CTDL"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 59, make_type_subtype("B")};
        // Only B has NSPHVV at first 100 AAs
        if (const auto pos = find_in_sequence(amino_acids, 100, {"NSPHVV"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 10, make_type_subtype("B")};
        // B specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"EHIRL"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 114, make_type_subtype("B")};
        // Only B (YAMAGATA?) has CPNATS in whole AA sequence
        if (const auto pos = find_in_sequence(amino_acids, 250, {"CPNATS"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 142, make_type_subtype("B")};
        if (const auto pos = find_in_sequence(amino_acids, 250, {"PNATSK"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 143, make_type_subtype("B")};
        // B specific (VICTORIA?)
        if (const auto pos = find_in_sequence(amino_acids, 150, {"NVTNG"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 144, make_type_subtype("B")};
        // --- below are redundant 2019-06-01
        // Only B has DRICT
        // if (const auto pos = find_in_sequence(amino_acids, 50, {"DRICT"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos), make_type_subtype("B")};
        // Only B has MKAIIVL at first 100 AAs
        // if (const auto pos = find_in_sequence(amino_acids, 50, {"MKAIIVL"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) + 15, make_type_subtype("B")};
        // // Only B has NVTGV at first 100 AAs
        // if (const auto pos = find_in_sequence(amino_acids, 50, {"NVTGV"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) - 24, make_type_subtype("B")};
    }

    // H2
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MTIT", "MAII"}); pos != std::string::npos && has_infix(amino_acids, pos + 14, "GDQIC"))
        return std::tuple{static_cast<int>(pos) + 15, make_type_subtype("A(H2)")};

    // H4
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MLS"}); pos != std::string::npos && (amino_acids[pos + 16] == 'Q' || has_infix(amino_acids, pos + 16, "SQNY")))
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H4)")};

    // H5
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MEKIV", "MERIV"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H5)")};

    // H6
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MIAIIV", "MIAIII"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H6)")};

    // H7
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MNIQ", "MNNQ", "MNTQ"});
        pos != std::string::npos && amino_acids[pos + 17] != 'S' && has_infix(amino_acids, pos + 18, "DKIC")) // SDKIC is H15 most probably
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H7)")};

    // H8
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MEKFIA"}); pos != std::string::npos && amino_acids[pos + 17] == 'D')
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H8)")};

    // H9
    if (const auto pos = find_in_sequence(amino_acids, 20, {"METIS", "MEIIS", "MEV"}); pos != std::string::npos && has_infix(amino_acids, pos + 17, "ADKIC"))
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H9)")};

    // H10
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MYK"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H10)")};

    // H11
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MK"}); pos != std::string::npos && has_infix(amino_acids, pos + 16, "DEIC"))
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H11)")};

    // H12
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MEK"}); pos != std::string::npos && has_infix(amino_acids, pos + 15, "AYDKIC"))
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H12)")};

    // H13
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MDI", "MAL", "MEV"}); pos != std::string::npos && has_infix(amino_acids, pos + 17, "ADRIC"))
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H13)")};

    // H14
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MIA"}); pos != std::string::npos && has_infix(amino_acids, pos + 14, "AYSQITN"))
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H14)")};

    // H15 - second stage only

    // H16
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MMVK", "MMIK"}); pos != std::string::npos && has_infix(amino_acids, pos + 19, "DKIC"))
        return std::tuple{static_cast<int>(pos) + 19, make_type_subtype("A(H16)")};

    // H17
    if (const auto pos = find_in_sequence(amino_acids, 20, {"MEL"}); pos != std::string::npos && has_infix(amino_acids, pos + 17, "GDRICI"))
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H17)")};

    // --------------------------------------------------
    // second stage

    // H4
    if (const auto pos = find_in_sequence(amino_acids, 100, {"QNYT"}); pos != std::string::npos && has_infix(amino_acids, pos + 11, "GHHA"))
        return std::tuple{static_cast<int>(pos), make_type_subtype("A(H4)")};

    // H11 (DEICIGYL is specific)
    if (const auto pos = find_in_sequence(amino_acids, 50, {"DEICIGYL"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos), make_type_subtype("A(H11)")};

    // H15
    if (const auto pos = find_in_sequence(amino_acids, 100, {"KSDKICLGHHA"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 2, make_type_subtype("A(H15)")};

    // --------------------------------------------------
    // third stage

    // H3
    {
        // Only H3 (and H0N0) has CTLID in the whole AA sequence
        if (const auto pos = find_in_sequence(amino_acids, 150, {"CTLID", "CTLMDALL", "CTLVD"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 63, make_type_subtype("A(H3)")};
        // Only H3 (and H0N0) has PNGTIVKTI in the whole AA sequence
        if (const auto pos = find_in_sequence(amino_acids, 100, {"PNGTIVKTI"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 20, make_type_subtype("A(H3)")};
        // Only H3 (and H0N0) has DKLYIWG in the whole AA sequence
        if (const auto pos = find_in_sequence(amino_acids, 200, {"DKLYIWG"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 174, make_type_subtype("A(H3)")};
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"SNCYPYDV"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 94, make_type_subtype("A(H3)")};
        // -- below is redundant
        // // Only H3 (and H0N0) has SNCY in the whole AA sequence
        // if (const auto pos = find_in_sequence(amino_acids, 150, {"SNCY"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) - 94, make_type_subtype("A(H3)")};
        // // Only H3 (and H0N0) has WTGVT in the whole AA sequence
        // if (const auto pos = find_in_sequence(amino_acids, 200, {"WTGVT"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) - 126, make_type_subtype("A(H3)")};
        // // Only H3 (and H0N0) has NWLTH in the whole AA sequence
        // if (const auto pos = find_in_sequence(amino_acids, 200, {"NWLTH"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) - 151, make_type_subtype("A(H3)")};
    }

    // H1
    {
        // VLEKN is H1 specific (whole AA sequence)
        if (const auto pos = find_in_sequence(amino_acids, 50, {"VLEKN"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 18, make_type_subtype("A(H1)")};
        // SSWSYI and ESWSYI are H1 specific (whole AA sequence)
        if (const auto pos = find_in_sequence(amino_acids, 150, {"SSWSYI", "ESWSYI"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 73, make_type_subtype("A(H1)")};
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"FERFEI"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 110, make_type_subtype("A(H1)")};
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 200, {"IWLVKKG"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 148, make_type_subtype("A(H1)")};
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 200, {"SSVSSF"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 105, make_type_subtype("A(H1)")};
        // (redundant) specific
        // if (const auto pos = find_in_sequence(amino_acids, 200, {"GVTAACPH"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) - 130, make_type_subtype("A(H1)")};
    }

    // H5
    {
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 100, {"GVKPLIL", "GVRPLIL"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 45, make_type_subtype("A(H5)")};
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 100, {"GWLLGNPMCDE"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 58, make_type_subtype("A(H5)")};
        // specific at first 150
        if (const auto pos = find_in_sequence(amino_acids, 150, {"NHFE"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 108, make_type_subtype("A(H5)")};
        // (redundant) specific
        // if (const auto pos = find_in_sequence(amino_acids, 50, {"GYHANNSTEQ"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) - 5, make_type_subtype("A(H5)")};
    }

    // H6
    {
        // QKEER is H6 specific
        if (const auto pos = find_in_sequence(amino_acids, 100, {"QKEER"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 35, make_type_subtype("A(H6)")};
        // EELKA is H6 specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"EELKA"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 98, make_type_subtype("A(H6)")};
    }

    // H7
    {
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 100, {"GQCGL"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 51, make_type_subtype("A(H7)")};
    }

    // H8
    {
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 200, {"FYRSINWL"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 141, make_type_subtype("A(H8)")};
    }

    // H9
    {
        // QSTN is H9 specific
        if (const auto pos = find_in_sequence(amino_acids, 50, {"QSTN"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 7, make_type_subtype("A(H9)")};
        // CDLLLGG, CDLLLEG are H9 specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"CDLLLGG", "CDLLLEG"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 66, make_type_subtype("A(H9)")};
        // LEELRS is H9 specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"LEELRS"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 97, make_type_subtype("A(H9)")};
        // SARSYQ is H9 specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"SARSYQ"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 106, make_type_subtype("A(H9)")};
        // SSYQRIQ is H9 specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"SSYQRIQ"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 108, make_type_subtype("A(H9)")};
    }

    // H10
    {
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 50, {"NGTIVKTLTNE"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 11, make_type_subtype("A(H10)")};
        // specific
        if (const auto pos = find_in_sequence(amino_acids, 150, {"QKIMESG"}); pos != std::string::npos)
            return std::tuple{static_cast<int>(pos) - 99, make_type_subtype("A(H10)")};
        // (redundant) specific
        // if (const auto pos = find_in_sequence(amino_acids, 50, {"LDKICLGHHA"}); pos != std::string::npos)
        //     return std::tuple{static_cast<int>(pos) + 1, make_type_subtype("A(H10)")};
    }

    // H11 (SSVEL is specific)
    if (const auto pos = find_in_sequence(amino_acids, 100, {"SSVEL"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) - 27, make_type_subtype("A(H11)")};

    // H13 (specific)
    if (const auto pos = find_in_sequence(amino_acids, 50, {"VGYLSTN"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) - 4, make_type_subtype("A(H13)")};

    // H16 (specific)
    if (const auto pos = find_in_sequence(amino_acids, 70, {"DTLTENGVP", "DTLIENGVP"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) - 16, make_type_subtype("A(H16)")};

    return std::nullopt;

} // acmacs::seqdb::v3::scan::align

// ----------------------------------------------------------------------

std::optional<int> local::Aligner::table_t::align(char start_aa, std::string_view amino_acids) const
{
    // fmt::print(stderr, "Aligner::table_t::align {}\n{}\n", debug_name, amino_acids);
    for (auto p_start = amino_acids.find(start_aa); p_start < (amino_acids.size() / 2); p_start = amino_acids.find(start_aa, p_start + 1)) {
        const auto all_pos = ranges::views::iota(0UL, std::min(max_sequence_length, amino_acids.size() - p_start));
        if (const auto failed_pos = ranges::find_if(all_pos, [this,amino_acids,p_start](size_t pos) -> bool { return data[number_of_symbols * pos + static_cast<size_t>(amino_acids[p_start + pos])]; }); failed_pos != ranges::end(all_pos)) {
            // if (*failed_pos > 10)
            //     fmt::print(stderr, "Aligner::table_t::align FAILED: shift:{} {}:{} -- {} -- {}\n", p_start, *failed_pos, amino_acids[p_start + *failed_pos], debug_name, amino_acids.substr(0, p_start + *failed_pos + 5));
        }
        else {
            // fmt::print(stderr, "Aligner::table_t::align good shift:{}\n", p_start);
            return static_cast<int>(p_start);
        }
    }

    return std::nullopt;

} // local::Aligner::table_t::align

// ----------------------------------------------------------------------

void local::Aligner::table_t::report(std::string prefix) const
{
    using iter_t = decltype(data.begin());
    const auto increment = [](iter_t first, iter_t last, iter_t value) -> iter_t {
        while (value != last) {
            ++value;
            if (*value == 0 && (value - first) != static_cast<ssize_t>('X') /* && (value - first) != static_cast<ssize_t>('-') */)
                break;
        }
        return value;
    };
    const auto begin = [this](size_t pos) -> iter_t { return data.begin() + number_of_symbols * pos; };
    const auto end = [begin](size_t pos) -> iter_t { return begin(pos) + static_cast<ssize_t>('Z' + 1); };

    std::array<iter_t, max_sequence_length> iters;
    std::array<bool, max_sequence_length> completed;
    size_t last_pos = 0;
    for (auto pos : ranges::views::iota(0UL, max_sequence_length)) {
        iters[pos] = increment(begin(pos), end(pos), begin(pos) + static_cast<ssize_t>('A') - 1);
        if (iters[pos] == end(pos)) {
            completed[pos] = true;
        }
        else {
            completed[pos] = false;
            last_pos = pos + 1;
        }
    }

    const auto print_line = [&iters, &completed, increment, begin, end, last_pos]() {
        for (auto pos : ranges::views::iota(0UL, last_pos)) {
            if (iters[pos] != end(pos)) {
                fmt::print(stderr, "{}", static_cast<char>(iters[pos] - begin(pos)));
                iters[pos] = increment(begin(pos), end(pos), iters[pos]);
                completed[pos] = iters[pos] == end(pos);
            }
            else
                fmt::print(stderr, " ");
        }
        fmt::print(stderr, "\n");
    };

    fmt::print(stderr, "{}", prefix);
    print_line();
    const std::string prefix_space(prefix.size(), ' ');
    while (!ranges::all_of(completed, [](auto val) { return val; })) {
        fmt::print(stderr, "{}", prefix_space);
        print_line();
    }

} // local::Aligner::table_t::report

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
