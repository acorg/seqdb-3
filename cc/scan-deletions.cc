#include <numeric>
#include <optional>

#include "acmacs-base/counter.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/date.hh"
#include "acmacs-base/range-v3.hh"
#include "seqdb-3/scan-deletions.hh"
#include "seqdb-3/hamming-distance.hh"

// ----------------------------------------------------------------------

namespace local
{
    struct not_verified : public std::runtime_error { using std::runtime_error::runtime_error; };

    using subtype_master_t = std::map<std::string, const acmacs::seqdb::scan::sequence_t*, std::less<>>;

    inline const acmacs::seqdb::scan::sequence_t* find_master(std::string_view subtype, const subtype_master_t& masters)
    {
        if (const auto it = masters.find(subtype); it != masters.end())
            return it->second;
        else
            return nullptr;
    }

    subtype_master_t masters_per_subtype(const std::vector<acmacs::seqdb::v3::scan::fasta::scan_result_t>& sequences);

#include "acmacs-base/global-constructors-push.hh"
    static const std::array master_sequences_for_insertions = {
        std::pair{std::string{"B"},  acmacs::seqdb::scan::sequence_t::from_aligned_aa(acmacs::virus::name_t{"B/BRISBANE/60/2008 VICTORIA (master_sequences_for_insertions)"},        "DRICTGITSSNSPHVVKTATQGEVNVTGVIPLTTTPTKSHFANLKGTETRGKLCPKCLNCTDLDVALGRPKCTGKIPSARVSILHEVRPVTSGCFPIMHDRTKIRQLPNLLRGYEHIRLSTHNVINAENAPGGPYKIGTSGSCPNITNGNGFFATMAWAVPKNDKNKTATNPLTIEVPYICTEGEDQITVWGFHSDNETQMAKLYGDSKPQKFTSSANGVTTHYVSQIGGFPNQTEDGGLPQSGRIVVDYMVQKSGKTGTITYQRGILLPQKVWCASGRSKVIKGSLPLIGEADCLHEKYGGLNKSKPYYTGEHAKAIGNCPIWVKTPLKLANGTKYRPPAKLLKERGFFGAIAGFLEGGWEGMIAGWHGYTSHGAHGVAVAADLKSTQEAINKITKNLNSLSELEVKNLQRLSGAMDELHNEILELDEKVDDLRADTISSQIELAVLLSNEGIINSEDEHLLALERKLKKMLGPSAVEIGNGCFETKHKCNQTCLDRIAAGTFDAGEFSLPTFDSLNITAASLNDDGLDNHTILLYYSTAASSLAVTLMIAIFVVYMVSRDNVSCSICL")},
        std::pair{std::string{"H1"}, acmacs::seqdb::scan::sequence_t::from_aligned_aa(acmacs::virus::name_t{"A(H1N1)/CALIFORNIA/7/2009 (master_sequences_for_insertions)"},          "DTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI")},
        std::pair{std::string{"H3"}, acmacs::seqdb::scan::sequence_t::from_aligned_aa(acmacs::virus::name_t{"A(H3N2)/HONG_KONG/1/1968_h1EAB4576 (master_sequences_for_insertions)"}, "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI")},
    };
#include "acmacs-base/diagnostics-pop.hh"

    inline bool is_whocc_subtype(const acmacs::virus::type_subtype_t& subtype)
    {
        const auto hn = subtype.hn_or_b();
        return hn == "B" || hn == "H3N2" || hn == "H1N1";
    }
}

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::detect_insertions_deletions(std::vector<fasta::scan_result_t>& sequence_data)
{
    const auto masters = local::masters_per_subtype(sequence_data);
    // fmt::print(stderr, "masters_per_subtype {}\n", masters.size());

// #pragma omp parallel for default(shared) schedule(static) // , 100)
    for (auto sc_p = sequence_data.begin(); sc_p != sequence_data.end(); ++sc_p) {
        if (!sc_p->reference && sc_p->sequence.aligned()) { //  && sc_p->sequence.type_subtype() == acmacs::virus::type_subtype_t{"B"}) {
            if (const auto* master = local::find_master(sc_p->sequence.type_subtype().h_or_b(), masters); master && master != &sc_p->sequence) {
                // AD_DEBUG("dels {}", sc_p->sequence.name());
                deletions_insertions(*master, sc_p->sequence);
            }
            else
                AD_WARNING(local::is_whocc_subtype(sc_p->sequence.type_subtype()), "no master for {}", sc_p->sequence.name());
        }
        // else
        //     AD_DEBUG("not aligned? {}", sc_p->sequence.name());
    }

} // detect_insertions_deletions

// ----------------------------------------------------------------------

local::subtype_master_t local::masters_per_subtype(const std::vector<acmacs::seqdb::v3::scan::fasta::scan_result_t>& sequences)
{
    std::map<std::string, acmacs::Counter<size_t>> aligned_lengths;
    for (const auto& sc : sequences | ranges::views::filter(acmacs::seqdb::v3::scan::fasta::is_good))
        aligned_lengths.try_emplace(std::string(sc.sequence.type_subtype().h_or_b())).first->second.count(sc.sequence.aa_aligned_length());

    subtype_master_t masters;
    for (const auto& [subtype, counter] : aligned_lengths) {
        // AD_DEBUG("masters_per_subtype {} ({})", subtype, counter.total());
        const acmacs::seqdb::v3::scan::sequence_t* master = nullptr;
        if (const auto found = std::find_if(std::begin(master_sequences_for_insertions), std::end(master_sequences_for_insertions), [subtype=subtype](const auto& en) { return en.first == subtype; }); found != std::end(master_sequences_for_insertions)) {
            master = &found->second;
        }
        else {
            const size_t threshold = counter.total() / 6;
            size_t master_length = 0;
            for (const auto& vt : counter.counter()) {
                if (vt.second > threshold)
                    master_length = vt.first;
            }
            size_t num_X = 0;
            for (const auto& sc : sequences | ranges::views::filter(acmacs::seqdb::v3::scan::fasta::is_good)) {
                if (sc.sequence.type_subtype().h_or_b() == subtype && sc.sequence.aa_aligned_length() == master_length) {
                    if (master == nullptr) {
                        master = &sc.sequence;
                        num_X = sc.sequence.aa_number_of_X();
                    }
                    else if (const auto seq_X = sc.sequence.aa_number_of_X(); seq_X < num_X) {
                        master = &sc.sequence;
                        num_X = seq_X;
                    }
                    if (num_X == 0)
                        break;
                }
            }
        }
        if (master == nullptr)
            throw std::runtime_error(AD_FORMAT("internal (no master for {}) in acmacs::seqdb::v3::scan::masters_per_subtype", subtype));
        masters.emplace(subtype, master);
        // AD_DEBUG("master_per_subtype {}: {} aas:{}", subtype, master->name(), master->aa_aligned_length());
    }

    return masters;

} // local::masters_per_subtype

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::deletions_insertions(const sequence_t& master, sequence_t& to_align)
{
    const acmacs::debug dbg = acmacs::debug::no;
    // const acmacs::debug dbg = local::is_whocc_subtype(to_align.type_subtype()) ? acmacs::debug::yes : acmacs::debug::no;
    // AD_DEBUG("master: {}  to-align: {}", master.name(), to_align.name());

    auto master_aa_aligned = master.aa_aligned(), to_align_aa_aligned = to_align.aa_aligned();
    try {
        to_align.deletions() = deletions_insertions(master_aa_aligned, to_align_aa_aligned, dbg);
    }
    catch (local::not_verified& err) {
        if (local::is_whocc_subtype(to_align.type_subtype())) {
            AD_WARNING("deletions_insertions NOT VERIFIED  master: \"{}\"   to-align: \"{}\"  err: {}", master.name(), to_align.name(), err.what());
            // try {
            //     deletions_insertions(master_aa_aligned, to_align_aa_aligned, acmacs::debug::yes);
            // }
            // catch (local::not_verified&) {
            // }
        }
    }
    if (local::is_whocc_subtype(to_align.type_subtype())) {
        if (master_aa_aligned[0] != to_align_aa_aligned[0] && to_align_aa_aligned[0] != 'X') {
            to_align.add_issue(sequence::issue::garbage_at_the_beginning);
            // AD_DEBUG(to_align.type_subtype() == acmacs::virus::type_subtype_t{"A(H3N2)"} ? acmacs::debug::yes : acmacs::debug::no, "{} NOT aligned\n{}", to_align.name(), to_align_aa_aligned);
            // AD_DEBUG("{} NOT aligned\n{}", to_align.name(), to_align_aa_aligned);
        }
        if (master.aa_aligned_length() > to_align.aa_aligned_length()) {
            to_align.add_issue(sequence::issue::too_short);
            // AD_DEBUG("too_short {} {} < {}\n{}\n{}", to_align.name(), to_align.aa_aligned_length(), master.aa_aligned_length(), master.aa_format(), to_align.aa_format());
        }
        else {
            const auto master_with_deletions = master.aa_format(), to_align_with_deletions = to_align.aa_format();
            constexpr const size_t tail_size = 10, mismatches_threshold = 2;
            if (hamming_distance(std::string_view{master_with_deletions}.substr(master_with_deletions.size() - tail_size),
                                 std::string_view{to_align_with_deletions}.substr(master_with_deletions.size() - tail_size)) > mismatches_threshold) {
                // AD_DEBUG("garbage_at_the_end {}\n{}\n{}", to_align.name(), master.aa_format(), to_align.aa_format());
                to_align.add_issue(sequence::issue::garbage_at_the_end);
            }
        }
    }
    AD_DEBUG(dbg, "deletions: {}", to_align.deletions());
    AD_PRINT_IF(dbg, "\n");

} // acmacs::seqdb::v3::scan::deletions_insertions

// ----------------------------------------------------------------------

namespace local
{
    constexpr const ssize_t common_threshold = 3; // assume the chunk is common after that number of consecutive common positions
    // constexpr const ssize_t not_common_threshold = 10; // assume rest is different after that number of consecutive different positions
    constexpr const size_t max_deletions_insertions = 200; // give up if this number of deletions/insertions does not help
    constexpr double verify_threshold = 0.6;              // if number of common is less than this fraction of non-X in shortest of to_align and master sequences, verification fails

    static constexpr const auto are_common = [](char a, char b) -> bool { return a == b && a != 'X' && a != '-'; };

    struct find_head_t
    {
        size_t head = 0;
        size_t common = 0;
    };

    inline std::string format(const find_head_t& head)
    {
        return fmt::format("head:{} common:{}", head.head, head.common);
    }

    template <typename Iter> find_head_t find_head(const Iter first1, const Iter last1, Iter first2, const Iter last2, acmacs::debug dbg)
    {
        // find the last part with common_f()==true that is not shorter than common_threshold
        // returns offset of the end of this part
        auto f1 = first1, common_start = last1, last_common_end = first1;
        size_t common = 0, common_at_last_common_end = 0, really_common_in_this_common_chunk = 0;
        const auto update_last_common_end = [&]() {
            if (common_start != last1 && really_common_in_this_common_chunk >= common_threshold) {
                last_common_end = f1;
                // fmt::print(stderr, "last_common_end:{}\n", f1 - first1);
                common_at_last_common_end = common;
            }
        };

        for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
            if (*f1 == *first2 || *f1 == 'X' || *first2 == 'X') {
                if (*f1 == *first2) {
                    ++common;
                    ++really_common_in_this_common_chunk;
                }
                if (common_start == last1)
                    common_start = f1;
                // AD_DEBUG(dbg, "common:{} common_start:{}", f1 - first1, common_start - first1);
            }
            else {
                // AD_DEBUG(dbg, "NOTcommon:{}", f1 - first1);
                update_last_common_end();
                // if (static_cast<size_t>(f1 - last_common_end) > not_common_threshold)
                //     break; // too many different, stop searching
                common_start = last1;
                really_common_in_this_common_chunk = 0;
            }
        }
        update_last_common_end();

        // AD_DEBUG(dbg, "find_head end last_common_end:{} common_at_last_common_end:{}", last_common_end - first1, common_at_last_common_end);
        if (const auto head = static_cast<size_t>(last_common_end - first1); common_at_last_common_end * 3 > head) {
            AD_DEBUG(dbg, "head: {}  common_at_last_common_end: {}", head, common_at_last_common_end);
            return {head, common_at_last_common_end};
        }
        else {
            AD_DEBUG(dbg, "too few common in the head, try more deletions");
            return {0, 0}; // too few common in the head, try more deletions
        }
    }

    inline find_head_t find_common_head(std::string_view s1, std::string_view s2, acmacs::debug dbg)
    {
        return find_head(s1.begin(), s1.end(), s2.begin(), s2.end(), dbg);
    }

    // inline size_t find_uncommon_head(std::string_view s1, std::string_view s2)
    // {
    //     auto p1 = s1.begin(), p2 = s2.begin();
    //     for (; p1 != s1.end() && p2 != s2.end(); ++p1, ++p2) {
    //         if (are_common(*p1, *p2))
    //             break;
    //     }
    //     return static_cast<size_t>(p1 - s1.begin());
    // }

    struct deletions_insertions_at_start_t
    {
        size_t deletions = 0;
        size_t insertions = 0;
        find_head_t head;
    };

    inline deletions_insertions_at_start_t deletions_insertions_at_start(std::string_view master, std::string_view to_align, acmacs::debug dbg)
    {
        deletions_insertions_at_start_t result;
        for (size_t dels = 1; dels < max_deletions_insertions; ++dels) {
            if (dels < master.size()) {
                result.head = find_common_head(master.substr(dels), to_align, acmacs::debug::no);
                AD_DEBUG(dbg, "dels:{} {}\n{}\n{}", dels, format(result.head), master.substr(dels), to_align);
                if (result.head.head > common_threshold) {
                    result.deletions = dels;
                    break;
                }
            }
            if (dels < to_align.size()) {
                result.head = find_common_head(master, to_align.substr(dels), acmacs::debug::no);
                if (result.head.head > common_threshold) {
                    result.insertions = dels;
                    break;
                }
            }
        }
        return result;
    }

    inline size_t number_of_common(std::string_view s1, std::string_view s2)
    {
        size_t common = 0;
        for (auto p1 = s1.begin(), p2 = s2.begin(); p1 != s1.end() && p2 != s2.end(); ++p1, ++p2) {
            if (are_common(*p1, *p2))
                ++common;
        }
        return common;
    }

    inline size_t number_of_common(std::string_view master, std::string_view to_align, acmacs::seqdb::v3::scan::deletions_insertions_t& deletions)
    {
        return number_of_common(acmacs::seqdb::v3::scan::format_aa(deletions.insertions, master), acmacs::seqdb::v3::scan::format_aa(deletions.deletions, to_align));
    }

}

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::deletions_insertions_t acmacs::seqdb::v3::scan::deletions_insertions(std::string_view master, std::string_view to_align, acmacs::debug dbg)
{
    AD_DEBUG(dbg, "initial:\n{}\n{}", master, to_align);
    deletions_insertions_t deletions;
    const auto initial_head = local::find_common_head(master, to_align, dbg);
    size_t master_offset = initial_head.head, to_align_offset = initial_head.head;
    std::string_view master_tail = master.substr(master_offset), to_align_tail = to_align.substr(to_align_offset);
    AD_DEBUG(dbg, "{} number_of_common:{}\n{}\n{}", local::format(initial_head), local::number_of_common(master.substr(0, initial_head.head), to_align.substr(0, initial_head.head)), master.substr(0, initial_head.head), to_align.substr(0, initial_head.head));

    const auto update_both = [](size_t dels, size_t head, std::string_view& s1, std::string_view& s2, size_t& offset1, size_t& offset2) {
        s1.remove_prefix(head + dels);
        s2.remove_prefix(head);
        offset1 += head + dels;
        offset2 += head;
    };

    size_t common = initial_head.common;
    while (!master_tail.empty() && !to_align_tail.empty()) {
        AD_DEBUG(dbg, "m-offset:{} a-offset:{} common:{}\n{}\n{}", master_offset, to_align_offset, common, master_tail, to_align_tail);
        const auto tail_deletions = local::deletions_insertions_at_start(master_tail, to_align_tail, acmacs::debug::no /* dbg */);
        AD_DEBUG(dbg, "dels:{} ins:{} {} number_of_common:{}", tail_deletions.deletions, tail_deletions.insertions, local::format(tail_deletions.head), local::number_of_common(master_tail.substr(tail_deletions.deletions, tail_deletions.head.head), to_align_tail.substr(tail_deletions.insertions, tail_deletions.head.head)));
        if (tail_deletions.head.head == 0) { // < local::common_threshold) {
            common += local::number_of_common(master_tail, to_align_tail); // to avoid common diff warning below in case tail contains common aas
            break; // tails are different, insertions/deletions do not help
        }
        if (tail_deletions.deletions) {
            deletions.deletions.push_back({pos0_t{to_align_offset}, tail_deletions.deletions});
            update_both(tail_deletions.deletions, tail_deletions.head.head, master_tail, to_align_tail, master_offset, to_align_offset);
        }
        else { // insertions or nothing (in some cases)
            if (tail_deletions.insertions)
                deletions.insertions.push_back({pos0_t{master_offset}, tail_deletions.insertions});
            update_both(tail_deletions.insertions, tail_deletions.head.head, to_align_tail, master_tail, to_align_offset, master_offset);
        }
        common += tail_deletions.head.common;
    }

    // // sanity check (remove)
    // if (const auto nc = local::number_of_common(master, to_align, deletions); nc != common)
    //     fmt::print(stderr, "common diff: {} vs. number_of_common:{} {}\n{}\n{}\n", common, nc, format(deletions), acmacs::seqdb::v3::scan::format(deletions.insertions, master, '.'), acmacs::seqdb::v3::scan::format(deletions.deletions, to_align, '.'));

    // verify
    const auto get_num_non_x = [](std::string_view seq) { return seq.size() - static_cast<size_t>(std::count(std::begin(seq), std::end(seq), 'X')); };
    const auto num_common_threshold = static_cast<double>(master.size() < to_align.size() ? get_num_non_x(master) : get_num_non_x(to_align)) * local::verify_threshold;
    if (static_cast<double>(common) < num_common_threshold) {
        throw local::not_verified(fmt::format("common:{} vs size:{} num_common_threshold:{:.2f}\n{}\n{}\n{}\n{}\n",
                                              common, to_align.size(), num_common_threshold, master, to_align,
                                              acmacs::seqdb::v3::scan::format_aa(deletions.insertions, master, '.'), acmacs::seqdb::v3::scan::format_aa(deletions.deletions, to_align, '.')));
    }

    return deletions;

} // acmacs::seqdb::v3::scan::deletions_insertions

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
