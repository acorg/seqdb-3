#include <numeric>

#include "acmacs-base/counter.hh"
#include "acmacs-base/fmt.hh"
#include "seqdb-3/insertions.hh"

// ----------------------------------------------------------------------

namespace local
{
    // static size_t detect_master_length(const std::vector<std::reference_wrapper<acmacs::seqdb::sequence_t>>& sequences);
    // static const acmacs::seqdb::sequence_t& get_master(const std::vector<std::reference_wrapper<acmacs::seqdb::sequence_t>>& sequences);
    // static void detect_deletions(acmacs::seqdb::sequence_t& sequence, const std::string_view master);

    static constexpr inline bool common(char a, char b) { return a == b && a != 'X' && a != '-'; }

    static inline size_t number_of_common(std::string_view s1, std::string_view s2)
    {
        size_t num = 0;
        for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
            if (common(*f1, *f2))
                ++num;
        }
        return num;
    }

    // static inline size_t number_of_common_before(std::string_view s1, std::string_view s2, size_t last) { return number_of_common(s1.substr(0, last), s2.substr(0, last)); }

} // namespace local

// ----------------------------------------------------------------------

std::map<std::string, const acmacs::seqdb::v3::sequence_t*> acmacs::seqdb::v3::masters_per_subtype(const std::vector<fasta::scan_result_t>& sequences)
{
    std::map<std::string, acmacs::Counter<size_t>> aligned_lengths;
    for (const auto& sc : sequences)
        aligned_lengths.try_emplace(std::string(sc.sequence.type_subtype().h_or_b())).first->second.count(sc.sequence.aa_aligned_length());

    std::map<std::string, const sequence_t*> masters;
    for (const auto& [subtype, counter] : aligned_lengths) {
        const size_t threshold = counter.total() / 6;
        size_t master_length = 0;
        for (const auto& vt : counter.counter()) {
            if (vt.second > threshold)
                master_length = vt.first;
        }
        for (const auto& sc : sequences) {
            if (sc.sequence.type_subtype().h_or_b() == subtype && sc.sequence.aa_aligned_length() == master_length) {
                masters.emplace(subtype, &sc.sequence);
                break;
            }
        }
    }

    return masters;

} // acmacs::seqdb::v3::masters_per_subtype

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::insertions_deletions(std::vector<std::reference_wrapper<seqdb::sequence_t>>& sequences)
{
    // const auto& master = local::get_master(sequences);
    // const auto master_seq = master.aa_aligned_fast();
    // for (auto& seq : sequences) {
    //     if (&seq.get() != &master)
    //         detect_deletions(seq.get(), master_seq);
    // }

} // acmacs::seqdb::v3::insertions_deletions

// ----------------------------------------------------------------------

// const acmacs::seqdb::sequence_t& local::get_master(const std::vector<std::reference_wrapper<acmacs::seqdb::sequence_t>>& sequences)
// {
//     const auto master_length = detect_master_length(sequences);
//     for (const auto& seq : sequences) {
//         if (seq.get().aa_aligned_length() == master_length)
//             return seq.get();
//     }
//     throw std::runtime_error("internal in insertions.cc get_master");

// } // get_master

// // ----------------------------------------------------------------------

// size_t local::detect_master_length(const std::vector<std::reference_wrapper<acmacs::seqdb::sequence_t>>& sequences)
// {
//     acmacs::Counter<size_t> aligned_lengths;
//     for (const auto& seq : sequences)
//         aligned_lengths.count(seq.get().aa_aligned_length());
//     size_t master_length = 0;
//     for (const auto& vt : aligned_lengths.counter()) {
//         // fmt::print("  {:3d} {:5d} {:5d}\n", vt.first, vt.second, sequences.size() / 6);
//         if (vt.second > (sequences.size() / 6))
//             master_length = vt.first;
//     }
//     // fmt::print("aligned_lengths {}: {} ::: {}\n", sequences.front().get().type_subtype().h_or_b(), master_length, aligned_lengths);
//     return master_length;

// } // detect_master_length

// ----------------------------------------------------------------------

// class adjust_pos
// {
//  public:
//     static inline adjust_pos begin(std::string_view to_align, std::string_view master, size_t pos = 0) { return {to_align, master, pos}; }
//     static inline adjust_pos end(std::string_view to_align, std::string_view master) { return {to_align, master}; }

//     bool operator==(const adjust_pos& an) const { return /* mToAlign == an.mToAlign && mMaster == an.mMaster && */ mPos == an.mPos; }
//     bool operator!=(const adjust_pos& an) const { return ! operator==(an); }

//     size_t operator*() const { return mPos; }
//       // const size_t *operator->() const { return &mPos; }
//     adjust_pos& operator++() { ++mPos; find(); return *this; }

//  private:
//    adjust_pos(std::string_view to_align, std::string_view master, size_t pos) : mToAlign(to_align), mMaster(master), mLastPos(std::min(to_align.size(), master.size())), mPos(pos) { find(); }
//    adjust_pos(std::string_view to_align, std::string_view master) : mToAlign(to_align), mMaster(master), mLastPos(std::min(to_align.size(), master.size())), mPos(mLastPos) {}
//    std::string_view mToAlign, mMaster;
//    size_t mLastPos, mPos;

//    void find()
//    {
//        while (mPos < mLastPos && (common(mToAlign[mPos], mMaster[mPos]) || mToAlign[mPos] == '-' || mMaster[mPos] == '-'))
//            ++mPos;
//     }

// }; // class adjust_pos

// // ----------------------------------------------------------------------

// struct DeletionPos
// {
//     DeletionPos(size_t aPos, size_t aNumDeletions, size_t aNumCommon) : pos(aPos), num_deletions(aNumDeletions), num_common(aNumCommon) {}
//     bool operator<(const DeletionPos& aNother) const { return num_common == aNother.num_common ? pos < aNother.pos : num_common > aNother.num_common; }
//     size_t pos, num_deletions, num_common;

//     void fix(size_t aPos, size_t aNumDeletions, size_t aNumCommon) { pos = aPos; num_deletions = aNumDeletions; num_common = aNumCommon; }

//     friend inline std::ostream& operator<<(std::ostream& out, const DeletionPos& aPos) { return out << "pos:" << aPos.pos << " num_deletions:" << aPos.num_deletions << " num_common:" << aPos.num_common; }
// };

// using DeletionPosSet = std::vector<DeletionPos>;

// static inline void update(DeletionPosSet& pos_set, std::string_view master, std::string_view to_align, size_t pos, size_t common_before)
// {
//     constexpr const size_t max_num_deletions = 15;
//     const size_t last_pos = std::min(master.size(), to_align.size());
//     if ((pos + max_num_deletions) < last_pos) {
//         for (size_t num_insert = 1; num_insert <= max_num_deletions; ++num_insert) {
//             if (common(master[pos + num_insert], to_align[pos])) {
//                 pos_set.emplace_back(pos, num_insert, common_before + number_of_common(master.substr(pos + num_insert), to_align.substr(pos)));
//             }
//         }
//     }
// }

// // ----------------------------------------------------------------------

// void detect_deletions(acmacs::seqdb::sequence_t& sequence, const std::string_view master)
// {
//     constexpr const bool yamagata_163_hack = true;
//     constexpr const bool victoria_tripledel2017_hack = true;
//     const bool b_type = sequence.type_subtype() == acmacs::virus::type_subtype_t{"B"};

//     auto to_align = sequence.aa_aligned();
//     size_t best_common = number_of_common(master, to_align);
//     std::vector<std::pair<size_t, size_t>> pos_number;
//     size_t start = 0;

//     while (start < to_align.size()) {
//         const size_t current_common = number_of_common(master, to_align);
//         DeletionPosSet pos_set;
//         adjust_pos pos = adjust_pos::begin(to_align, master, start), pos_end = adjust_pos::end(to_align, master);
//         start = to_align.size();
//         for (; pos != pos_end; ++pos)
//             update(pos_set, master, to_align, *pos, number_of_common_before(master, to_align, *pos));

//         if (!pos_set.empty()) {
//             auto& del_pos = *std::min_element(pos_set.begin(), pos_set.end());
//             // dbg << "del_pos: " << del_pos << '\n';
//             if (del_pos.num_common > current_common) {
//                 if (yamagata_163_hack && b_type && del_pos.num_deletions == 1 && del_pos.pos > (163 - 1) && del_pos.pos <= (166 - 1)) {
//                       // std::cout << "INFO: yamagata_163_hack applied for " << entry_seq.make_name() << '\n';
//                     // yamagata deletion must be at 163
//                     // David Burke 2017-08-17: deletions ( and insertions) of amino acids usually occur in regions of the protein structure where it changes direction ( loops ).
//                     // In the case of HA, this is after VPK and before NKTAT/YKNAT.
//                     to_align.insert(163 - 1, 1, '-');
//                     del_pos.fix(163 - 1, 1, number_of_common(master, to_align)); // -1 because we count from zero here
//                 }
//                 else if (victoria_tripledel2017_hack && b_type && del_pos.num_deletions == 3 && del_pos.pos == (164 - 1)) {
//                     // The triple deletion is 162, 163 and 164 (pos 1 based). this is the convention that has been chosen (Sarah 2018-08-16 08:31)
//                     to_align.insert(162 - 1, del_pos.num_deletions, '-');
//                     del_pos.fix(162 - 1, del_pos.num_deletions, number_of_common(master, to_align)); // -1 because we count from zero here
//                 }
//                 else {
//                     to_align.insert(del_pos.pos, del_pos.num_deletions, '-');
//                 }
//                 start = del_pos.pos + del_pos.num_deletions + 1;
//                 pos_number.emplace_back(del_pos.pos, del_pos.num_deletions);
//                 if (best_common < del_pos.num_common)
//                     best_common = del_pos.num_common;
//             }
//         }
//     }
//     if (best_common < static_cast<size_t>(to_align.size() * 0.7)) {
//         fmt::print(stderr, "insertions? master:{} seq:{} {}\n{}\n{}\n{}\n\n", master.size(), sequence.aa_aligned_length(), *sequence.name(), master, to_align, sequence.aa_aligned());
//         // throw std::runtime_error("detect_deletions switch master");
//         // // dbg << "Too bad matching (common:" << best_common << " threshold:" << (master.size() * 0.7) << "), should we switch master?" << '\n';
//         // to_align = to_align_orig;
//         // throw SwitchMaster{};
//     }

//     if (!pos_number.empty() && (pos_number.size() > 1 || pos_number[0].second > 1))
//         fmt::print(stderr, "INSERTIONS {} master:{} seq:{} {}\n{}\n{}\n{}\n\n", *sequence.type_subtype(), master.size(), sequence.aa_aligned_length(), *sequence.name(), master, to_align, sequence.aa_aligned());
//     // return pos_number;

// } // detect

// // ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
