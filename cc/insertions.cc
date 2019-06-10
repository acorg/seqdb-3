#include <numeric>

#include "acmacs-base/counter.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/range-v3.hh"
#include "seqdb-3/insertions.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subtype_master_t acmacs::seqdb::v3::masters_per_subtype(const std::vector<fasta::scan_result_t>& sequences)
{
    std::map<std::string, acmacs::Counter<size_t>> aligned_lengths;
    for (const auto& sc : sequences | ranges::view::filter(fasta::is_aligned))
        aligned_lengths.try_emplace(std::string(sc.sequence.type_subtype().h_or_b())).first->second.count(sc.sequence.aa_aligned_length());

    subtype_master_t masters;
    for (const auto& [subtype, counter] : aligned_lengths) {
        const size_t threshold = counter.total() / 6;
        size_t master_length = 0;
        for (const auto& vt : counter.counter()) {
            if (vt.second > threshold)
                master_length = vt.first;
        }
        for (const auto& sc : sequences | ranges::view::filter(fasta::is_aligned)) {
            if (sc.sequence.type_subtype().h_or_b() == subtype && sc.sequence.aa_aligned_length() == master_length) {
                masters.emplace(subtype, &sc.sequence);
                break;
            }
        }
    }

    return masters;

} // acmacs::seqdb::v3::masters_per_subtype

// ----------------------------------------------------------------------

namespace local
{
    struct deletions_t
    {
        struct pos_num_t
        {
            size_t pos;
            size_t num;
        };

        std::vector<pos_num_t> data;

        void add(size_t pos, size_t num) { data.push_back({pos, num}); }
        bool empty() const { return data.empty(); }
        size_t number_of_deletions() const { return std::accumulate(std::begin(data), std::end(data), 0UL, [](size_t acc, const auto& en) { return acc + en.num; }); }

        std::string format(std::string_view sequence) const
        {
            fmt::memory_buffer out;
            size_t pos = 0;
            for (const auto& en : data) {
                fmt::format_to(out, "{}{}", sequence.substr(pos, en.pos - pos), std::string(en.num, '-'));
                pos = en.pos;
            }
            fmt::format_to(out, "{}", sequence.substr(pos));
            return fmt::to_string(out);
        }
    };

    // returns if deletion was found (and added)
    bool find_deletion(std::string_view to_align, std::string_view master, deletions_t& deletions);

    inline bool find_deletion(acmacs::seqdb::v3::sequence_t& to_align, const acmacs::seqdb::v3::sequence_t& master, deletions_t& deletions)
    {
        if (const auto [aligned, shift] = to_align.aa_shifted(); shift == 0)
            return find_deletion(aligned, master.aa_aligned_fast(), deletions);
        else
            return find_deletion(to_align.aa_aligned(), master.aa_aligned_fast(), deletions);
    }

} // namespace local

template <> struct fmt::formatter<local::deletions_t>
{
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }
    template <typename FormatContext> auto format(const local::deletions_t& deletions, FormatContext& ctx)
    {
        auto out = format_to(ctx.out(), "deletions({})[", deletions.data.size());
        for (const auto& en : deletions.data)
            out = format_to(out, " {}:{}", en.pos, en.num);
        return format_to(out, "]");
    }
};

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::insertions_deletions(sequence_t& to_align, const sequence_t& master)
{
    local::deletions_t deletions;
    if (to_align.aa_aligned_length() <= master.aa_aligned_length()) {
        find_deletion(to_align, master, deletions);
        if (!deletions.empty()) {
            fmt::print(stderr, "DEL {}\n{}\n{}\n\n", deletions, master.aa_aligned_fast(), deletions.format(to_align.aa_aligned()));
        }
        else if (to_align.aa_aligned_length() != master.aa_aligned_length()) {
            fmt::print(stderr, "!=\n{}\n{}\n", master.aa_aligned_fast(), to_align.aa_aligned());
        }
    }
    else {
        // fmt::print(stderr, "insertions_deletions {} > {} ::: {} {}\n{}\n{}\n", to_align.aa_aligned_length(), master.aa_aligned_length(), to_align.type_subtype(), to_align.full_name(),
        // master.aa_aligned_fast(), to_align.aa_aligned());
    }

} // acmacs::seqdb::v3::insertions_deletions

// ----------------------------------------------------------------------

namespace local
{
    static constexpr const auto are_common = [](char a, char b) -> bool { return a == b && a != 'X' && a != '-'; };
    // static constexpr const auto arenot_common = [](char a, char b) -> bool { return !are_common(a, b); };

    // static inline size_t number_of_common(std::string_view s1, std::string_view s2)
    // {
    //     size_t num = 0;
    //     for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
    //         if (are_common(*f1, *f2))
    //             ++num;
    //     }
    //     return num;
    // }

    template <typename Iter, typename Common> size_t find_head_tail_1(const Iter first1, const Iter last1, Iter first2, const Iter last2, const ssize_t threshold, Common common_f)
    {
        auto f1 = first1, last_common = first1 - 1;
        for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
            // fmt::print(stderr, " {}{}", *f1, *first2);
            if (common_f(*f1, *first2))
                last_common = f1;
            else if ((f1 - last_common) >= threshold)
                break;
        }
        // fmt::print(stderr, " {}\n", last_common - first1 + 1);
        return static_cast<size_t>(last_common - first1 + 1);
    }

    template <typename Iter, typename Common> size_t find_head_tail_2(const Iter first1, const Iter last1, Iter first2, const Iter last2, const ssize_t threshold, Common common_f)
    {
        auto f1 = first1, common_first = first1;
        for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
            // fmt::print(stderr, " {}{}", *f1, *first2);
            if (common_f(*f1, *first2)) {
                if ((f1 - common_first) >= (threshold - 1))
                    break;
            }
            else {
                common_first = f1 + 1;
            }
        }
        // fmt::print(stderr, " {}\n", common_first - first1);
        return static_cast<size_t>(common_first - first1);
    }

    // static inline size_t find_common_head(std::string_view s1, std::string_view s2, ssize_t threshold)
    // {
    //     return find_head_tail(s1.begin(), s1.end(), s2.begin(), s2.end(), threshold, are_common);
    // }

    // static inline size_t find_uncommon_head(std::string_view s1, std::string_view s2)
    // {
    //     return find_head_tail(s1.begin(), s1.end(), s2.begin(), s2.end(), 0, arenot_common);
    // }

    static inline size_t find_common_tail(std::string_view s1, std::string_view s2, ssize_t threshold)
    {
        return find_head_tail_1(s1.rbegin(), s1.rend(), s2.rbegin(), s2.rend(), threshold, are_common);
    }

    static inline size_t find_uncommon_tail(std::string_view s1, std::string_view s2, ssize_t threshold)
    {
        return find_head_tail_2(s1.rbegin(), s1.rend(), s2.rbegin(), s2.rend(), threshold, are_common);
    }

    // static inline size_t number_of_common_before(std::string_view s1, std::string_view s2, size_t last) { return number_of_common(s1.substr(0, last), s2.substr(0, last)); }

    bool find_deletion(std::string_view to_align, std::string_view master, deletions_t& deletions)
    {
        constexpr const ssize_t head_tail_threshold = 3;

        const auto tail = find_common_tail(master, to_align, head_tail_threshold);
        if (tail == to_align.size()) {
            if (master.size() == to_align.size())
                return false;
            else {
                // deletion at the beginning
                deletions.add(0, master.size() - to_align.size());
                return true;
            }
        }

        const auto front_aligned_head_size = to_align.size() - tail;
        // fmt::print(stderr, "{}\n{}\n   tail:{} front_aligned_head_size:{}\n", master, to_align, tail, front_aligned_head_size);
        // fmt::print(stderr, "{} {}\n{} {}\n:\n", master.substr(0, master.size() - tail), master.substr(master.size() - tail), to_align.substr(0, to_align.size() - tail),
        // to_align.substr(to_align.size() - tail));

        const auto middle = find_uncommon_tail(master.substr(0, front_aligned_head_size), to_align.substr(0, front_aligned_head_size), head_tail_threshold);
        const auto del_pos = front_aligned_head_size - middle;
        const auto del_num = 1; // master.size() - tail - front_aligned_head_size + middle;
        deletions.add(del_pos, del_num);
        if (del_num > 1) {
            fmt::print(stderr, "tail:{} middle:{} front_aligned_head_size:{} head:{} dels:{}\n", tail, middle, front_aligned_head_size, to_align.size() - tail - middle, del_num);
            fmt::print(stderr, "{} {} {}\n{} {} {}\n", master.substr(0, del_pos), master.substr(del_pos, del_num), master.substr(del_pos + del_num), to_align.substr(0, del_pos),
                       std::string(del_num, ' '), to_align.substr(del_pos));
        }
        return true;
    }

} // namespace local

// ----------------------------------------------------------------------

// 2019-06-09 ****************************************************************************************************

// local::deletions_t local::find_deletions(std::string_view to_align, std::string_view master)
// {
//     constexpr const ssize_t head_tail_threshold = 3;
//     auto head = find_common_head(master, to_align, head_tail_threshold);
//     if (head == to_align.size())
//         return deletions_t{}; // to_align truncated?

//     deletions_t result;
//     auto tail = find_common_tail(master.substr(head), to_align.substr(head), head_tail_threshold);
//     if (to_align.size() != master.size()) {

//         struct parts_t
//         {
//             std::vector<size_t> uncommon, common;
//         };
//         auto start = head;
//         auto master_chunk_size = master.size() - head - tail, to_align_chunk_size = to_align.size() - head - tail;
//         parts_t parts;
//         for (auto find_common = false; master_chunk_size > 0 && to_align_chunk_size > 0; find_common = !find_common) {
//             size_t part = 0;
//             if (find_common) {
//                 part = find_common_head(master.substr(start, master_chunk_size), to_align.substr(start, to_align_chunk_size), 0);
//                 parts.common.push_back(part);
//             }
//             else {
//                 part = find_uncommon_head(master.substr(start, master_chunk_size), to_align.substr(start, to_align_chunk_size));
//                 parts.uncommon.push_back(part);
//             }
//             if (part == 0)
//                 throw std::runtime_error("internal local::find_deletions");
//             start += part;
//             master_chunk_size -= part;
//             to_align_chunk_size -= part;
//         }

//         const auto show = [&](auto&& prefix) {
//             fmt::memory_buffer master_middle, to_align_middle;
//             size_t offset = 0;
//             for (size_t ind = 0; ind < std::max(parts.common.size(), parts.uncommon.size()); ++ind) {
//                 if (ind < parts.uncommon.size()) {
//                     fmt::format_to(master_middle, ".{}", master.substr(head + offset, parts.uncommon[ind]));
//                     fmt::format_to(to_align_middle, ".{}", to_align.substr(head + offset, parts.uncommon[ind]));
//                     offset += parts.uncommon[ind];
//                 }
//                 if (ind < parts.common.size()) {
//                     fmt::format_to(master_middle, "_{}", master.substr(head + offset, parts.common[ind]));
//                     fmt::format_to(to_align_middle, "_{}", to_align.substr(head + offset, parts.common[ind]));
//                     offset += parts.common[ind];
//                 }
//             }
//             fmt::print(stderr, "{}:\nuncommon: {:3d} {}\n  common: {:3d} {}\nhead:{} tail:{} to_align.size:{} master.size:{}\n{} {} {}\n{} {} {}\n\n", prefix, parts.uncommon.size(), parts.uncommon,
//                        parts.common.size(), parts.common, head, tail, to_align.size(), master.size(),
//                        master.substr(0, head), fmt::to_string(master_middle) /* master.substr(head, master.size() - head - tail) */, master.substr(master.size() - tail),
//                        to_align.substr(0, head), fmt::to_string(to_align_middle) /* to_align.substr(head, to_align.size() - head - tail) */, to_align.substr(to_align.size() - tail));
//         };

//         if (parts.common.size() != parts.uncommon.size()) {
//             show("?parts-size-differ");
//         }
//         else if (parts.common.size() == 1 && parts.uncommon[0] <= 3 && parts.common[0] > 2) {
//             head += parts.uncommon[0] + parts.common[0];
//             if ((head + tail) == to_align.size()) {
//                 result.push_back({head, master.size() - head - tail});
//             }
//             else {
//                 show("?head+tail");
//             }
//         }
//         else {
//             show("");
//         }
//     }
//     else {
//         const auto num_common = number_of_common(master.substr(head, master.size() - head - tail), to_align.substr(head, to_align.size() - head - tail));
//         fmt::print(stderr, "equal size common {} of {}\n", num_common, to_align.size() - head - tail);
//     }
//     return result;

// } // local::find_deletions

// ****************************************************************************************************

// void acmacs::seqdb::v3::insertions_deletions(std::vector<std::reference_wrapper<seqdb::sequence_t>>& sequences)
// {
//     // const auto& master = local::get_master(sequences);
//     // const auto master_seq = master.aa_aligned_fast();
//     // for (auto& seq : sequences) {
//     //     if (&seq.get() != &master)
//     //         detect_deletions(seq.get(), master_seq);
//     // }

// } // acmacs::seqdb::v3::insertions_deletions

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
