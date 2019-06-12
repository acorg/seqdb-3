#include <numeric>
#include <optional>

#include "acmacs-base/counter.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/date.hh"
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
        const sequence_t* master = nullptr;
        size_t num_X = 0;
        for (const auto& sc : sequences | ranges::view::filter(fasta::is_aligned)) {
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
        if (master == nullptr)
            throw std::runtime_error("internal in acmacs::seqdb::v3::masters_per_subtype");
        masters.emplace(subtype, master);
    }

    return masters;

} // acmacs::seqdb::v3::masters_per_subtype

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::deletions_insertions(const sequence_t& master, sequence_t& to_align)
{
    deletions_insertions_t deletions;
    const auto [master_aligned, master_shift] = master.aa_shifted();
    const auto [to_align_aligned, to_align_shift] = to_align.aa_shifted();
    if (master_shift == 0) {
        if (to_align_shift == 0)
            deletions = deletions_insertions(master_aligned, to_align_aligned);
        else
            deletions = deletions_insertions(master_aligned, to_align.aa_aligned());
    }
    else {
        if (to_align_shift == 0)
            deletions = deletions_insertions(master.aa_aligned(), to_align_aligned);
        else
            deletions = deletions_insertions(master.aa_aligned(), to_align.aa_aligned());
    }

} // acmacs::seqdb::v3::deletions_insertions

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::format(const std::vector<deletions_insertions_t::pos_num_t>& pos_num, std::string_view sequence, char deletion_symbol)
{
    fmt::memory_buffer out;
    size_t pos = 0;
    for (const auto& en : pos_num) {
        fmt::format_to(out, "{}{}", sequence.substr(pos, en.pos - pos), std::string(en.num, deletion_symbol));
        pos = en.pos;
    }
    fmt::format_to(out, "{}", sequence.substr(pos));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::format

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::format(const deletions_insertions_t& deletions)
{
    fmt::memory_buffer out;
    const auto frmt = [&out](const char* prefix, const auto& num_pos) {
        if (!num_pos.empty()) {
            fmt::format_to(out, "{}[{}](", prefix, num_pos.size());
            for (const auto& en : num_pos)
                fmt::format_to(out, " {}:{}", en.pos, en.num);
            fmt::format_to(out, ")");
        }
    };

    frmt("DEL", deletions.deletions);
    frmt(" INS", deletions.insertions);
    return fmt::to_string(out);

} // acmacs::seqdb::v3::format

// ----------------------------------------------------------------------

namespace local
{
    constexpr const ssize_t common_threshold = 3; // assume the chunk is common after that number of consecutive common positions
    constexpr const ssize_t different_threshold = 20; // assume rest is different after that number of consecutive different positions

    static constexpr const auto are_common = [](char a, char b) -> bool { return a == b && a != 'X' && a != '-'; };

    struct find_head_t
    {
        size_t head = 0;
        size_t common = 0;
    };

    template <typename Iter, typename Common> find_head_t find_head(const Iter first1, const Iter last1, Iter first2, const Iter last2, Common common_f)
    {
        // find the last part with common_f()==true that is not shorter than common_threshold
        // returns offset of the end of this part
        auto f1 = first1, common_start = last1, last_common_end = first1;
        size_t common = 0, common_at_last_common_end = 0;
        const auto update_last_common_end = [&]() {
            if (common_start != last1 && (f1 - common_start) >= common_threshold) {
                last_common_end = f1;
                common_at_last_common_end = common;
                common_start = last1; // reset
            }
        };

        for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
            if (!common_f(*f1, *first2)) {
                update_last_common_end();
                if (static_cast<size_t>(f1 - last_common_end) > different_threshold)
                    break; // too many different, stop searching
            }
            else {
                ++common;
                if (common_start == last1)
                    common_start = f1;
            }
        }
        update_last_common_end();

        if (const auto head = static_cast<size_t>(last_common_end - first1); common_at_last_common_end * 3 > head)
            return {head, common_at_last_common_end};
        else
            return {0, 0};      // too few common in the head, try more deletions
    }

    static inline find_head_t find_common_head(std::string_view s1, std::string_view s2)
    {
        return find_head(s1.begin(), s1.end(), s2.begin(), s2.end(), are_common);
    }

    struct deletions_insertions_at_start_t
    {
        size_t deletions = 0;
        size_t insertions = 0;
        find_head_t head;
    };

    inline deletions_insertions_at_start_t deletions_insertions_at_start(std::string_view master, std::string_view to_align)
    {
        deletions_insertions_at_start_t result;
        for (size_t dels = 1; dels < different_threshold; ++dels) {
            if (dels < master.size()) {
                result.head = find_common_head(master.substr(dels), to_align);
                // fmt::print(stderr, "dels:{} head:{} common:{}\n{}\n{}\n", dels, result.head.head, result.head.common, master.substr(dels), to_align);
                if (result.head.head > common_threshold) {
                    result.deletions = dels;
                    break;
                }
            }
            if (dels < to_align.size()) {
                result.head = find_common_head(master, to_align.substr(dels));
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

    inline size_t number_of_common(std::string_view master, std::string_view to_align, acmacs::seqdb::v3::deletions_insertions_t& deletions)
    {
        return number_of_common(acmacs::seqdb::v3::format(deletions.insertions, master), acmacs::seqdb::v3::format(deletions.deletions, to_align));
    }

}

// ----------------------------------------------------------------------

acmacs::seqdb::v3::deletions_insertions_t acmacs::seqdb::v3::deletions_insertions(std::string_view master, std::string_view to_align, debug dbg)
{
    if (dbg == debug::yes)
        fmt::print(stderr, "initial:\n{}\n{}\n\n", master, to_align);

    deletions_insertions_t deletions;
    const auto initial_head = local::find_common_head(master, to_align);
    size_t master_offset = initial_head.head, to_align_offset = initial_head.head;
    std::string_view master_tail = master.substr(master_offset), to_align_tail = to_align.substr(to_align_offset);
    if (dbg == debug::yes)
        fmt::print(stderr, "initial_head:{} common:{} number_of_common:{}\n", initial_head.head, initial_head.common, local::number_of_common(master.substr(0, initial_head.head), to_align.substr(0, initial_head.head)));

    const auto update_both = [](size_t dels, size_t head, std::string_view& s1, std::string_view& s2, size_t& offset1, size_t& offset2) {
        s1.remove_prefix(head + dels);
        s2.remove_prefix(head);
        offset1 += head + dels;
        offset2 += head;
    };

    size_t common = initial_head.common;
    while (!master_tail.empty() && !to_align_tail.empty()) {
        if (dbg == debug::yes)
            fmt::print(stderr, "m-offset:{} a-offset:{} common:{}\n{}\n{}\n", master_offset, to_align_offset, common, master_tail, to_align_tail);
        const auto tail_deletions = local::deletions_insertions_at_start(master_tail, to_align_tail);
        if (dbg == debug::yes)
            fmt::print(stderr, "dels:{} ins:{} head:{}  common:{} number_of_common:{}\n", tail_deletions.deletions, tail_deletions.insertions, tail_deletions.head.head, tail_deletions.head.common, local::number_of_common(master_tail.substr(tail_deletions.deletions, tail_deletions.head.head), to_align_tail.substr(tail_deletions.insertions, tail_deletions.head.head)));
        if (tail_deletions.head.head == 0) { // < local::common_threshold) {
            common += local::number_of_common(master_tail, to_align_tail); // to avoid common diff warning below in case tail contains common aas
            break; // tails are different, insertions/deletions do not help
        }

        if (tail_deletions.deletions) {
            deletions.deletions.push_back({to_align_offset, tail_deletions.deletions});
            update_both(tail_deletions.deletions, tail_deletions.head.head, master_tail, to_align_tail, master_offset, to_align_offset);
        }
        else if (tail_deletions.insertions) {
            deletions.insertions.push_back({master_offset, tail_deletions.insertions});
            update_both(tail_deletions.insertions, tail_deletions.head.head, to_align_tail, master_tail, to_align_offset, master_offset);
        }

        common += tail_deletions.head.common;
    }

    // sanity check (remove)
    if (const auto nc = local::number_of_common(master, to_align, deletions); nc != common)
        fmt::print(stderr, "common diff: {} vs. number_of_common:{} {}\n{}\n{}\n", common, nc, format(deletions), acmacs::seqdb::v3::format(deletions.insertions, master, '.'), acmacs::seqdb::v3::format(deletions.deletions, to_align, '.'));

    // verify
    constexpr double diff_threshold = 0.7;
    const auto num_common_threshold = (to_align.size() - static_cast<size_t>(std::count(std::begin(to_align), std::end(to_align), 'X'))) * diff_threshold;
    if (common < num_common_threshold) {
        fmt::print(stderr, "------ NOT VERIFIED common:{} vs size:{} num_common_threshold:{:.2f} ----------\n{}\n{}\n{}\n{}\n", common, to_align.size(), num_common_threshold, master, to_align,
                   acmacs::seqdb::v3::format(deletions.insertions, master, '.'), acmacs::seqdb::v3::format(deletions.deletions, to_align, '.'));
    }

    return deletions;

} // acmacs::seqdb::v3::deletions_insertions

// ----------------------------------------------------------------------

// 2019-06-11 ****************************************************************************************************

// namespace local
// {
//     static constexpr const auto are_common = [](char a, char b) -> bool { return a == b && a != 'X' && a != '-'; };

//     struct deletions_t
//     {
//         struct pos_num_t
//         {
//             size_t pos;
//             size_t num;
//         };

//         std::vector<pos_num_t> data;

//         void add(size_t pos, size_t num) { data.push_back({pos, num}); }
//         bool empty() const { return data.empty(); }
//         size_t size() const { return data.size(); }
//         const auto& back() const { return data.back(); }
//         void increment_last_num() { ++data.back().num; }
//         void pop_back() { data.pop_back(); }
//         size_t number_of_deletions() const { return std::accumulate(std::begin(data), std::end(data), 0UL, [](size_t acc, const auto& en) { return acc + en.num; }); }
//         size_t number_of_deletions_without_last() const
//         {
//             if (data.size() < 2)
//                 return 0;
//             return std::accumulate(std::begin(data), std::end(data) - 1, 0UL, [](size_t acc, const auto& en) { return acc + en.num; });
//         }

//         std::string format(std::string_view sequence) const
//         {
//             fmt::memory_buffer out;
//             size_t pos = 0;
//             for (const auto& en : data) {
//                 fmt::format_to(out, "{}{}", sequence.substr(pos, en.pos - pos), std::string(en.num, '-'));
//                 pos = en.pos;
//             }
//             fmt::format_to(out, "{}", sequence.substr(pos));
//             return fmt::to_string(out);
//         }
//     };

//     struct find_deletions_t
//     {
//         size_t num_dels = 0;
//         size_t num_common = 0;
//     };

//     enum class debug { no, yes };

//     find_deletions_t find_deletions(std::string_view master, std::string_view to_align, deletions_t& deletions, debug dbg);

//     inline find_deletions_t find_deletions(const acmacs::seqdb::v3::sequence_t& master, const acmacs::seqdb::v3::sequence_t& to_align, deletions_t& deletions, debug dbg)
//     {
//         const auto [master_aligned, master_shift] = master.aa_shifted();
//         const auto [to_align_aligned, to_align_shift] = to_align.aa_shifted();
//         if (master_shift == 0) {
//             if (to_align_shift == 0)
//                 return find_deletions(master_aligned, to_align_aligned, deletions, dbg);
//             else
//                 return find_deletions(master_aligned, to_align.aa_aligned(), deletions, dbg);
//         }
//         else {
//             if (to_align_shift == 0)
//                 return find_deletions(master.aa_aligned(), to_align_aligned, deletions, dbg);
//             else
//                 return find_deletions(master.aa_aligned(), to_align.aa_aligned(), deletions, dbg);
//         }
//     }

//     inline size_t number_of_common(std::string_view s1, std::string_view s2)
//     {
//         size_t num = 0;
//         for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
//             if (are_common(*f1, *f2))
//                 ++num;
//         }
//         return num;
//     }

//     inline size_t number_of_common(std::string_view master, std::string_view to_align, deletions_t& deletions)
//     {
//         size_t num = 0, offset = 0, num_dels = 0;
//         for (const auto& del : deletions.data) {
//             num += number_of_common(master.substr(offset + num_dels, del.pos - offset), to_align.substr(offset, del.pos - offset));
//             num_dels += del.num;
//             offset = del.pos;
//         }
//         num += number_of_common(master.substr(offset + num_dels, master.size() - offset), to_align.substr(offset, to_align.size() - offset));
//         return num;
//     }

//     inline size_t number_of_common(const acmacs::seqdb::v3::sequence_t& master, const acmacs::seqdb::v3::sequence_t& to_align, deletions_t& deletions)
//     {
//         const auto [master_aligned, master_shift] = master.aa_shifted();
//         const auto [to_align_aligned, to_align_shift] = to_align.aa_shifted();
//         if (master_shift == 0) {
//             if (to_align_shift == 0)
//                 return number_of_common(master_aligned, to_align_aligned, deletions);
//             else
//                 return number_of_common(master_aligned, to_align.aa_aligned(), deletions);
//         }
//         else {
//             if (to_align_shift == 0)
//                 return number_of_common(master.aa_aligned(), to_align_aligned, deletions);
//             else
//                 return number_of_common(master.aa_aligned(), to_align.aa_aligned(), deletions);
//         }
//     }

// } // namespace local

// template <> struct fmt::formatter<local::deletions_t>
// {
//     template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }
//     template <typename FormatContext> auto format(const local::deletions_t& deletions, FormatContext& ctx)
//     {
//         auto out = format_to(ctx.out(), "deletions({})[", deletions.data.size());
//         for (const auto& en : deletions.data)
//             out = format_to(out, " {}:{}", en.pos, en.num);
//         return format_to(out, "]");
//     }
// };

// // ----------------------------------------------------------------------

// void acmacs::seqdb::v3::insertions_deletions(const sequence_t& master, sequence_t& to_align)
// {
//     constexpr double diff_threshold = 0.7;
//     const auto num_common_threshold = to_align.aa_number_of_not_X() * diff_threshold;

//     local::deletions_t deletions;
//     const auto fdr = local::find_deletions(master, to_align, deletions, local::debug::no);
//     if (fdr.num_common < num_common_threshold) {
//         local::deletions_t insertions;
//         const auto fir = local::find_deletions(to_align, master, insertions, local::debug::no);
//         if (fir.num_common < num_common_threshold) {
//             fmt::print(stderr, "------ NOT VERIFIED {} vs {} ({:.2f}) dels:{} ins:{} ----------\n{}\n{}\n{}\n{}\n",
//                        fdr.num_common, to_align.aa_aligned_length(), double(fdr.num_common) / to_align.aa_aligned_length(), fdr.num_dels, fir.num_dels,
//                        master.full_name(), to_align.full_name(),
//                        master.aa_aligned_fast(), to_align.aa_aligned());
//             local::deletions_t deletions2, insertions2;
//             fmt::print(stderr, "deletions\n");
//             local::find_deletions(master, to_align, deletions2, local::debug::yes);
//             fmt::print(stderr, "insertions\n");
//             local::find_deletions(to_align, master, insertions2, local::debug::yes);
//             fmt::print(stderr, "\n");
//         }
//         else {
//             // good, set insertions
//         }
//     }
//     else {
//         // good, set deletions
//     }

//     // if (to_align.aa_aligned_length() <= master.aa_aligned_length()) {
//     // }
//     // else {
//     //     // fmt::print(stderr, "insertions_deletions {} > {} ::: {} {}\n{}\n{}\n", to_align.aa_aligned_length(), master.aa_aligned_length(), to_align.type_subtype(), to_align.full_name(),
//     //     // master.aa_aligned_fast(), to_align.aa_aligned());
//     // }

// } // acmacs::seqdb::v3::insertions_deletions

// // ----------------------------------------------------------------------

// namespace local
// {
//     constexpr const ssize_t common_threshold = 3;

//     struct find_head_t
//     {
//         size_t head;
//         size_t num_common;
//     };

//     template <typename Iter, typename Common> find_head_t find_head(const Iter first1, const Iter last1, Iter first2, const Iter last2, Common common_f)
//     {
//         // find the last part with common_f()==true that is not shorter than common_threshold
//         // returns offset of the end of this part
//         auto f1 = first1, common_start = last1, last_common_end = first1;
//         size_t num_common = 0;
//         const auto update_last_common_end = [&]() {
//             if (common_start != last1 && (f1 - common_start) >= common_threshold)
//                 last_common_end = f1;
//         };

//         for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
//             if (!common_f(*f1, *first2)) {
//                 update_last_common_end();
//                 common_start = last1; // reset
//             }
//             else {
//                 ++num_common;
//                 if (common_start == last1)
//                     common_start = f1;
//             }
//         }
//         update_last_common_end();
//         const auto head = static_cast<size_t>(last_common_end - first1);
//         if (num_common * 3 > head)
//             return {head, num_common};
//         else
//             return {0, 0};      // too few common in the head, try more deletions
//     }

//     static inline find_head_t find_common_head(std::string_view s1, std::string_view s2)
//     {
//         return find_head(s1.begin(), s1.end(), s2.begin(), s2.end(), are_common);
//     }

//     find_deletions_t find_deletions(std::string_view master, std::string_view to_align, deletions_t& deletions, debug dbg)
//     {
//         size_t offset = 0;
//         find_deletions_t fd;
//         for (auto loop = 0; offset < to_align.size() && (offset + fd.num_dels) < master.size(); ++loop) {
//             try {
//                 const auto head_common = find_common_head(master.substr(offset + fd.num_dels), to_align.substr(offset));
//                 if (dbg == debug::yes)
//                     fmt::print(stderr, "head:{} common:{}\n{}\n{}\n", head_common.head, head_common.num_common, master.substr(offset + fd.num_dels), to_align.substr(offset));
//                 offset += head_common.head;
//                 if (head_common.head == 0 && !deletions.empty()) {
//                     deletions.increment_last_num();
//                     ++fd.num_dels;
//                 }
//                 else {
//                     if (head_common.head != 0)
//                         fd.num_common += head_common.num_common;
//                     if (offset < to_align.size()) {
//                         deletions.add(offset, 1);
//                         ++fd.num_dels;
//                     }
//                 }
//                 if (fd.num_dels && offset < to_align.size() && (offset + fd.num_dels) >= master.size()) {
//                     // master and to_align tails differ and deletions do not help
//                     deletions.pop_back();
//                 }
//             }
//             catch (...) {
//                 fmt::print(stderr, "------ EXCEPTION offset:{} num_dels:{} loop:{} ----------\n{}\n{}\n\n", offset, fd.num_dels, loop, master, deletions.format(to_align));
//                 throw;
//             }
//         }
//         return fd;
//     }

// } // namespace local

// 2019-06-11 ****************************************************************************************************

// namespace local
// {
//     static constexpr const auto are_common = [](char a, char b) -> bool { return a == b && a != 'X' && a != '-'; };
//     // static constexpr const auto arenot_common = [](char a, char b) -> bool { return !are_common(a, b); };

//     // static inline size_t number_of_common(std::string_view s1, std::string_view s2)
//     // {
//     //     size_t num = 0;
//     //     for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
//     //         if (are_common(*f1, *f2))
//     //             ++num;
//     //     }
//     //     return num;
//     // }

//     template <typename Iter, typename Common> size_t find_head_tail_1(const Iter first1, const Iter last1, Iter first2, const Iter last2, const ssize_t threshold, Common common_f)
//     {
//         auto f1 = first1, last_common = first1 - 1;
//         for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
//             // fmt::print(stderr, " {}{}", *f1, *first2);
//             if (common_f(*f1, *first2))
//                 last_common = f1;
//             else if ((f1 - last_common) >= threshold)
//                 break;
//         }
//         // fmt::print(stderr, " {}\n", last_common - first1 + 1);
//         return static_cast<size_t>(last_common - first1 + 1);
//     }

//     template <typename Iter, typename Common> size_t find_head_tail_2(const Iter first1, const Iter last1, Iter first2, const Iter last2, const ssize_t threshold, Common common_f)
//     {
//         auto f1 = first1, common_first = first1;
//         for (; f1 != last1 && first2 != last2; ++f1, ++first2) {
//             // fmt::print(stderr, " {}{}", *f1, *first2);
//             if (common_f(*f1, *first2)) {
//                 if ((f1 - common_first) >= (threshold - 1))
//                     break;
//             }
//             else {
//                 common_first = f1 + 1;
//             }
//         }
//         // fmt::print(stderr, " {}\n", common_first - first1);
//         return static_cast<size_t>(common_first - first1);
//     }

//     // static inline size_t find_common_head(std::string_view s1, std::string_view s2, ssize_t threshold)
//     // {
//     //     return find_head_tail_1(s1.begin(), s1.end(), s2.begin(), s2.end(), threshold, are_common);
//     // }

//     // static inline size_t find_uncommon_head(std::string_view s1, std::string_view s2)
//     // {
//     //     return find_head_tail(s1.begin(), s1.end(), s2.begin(), s2.end(), 0, arenot_common);
//     // }

//     static inline std::pair<size_t, size_t> find_common_tail(std::string_view s1, std::string_view s2, ssize_t threshold)
//     {
//         for (ssize_t s1_offset = 0; static_cast<size_t>(s1_offset) < s1.size(); ++s1_offset) {
//             if (const auto tail = find_head_tail_1(s1.rbegin() + s1_offset, s1.rend(), s2.rbegin(), s2.rend(), threshold, are_common); tail >= static_cast<size_t>(threshold))
//                 return {tail + static_cast<size_t>(s1_offset), tail};
//         }
//         throw std::runtime_error("insertions.cc local::find_head_tail internal");
//     }

//     static inline size_t find_uncommon_tail(std::string_view s1, std::string_view s2, ssize_t threshold)
//     {
//         return find_head_tail_2(s1.rbegin(), s1.rend(), s2.rbegin(), s2.rend(), threshold, are_common);
//     }

//     // static inline size_t number_of_common_before(std::string_view s1, std::string_view s2, size_t last) { return number_of_common(s1.substr(0, last), s2.substr(0, last)); }

//     void find_deletions(std::string_view to_align, std::string_view master, deletions_t& deletions)
//     {
//         // fmt::memory_buffer DEBUG;

//         constexpr const ssize_t head_tail_threshold = 3;

//         // const auto find_deletion = [](std::string_view to_align_3, std::string_view master_3) -> std::optional<size_t> {
//         //     const auto middle = find_uncommon_tail(master_3.substr(0, to_align_3.size()), to_align_3, head_tail_threshold);
//         //     if (middle == 0)
//         //         return std::nullopt;
//         //     return to_align_3.size() - middle;
//         // };

//         // fmt::format_to(DEBUG, "------------------------------\n{}\n{}\n", master, to_align);

//         try {
//             const auto [master_tail, to_align_tail] = find_common_tail(master, to_align, head_tail_threshold);
//             const auto master_extra_tail = master_tail - to_align_tail;
//             if (to_align_tail == to_align.size()) {
//                 if (master.size() > (to_align.size() + master_extra_tail))
//                     deletions.add(0, master.size() - master_extra_tail - to_align.size()); // deletion at the beginning
//                 return;
//             }
//             const auto to_align_without_tail = to_align.substr(0, to_align.size() - to_align_tail), master_without_tail = master.substr(0, master.size() - master_tail);

//             const auto middle = find_uncommon_tail(master_without_tail.substr(0, to_align_without_tail.size()), to_align_without_tail, head_tail_threshold - 1);
//             const auto head = to_align_without_tail.size() - middle;
//             if ((master.size() - master_extra_tail - to_align.size()) == 1) {
//                 deletions.add(head, 1);
//                 // if (head != 162) {
//                 //     fmt::print(stderr, "------\n{}\n{}\nhead:{} middle:{}\n{}\n{}\n", master, to_align, head, middle, master, deletions.format(to_align));
//                 // }
//                 return;
//             }

//             fmt::print(stderr, "------\n{}\n{}\nmiddle:{} m-tail:{} a-tail:{}\n{} {} {}\n{} {} {}\n", master, to_align, middle, master_tail, to_align_tail,
//                            master_without_tail.substr(0, head), master_without_tail.substr(head), master.substr(master.size() - master_tail),
//                        to_align_without_tail.substr(0, head), to_align_without_tail.substr(head, middle), to_align.substr(to_align.size() - to_align_tail)
//                        );
//             // if (middle == to_align_without_tail.size()) {
//             // }

//             // // ----------------------------------------------------------------------

//             // const auto head = find_common_head(master_without_tail, to_align_without_tail, head_tail_threshold);
//             // if (head == to_align_without_tail.size()) {
//             //     if (master_without_tail.size() > to_align_without_tail.size())
//             //         deletions.add(0, master_without_tail.size() - to_align_without_tail.size()); // deletion at the end of head
//             //     return;
//             // }

//             // fmt::format_to(DEBUG, "\n{} {}\n{} {}\n", master_without_tail.substr(0, head), master_without_tail.substr(head), to_align_without_tail.substr(0, head), to_align_without_tail.substr(head));

//             // const auto to_align_middle = to_align_without_tail.substr(head), master_middle = master_without_tail.substr(head);
//             // if (to_align_middle.size() < master_middle.size())
//             //     deletions.add(head, master_middle.size() - to_align_middle.size());

//             // if (to_align_middle.size() > 3)
//             //     fmt::print(stderr, "{}DEL {}\n{}\n{}\n\n", fmt::to_string(DEBUG), deletions, master, deletions.format(to_align));

//             // // remove deletions at the end
//             // if (!deletions.empty() && master.size() == (deletions.back().pos + deletions.back().num + deletions.number_of_deletions_without_last()))
//             //     deletions.pop_back();

//             // // ----------------------------------------------------------------------

//             //     size_t to_align_start = head;
//             // std::optional<size_t> del_pos;
//             // while ((del_pos = find_deletion(to_align_without_tail.substr(to_align_start), master_without_tail.substr(to_align_start + deletions.number_of_deletions())))
//             //            .has_value()) {
//             //     to_align_start += *del_pos;
//             //     if (!deletions.empty() && *del_pos == 0)
//             //         deletions.increment_last_num();
//             //     else
//             //         deletions.add(to_align_start, 1);
//             //     fmt::format_to(DEBUG, "{}\n{}\n", deletions, deletions.format(to_align));
//             //     if (const auto nd = deletions.number_of_deletions(); nd > 5 && (to_align_without_tail.size() + nd) > master_without_tail.size()) {
//             //         break;
//             //     }
//             // }

//             // if (master.size() != (to_align.size() + deletions.number_of_deletions())) {
//             //     fmt::format_to(DEBUG, "WARNING: != {}\n{}\n{}\n", deletions, master, deletions.format(to_align));
//             // }

//             // if (!deletions.empty() && (deletions.size() > 1 || deletions.back().pos != 162 || deletions.back().num > 1)) {
//             //     fmt::format_to(DEBUG, "DEL {}\n{}\n{}\n", deletions, master, deletions.format(to_align));
//             //     fmt::print(stderr, "{}\n", fmt::to_string(DEBUG));
//             // }
//         }
//         catch (std::exception& err) {
//             fmt::print(stderr, "------\nERROR: {}\n{}\n{}\n", err, master, to_align);
//             throw;
//         }
//     }

// } // namespace local

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
