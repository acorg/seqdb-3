#include "acmacs-base/named-type.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/fmt.hh"
#include "seqdb-3/align.hh"
#include "seqdb-3/hamming-distance.hh"

// http://signalpeptide.com

// ----------------------------------------------------------------------

namespace align_detail
{
    using max_offset_t = acmacs::named_size_t<struct max_offset_t_tag>;
    using hdth_t = acmacs::named_size_t<struct hdth_t_tag>; // hamming_distance_threshold
    using shift_t = acmacs::named_size_t<struct shift_t_tag>;
    constexpr const shift_t shift_is_pattern_size{-1};

    struct pat_t
    {
        std::string_view type_subtype_prefix;
        std::string_view type_subtype;
        std::string_view pattern;
        max_offset_t max_offset;
        hdth_t hamming_distance_threshold;
        shift_t shift;
    };

    inline std::string_view type_subtype_hint(std::string_view type_subtype)
    {
        return type_subtype.size() > 5 && type_subtype[0] == 'A' ? type_subtype.substr(0, 5) : type_subtype;
    }

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

    static const std::array patterns{
        pat_t{"A(H3N", "A(H3N2)", "MKTIIALSYILCLVFA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size},
        // pat_t{"A(H3N", "A(H3N2)", "MKTIIALSYIFCLALG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size}, // A/Aichi/2/1968
        // pat_t{"A(H3N", "A(H3N2)", "MKTLIALSYIFCLVLG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size},
        // pat_t{"A(H3N", "A(H3N2)", "MKTIIALSYIFCLALG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size}, // A/Hong Kong/1/1968
        // pat_t{"A(H3N", "A(H3N2)", "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENC", max_offset_t{100}, hdth_t{6}, shift_t{0}},
        // pat_t{"A(H3N", "A(H3N2)", "QKLPGNNNSTATLCLGHHAVPNGTIVKTI",                                    max_offset_t{100}, hdth_t{6}, shift_t{0}},
    };

#pragma GCC diagnostic pop

} // namespace align_detail

std::optional<std::tuple<int, std::string_view>> acmacs::seqdb::v3::align(std::string_view amino_acids, std::string_view type_subtype_hint, std::string_view /*debug_name*/)
{
    const auto hint = align_detail::type_subtype_hint(type_subtype_hint);
    const auto make_type_subtype = [hint, type_subtype_hint](const auto& pattern) -> std::string_view {
        if (hint == pattern.type_subtype_prefix)
            return type_subtype_hint;
        else
            return pattern.type_subtype;
    };
    const auto suitable_subtype = [hint](const auto& pattern) -> bool { return pattern.type_subtype_prefix == hint || hint.empty() || hint == "A(H0N"; };

    for (const auto& pattern : align_detail::patterns | ranges::view::filter(suitable_subtype)) {
        for (auto p_start = amino_acids.find(pattern.pattern[0]); p_start < *pattern.max_offset; p_start = amino_acids.find(pattern.pattern[0], p_start + 1)) {
            if (hamming_distance(pattern.pattern, amino_acids.substr(p_start, pattern.pattern.size())) < *pattern.hamming_distance_threshold) {
                if (pattern.shift == align_detail::shift_is_pattern_size)
                    return std::tuple{static_cast<int>(p_start + pattern.pattern.size()), make_type_subtype(pattern)};
                else
                    return std::tuple{static_cast<int>(p_start + *pattern.shift), make_type_subtype(pattern)};
            }
        }
    }
    return std::nullopt;

} // acmacs::seqdb::v3::align

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Aligner::update(std::string_view amino_acids, int shift, std::string_view type_subtype)
{
    const std::string hint{align_detail::type_subtype_hint(type_subtype)};
    tables_.try_emplace(hint).first->second.update(amino_acids, shift);

} // acmacs::seqdb::v3::Aligner::update

// ----------------------------------------------------------------------

std::optional<std::tuple<int, std::string_view>> acmacs::seqdb::v3::Aligner::align(std::string_view amino_acids, std::string_view type_subtype_hint, std::string_view debug_name) const
{
    const auto hint = align_detail::type_subtype_hint(type_subtype_hint);

} // acmacs::seqdb::v3::Aligner::align

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Aligner::report() const
{
    fmt::print(stderr, "Aligner {}\n", tables_.size());
    for (const auto& [type_subtype, table] : tables_)
        table.report(fmt::format(" {:<8s} ", type_subtype));
    fmt::print(stderr, "\n");

} // acmacs::seqdb::v3::Aligner::report

// ----------------------------------------------------------------------

acmacs::seqdb::v3::Aligner::table_t::table_t()
{
    ranges::fill(data, 1);

    // X and - do not contribute at any position
    for (auto pos : ranges::view::iota(0UL, max_sequence_length)) {
        data[number_of_symbols * pos + static_cast<size_t>('X')] = 0;
        data[number_of_symbols * pos + static_cast<size_t>('-')] = 0;
    }

} // acmacs::seqdb::v3::Aligner::table_t::table_t

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Aligner::table_t::update(std::string_view amino_acids, int shift)
{
    auto pos = static_cast<size_t>(-shift);
    for (char aa : amino_acids) {
        data[number_of_symbols * pos + static_cast<size_t>(aa)] = 0;
        ++pos;
    }

} // acmacs::seqdb::v3::Aligner::table_t::update

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Aligner::table_t::report(std::string prefix) const
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
    for (auto pos : ranges::view::iota(0UL, max_sequence_length)) {
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
        for (auto pos : ranges::view::iota(0UL, last_pos)) {
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

} // acmacs::seqdb::v3::Aligner::table_t::report

// ----------------------------------------------------------------------

//     // if (aa_.empty())
//     //     return false;
//     // if (type_subtype_hint == "A(H3N2)")
//     //     return align_h3n2(debug_name) || align_any(debug_name, type_subtype_hint);
//     // else
//     //     return align_any(debug_name, type_subtype_hint);

// } // acmacs::seqdb::v3::sequence_t::align

// // ----------------------------------------------------------------------


// #pragma GCC diagnostic push
// #ifdef __clang__
// #pragma GCC diagnostic ignored "-Wglobal-constructors"
// #pragma GCC diagnostic ignored "-Wexit-time-destructors"
// #endif

// using max_offset_t = acmacs::named_size_t<struct max_offset_t_tag>;
// using hdth_t = acmacs::named_size_t<struct hdth_t_tag>; // hamming_distance_threshold
// using shift_t = acmacs::named_size_t<struct shift_t_tag>;
// constexpr const shift_t shift_is_pattern_size{-1};

// struct pat_t
// {
//     std::string_view pattern;
//     max_offset_t max_offset;
//     hdth_t hamming_distance_threshold;
//     shift_t shift;
// };

// static const std::array sH3patterns{
//     //pat_t{"MKTIIALSYIFCLALG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size}, // A/Aichi/2/1968
//     pat_t{"MKTIIALSYILCLVFA",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size},
//     // pat_t{"MKTLIALSYIFCLVLG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size},
//     // pat_t{"MKTIIALSYIFCLALG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size}, // A/Hong Kong/1/1968
//     // pat_t{"QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENC", max_offset_t{100}, hdth_t{6}, shift_t{0}},
//     // pat_t{"QKLPGNNNSTATLCLGHHAVPNGTIVKTI",                                    max_offset_t{100}, hdth_t{6}, shift_t{0}},
// };

// #pragma GCC diagnostic pop

// bool acmacs::seqdb::v3::sequence_t::align_h3n2(std::string_view debug_name)
// {
//     const std::string_view data{aa_};
//     std::vector<size_t> hamds;

//     for (const auto& pattern : sH3patterns) {
//         const std::string_view look_for = pattern.pattern.substr(0, 2);
//         for (auto p1_start = data.find(look_for); p1_start < *pattern.max_offset; p1_start = data.find(look_for, p1_start + look_for.size())) {
//             if (const auto hamd = hamming_distance(pattern.pattern, data.substr(p1_start, pattern.pattern.size())); hamd < *pattern.hamming_distance_threshold) {
//                 if (pattern.shift == shift_is_pattern_size)
//                     set_shift_aa(p1_start + pattern.pattern.size());
//                 else
//                     set_shift_aa(p1_start + *pattern.shift);
//                 type_subtype_ = "A(H3N2)";
//                 // fmt::print(stderr, "H3 ({}) {}\n{}\n", hamd, debug_name, std::string_view(aa_.data(), 200));
//                 return true;
//             }
//             else
//                 hamds.push_back(hamd);
//         }
//     }

//     // if (!hamds.empty())
//     //     fmt::print(stderr, "NOT H3? ({}) {}\n{}\n", hamds, debug_name, std::string_view(aa_.data(), 200));
//     return false;

//     // const std::string_view pattern_1{sH3N2_align_1};

//     // const std::string_view look_for = pattern_1.substr(0, 2);
//     // std::vector<size_t> hamds;
//     // for (auto p1_start = data.find(look_for); p1_start < (data.size() - pattern_1.size()); p1_start = data.find(look_for, p1_start + look_for.size())) {
//     //     if (const auto hamd = ::string::hamming_distance(pattern_1, data.substr(p1_start, pattern_1.size())); hamd < 6) {
//     //         shift_aa_ = p1_start;
//     //         shift_nuc_ = nuc_translation_offset_ + p1_start * 3;
//     //         return true;
//     //     }
//     //     else
//     //         hamds.push_back(hamd);
//     // }

//     // fmt::print(stderr, "NOT H3? ({}) {}\n{}\n", hamds, debug_name, aa_); // std::string_view(aa_.data(), 100));
//     // return false;

//     // if (aa_start == std::string::npos) {
//     //     // fmt::print(stderr, "NOT H3? {}\n{}\n", debug_name, std::string_view(aa_.data(), 100));
//     //     return false;
//     // }
//     // const auto hamd = ::string::hamming_distance(pattern_1, std::string_view(aa_.data() + aa_start, pattern_1.size()));
//     // if (hamd > 10)
//     //     fmt::print(stderr, "hamm {} {}\n{}\n", hamd, debug_name, std::string_view(aa_.data(), 120));


// } // acmacs::seqdb::v3::sequence_t::align_h3n2

// // ----------------------------------------------------------------------

// bool acmacs::seqdb::v3::sequence_t::align_any(std::string_view debug_name, std::string_view except)
// {
//     // if (except != "A(H3N2)" && align_h3n2(debug_name))
//     //     return true;
//     return false;

// } // acmacs::seqdb::v3::sequence_t::align_any

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
