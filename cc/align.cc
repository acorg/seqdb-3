#include "acmacs-base/named-type.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/fmt.hh"
#include "seqdb-3/align.hh"
#include "seqdb-3/hamming-distance.hh"

// http://signalpeptide.com

// ----------------------------------------------------------------------

namespace align_detail
{
    inline std::string_view type_subtype_hint(std::string_view type_subtype)
    {
        return type_subtype.size() > 5 && type_subtype[0] == 'A' ? type_subtype.substr(0, 5) : type_subtype;
    }

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
        const char* type_subtype_hint;
        char start_aa;
    };

    static constexpr const std::array start_aa_table{
        start_aa_t{"A(H1N", 'D'}, // DTIC, DTLC
        start_aa_t{"A(H2N", 'D'}, // DQIC
        start_aa_t{"A(H3N", 'Q'},
        start_aa_t{"A(H4N", 'Q'},
        start_aa_t{"A(H5N", 'D'},
        start_aa_t{"A(H6N", 'D'},
        start_aa_t{"A(H7N", 'D'},
        start_aa_t{"A(H8N", 'D'}, // DRIC
        start_aa_t{"A(H9N", 'D'}, // DKIC
        start_aa_t{"A(H10", 'D'}, // DRIC
        start_aa_t{"A(H11", 'D'}, // DEIC
        start_aa_t{"A(H12", 'D'}, // DKIC
        start_aa_t{"A(H13", 'D'}, // DRIC
        start_aa_t{"A(H14", 'Q'}, // QITN
        start_aa_t{"A(H15", 'D'}, // DKIC
        start_aa_t{"A(H16", 'D'}, // DKIC
        start_aa_t{"A(H17", 'D'}, // DRIC
        start_aa_t{"B",     'D'}, // DRIC
    };

    inline char start_aa(std::string_view hint)
    {
        if (const auto found = ranges::find_if(start_aa_table, [hint](const auto& entry) { return entry.type_subtype_hint == hint; }); found != ranges::end(start_aa_table))
            return found->start_aa;
        throw std::runtime_error(fmt::format("align_detail::start_aa: unsupported type_subtype: {}", hint));
    }

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

#pragma GCC diagnostic pop

}

// ----------------------------------------------------------------------

std::optional<std::tuple<int, std::string_view>> acmacs::seqdb::v3::align(std::string_view amino_acids, std::string_view type_subtype_hint)
{
    const auto make_type_subtype = [type_subtype_hint](std::string_view detected_type_subtype) -> std::string_view {
        if (type_subtype_hint.substr(0, 4) == detected_type_subtype.substr(0, 4))
            return type_subtype_hint;
        else
            return detected_type_subtype;
    };

    // --------------------------------------------------
    // first stage

    // H1
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MKV", "MKA"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 17, "DTLC"))
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H1)")};

    // H2
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MTIT", "MAII"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 14, "GDQIC"))
        return std::tuple{static_cast<int>(pos) + 15, make_type_subtype("A(H2)")};

    // H3
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MKTII"}); pos != std::string::npos && (amino_acids[pos + 16] == 'Q' || amino_acids[pos + 15] == 'A')) // amino_acids.substr(pos + 15, 2) != "DR") { // DR[ISV]C - start of the B sequence (signal peptide is 15 aas!)
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H3)")};

    // H4
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MLS"}); pos != std::string::npos && (amino_acids[pos + 16] == 'Q' || align_detail::has_infix(amino_acids, pos + 16, "SQNY")))
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H4)")};

    // H5
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MEKIV"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H5)")};

    // H6
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MIAIIV"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H6)")};

    // H7
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MNIQ", "MNNQ", "MNTQ"}); pos != std::string::npos && amino_acids[pos + 17] != 'S' && align_detail::has_infix(amino_acids, pos + 18, "DKIC")) // SDKIC is H15 most probably
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H7)")};

    // H8
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MEKFIA"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H8)")};

    // H9
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"METIS", "MEIIS", "MEV"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 17, "ADKIC"))
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H9)")};

    // H10
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MYK"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H10)")};

    // // H11
    // if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MEK"}); pos != std::string::npos)
    //     return std::tuple{static_cast<int>(pos) + 16, make_type_subtype("A(H11)")};

    // H12
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MEK"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 15, "AYDKIC"))
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H12)")};

    // H13
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MDI", "MAL", "MEV"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 17, "ADRIC"))
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H13)")};

    // H14
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MIA"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 14, "AYSQITN"))
        return std::tuple{static_cast<int>(pos) + 17, make_type_subtype("A(H14)")};

    // H15 - second stage only

    // H16
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MMVK", "MMIK"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 19, "DKIC"))
        return std::tuple{static_cast<int>(pos) + 19, make_type_subtype("A(H16)")};

    // H17
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"MEL"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 17, "GDRICI"))
        return std::tuple{static_cast<int>(pos) + 18, make_type_subtype("A(H17)")};

    // --------------------------------------------------
    // second stage

    // H4
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 100, {"QNYT"}); pos != std::string::npos && align_detail::has_infix(amino_acids, pos + 11, "GHHA"))
        return std::tuple{static_cast<int>(pos), make_type_subtype("A(H4)")};

    // H10
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 100, {"LDKICLGHHA"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 1, make_type_subtype("A(H10)")};

    // H15
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 100, {"KSDKICLGHHA"}); pos != std::string::npos)
        return std::tuple{static_cast<int>(pos) + 2, make_type_subtype("A(H15)")};

    // --------------------------------------------------
    // third stage

    // H1
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 100, {"DTLC"}); pos != std::string::npos && pos >= 17 && amino_acids[pos - 17] == 'M')
        return std::tuple{static_cast<int>(pos) - 5, make_type_subtype("A(H1)")};

    // H5
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 100, {"GYHA"}); pos != std::string::npos && pos >= 21 && align_detail::has_infix(amino_acids, pos - 5, "DQ") && amino_acids[pos - 21] == 'M')
        return std::tuple{static_cast<int>(pos) - 5, make_type_subtype("A(H5)")};

    // H7
    if (const auto pos = align_detail::find_in_sequence(amino_acids, 20, {"DKICLGHHAV"}); pos == 0) // sequence start
        return std::tuple{static_cast<int>(pos), make_type_subtype("A(H7)")};

    // // H3
    // if (const auto pos = align_detail::find_in_sequence(amino_acids, 100, {"GHHA"}); pos != std::string::npos &&
    //     // (pos < 16 || amino_acids[pos - 16] == 'Q') &&
    //     amino_acids[pos + 5] == 'P' && amino_acids[pos + 7] == 'G') //  && amino_acids[pos + 15] == 'D')
    //     return std::tuple{static_cast<int>(pos) - 16, make_type_subtype("A(H3)")};

    return std::nullopt;

} // acmacs::seqdb::v3::align

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::Aligner::update(std::string_view amino_acids, int shift, std::string_view type_subtype)
{
    const std::string hint{align_detail::type_subtype_hint(type_subtype)};
    tables_.try_emplace(hint).first->second.update(amino_acids, shift);

} // acmacs::seqdb::v3::Aligner::update

// ----------------------------------------------------------------------

std::optional<std::tuple<int, std::string_view>> acmacs::seqdb::v3::Aligner::align(std::string_view amino_acids, std::string_view type_subtype_hint) const
{
    const auto hint = align_detail::type_subtype_hint(type_subtype_hint);
    if (const auto found = tables_.find(hint); found != tables_.end()) {
        if (const auto res = found->second.align(align_detail::start_aa(hint), amino_acids); res.has_value())
            return std::tuple(*res, type_subtype_hint);
    }

    // for (const auto& [ts, table] : tables_) {
    //     if (ts != hint) {
    //         if (const auto res = found->second.align(align_detail::start_aa(ts), amino_acids, debug_name); res.has_value())
    //             return std::tuple(*res, ts);
    //     }
    // }

    return std::nullopt;

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

std::optional<int> acmacs::seqdb::v3::Aligner::table_t::align(char start_aa, std::string_view amino_acids) const
{
    // fmt::print(stderr, "Aligner::table_t::align {}\n{}\n", debug_name, amino_acids);
    for (auto p_start = amino_acids.find(start_aa); p_start < (amino_acids.size() / 2); p_start = amino_acids.find(start_aa, p_start + 1)) {
        const auto all_pos = ranges::view::iota(0UL, std::min(max_sequence_length, amino_acids.size() - p_start));
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

} // acmacs::seqdb::v3::Aligner::table_t::align

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

// ****************************************************************************************************

// namespace align_detail
// {
//     using max_offset_t = acmacs::named_size_t<struct max_offset_t_tag>;
//     using hdth_t = acmacs::named_size_t<struct hdth_t_tag>; // hamming_distance_threshold
//     using shift_t = acmacs::named_size_t<struct shift_t_tag>;
//     constexpr const shift_t shift_is_pattern_size{-1};

//     struct pat_t
//     {
//         std::string_view type_subtype_prefix;
//         std::string_view type_subtype;
//         char start_aa;
//         std::string_view pattern;
//         max_offset_t max_offset;
//         hdth_t hamming_distance_threshold;
//         shift_t shift;
//     };

// #pragma GCC diagnostic push
// #ifdef __clang__
// #pragma GCC diagnostic ignored "-Wglobal-constructors"
// #pragma GCC diagnostic ignored "-Wexit-time-destructors"
// #endif

//     static const std::array patterns{
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSYILCLVFA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com [50k]
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSCILCLVFA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A/New York/47/2003(H3N2) [2.2k]
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIAFSCILCLIFA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // A/NEW JERSEY/53/2015 (CDC) [1.1k]
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSHILCLVFA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A/New York/440/2000(H3N2) [745]
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSYILCMVFA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A/New York/152/2000(H3N2) [692]

//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIVLSCFFCLAFC", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N2)/blue-winged teal/Ohio/31/1999
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSYIFCLAFG", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N2)/Duck/Hong Kong/7/1975, A(H3N8)/duck/Chabarovsk/1610/1972
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSYVFCLAFG", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N6)/duck/Nanchang/8-174/2000
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIALSYIFCLALG", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N2)/Hong Kong/1-1-MA-12/1968 + 63 other results

//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTTIILILLTHWVYS", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N8)/equine/Idaho/37875/1991 + 70 other results
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTVIALSYILCLTFG", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N8)/Duck/Ukraine/1/1963
//         pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTTIVLILLTHWVYS", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com A(H3N8)/Equine/Kentucky/1/1987

//         // pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIAFSCILCQISA", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // A/SWINE/MANITOBA/D0083/2013
//         // pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTIIAFSCILCQISS", max_offset_t{50}, hdth_t{2}, shift_is_pattern_size}, // A/SWINE/MANITOBA/D0180/2012

//         // pat_t{"A(H3N", "A(H3N2)", 'Q', "MKTLIALSYIFCLVLG",                                                 max_offset_t{ 50}, hdth_t{2}, shift_is_pattern_size}, // signalpeptide.com
//         // A(H3N2)/Swine/Ukkel/1/1984 pat_t{"A(H3N", "A(H3N2)", 'Q', "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENC", max_offset_t{100}, hdth_t{6}, shift_t{0}}, pat_t{"A(H3N",
//         // "A(H3N2)", 'Q', "QKLPGNNNSTATLCLGHHAVPNGTIVKTI",                                    max_offset_t{100}, hdth_t{6}, shift_t{0}},
//     };

//     // signalpeptide H3
//     // MKTIIALCYILCLVFA
//     // MKTIIALSHIFCLVLG
//     // MKTIIALSYIFCLAFA
//     // MKTIIALSYIFCLAFS
//     // MKTIIALSYIFCLVFA
//     // MKTIIALSYIFCLVLG
//     // MKTIIALSYIFCQVFA
//     // MKTIIALSYIFCQVLA
//     // MKTIIALSYILCLVFA
//     // MKTIIALSYISCLVFA
//     // MKTIIVLSCFFCLAFS
//     // MKTIIVLSYFFCLALS
//     // MKTTIILILLIHWVHS
//     // MKTTTILILLTHWVHS

// #pragma GCC diagnostic pop

//     inline char start_aa(std::string_view type_subtype)
//     {
//         for (const auto& pattern : patterns) {
//             if (pattern.type_subtype_prefix == type_subtype || pattern.type_subtype == type_subtype)
//                 return pattern.start_aa;
//         }
//         throw std::runtime_error(fmt::format("align_detail::start_aa: unsupported type_subtype: {}", type_subtype));
//     }

// } // namespace align_detail

// std::optional<std::tuple<int, std::string_view>> acmacs::seqdb::v3::align(std::string_view amino_acids, std::string_view type_subtype_hint)
// {
//     const auto hint = align_detail::type_subtype_hint(type_subtype_hint);
//     const auto make_type_subtype = [hint, type_subtype_hint](const auto& pattern) -> std::string_view {
//         if (hint == pattern.type_subtype_prefix)
//             return type_subtype_hint;
//         else
//             return pattern.type_subtype;
//     };
//     const auto suitable_subtype = [hint](const auto& pattern) -> bool { return pattern.type_subtype_prefix == hint || hint.empty() || hint == "A(H0N"; };

//     for (const auto& pattern : align_detail::patterns | ranges::view::filter(suitable_subtype)) {
//         for (auto p_start = amino_acids.find(pattern.pattern[0]); p_start < *pattern.max_offset; p_start = amino_acids.find(pattern.pattern[0], p_start + 1)) {
//             if (hamming_distance(pattern.pattern, amino_acids.substr(p_start, pattern.pattern.size())) < *pattern.hamming_distance_threshold) {
//                 if (pattern.shift == align_detail::shift_is_pattern_size)
//                     return std::tuple{static_cast<int>(p_start + pattern.pattern.size()), make_type_subtype(pattern)};
//                 else
//                     return std::tuple{static_cast<int>(p_start + *pattern.shift), make_type_subtype(pattern)};
//             }
//         }
//     }
//     return std::nullopt;

// } // acmacs::seqdb::v3::align

// ****************************************************************************************************

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
