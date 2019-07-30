#include "acmacs-base/enumerate.hh"
#include "acmacs-base/range.hh"
#include "acmacs-base/range-v3.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

namespace local
{
    static inline char seq_aa(char aa, std::string_view seq, size_t pos)
    {
        if (seq.size() <= pos)
            return '-';
        else if (seq[pos] == aa)
            return '.';
        else
            return seq[pos];
    }

    // ----------------------------------------------------------------------

    template <typename R1, typename R2> inline std::vector<size_t> positions_with_differences(std::string_view master_sequence, R1&& r1, R2&& r2)
    {
        std::vector<size_t> result;
        for (size_t pos = 0; pos < master_sequence.size(); ++pos) {
            const auto aa_same = [aa = master_sequence[pos], pos](std::string_view seq) { return seq.size() <= pos || seq[pos] == aa; };
            if (!ranges::all_of(r1, aa_same) || !ranges::all_of(r2, aa_same))
                result.push_back(pos);
        }
        return result;
    }

    template <typename R> inline std::vector<size_t> positions_with_differences(std::string_view master_sequence, R&& r1)
    {
        return positions_with_differences(master_sequence, std::forward<R>(r1), ranges::view::empty<R>());
    }

    // ----------------------------------------------------------------------

    constexpr const char* compare_report_text_split_space = "  ";

    template <typename R1, typename R2> std::string compare_report_text_header(const acmacs::seqdb::v3::ref& master, R1&& r1, R2&& r2)
    {
        fmt::memory_buffer out, col_labels;
        fmt::format_to(col_labels, "     0");
        fmt::format_to(out, " 0 {}\n", master.seq_id());
        size_t col_no = 1;
        for (const auto& ref : r1) {
            fmt::format_to(out, "{:2d} {}\n", col_no, ref.seq_id());
            fmt::format_to(col_labels, " {:2d}", col_no);
            ++col_no;
        }
        fmt::format_to(out, "\n");
        fmt::format_to(col_labels, compare_report_text_split_space);
        for (const auto& ref : r2) {
            fmt::format_to(out, "{:2d} {}\n", col_no, ref.seq_id());
            fmt::format_to(col_labels, " {:2d}", col_no);
            ++col_no;
        }
        fmt::format_to(out, "\n{}\n", fmt::to_string(col_labels));
        return fmt::to_string(out);
    }

    template <typename R1, typename R2> std::string compare_report_text(const acmacs::seqdb::v3::ref& master, R1&& r1, R2&& r2)
    {
        fmt::memory_buffer out;

        const auto aligned = [](const auto& ref) -> std::string_view { return ref.seq().aa_aligned(); };
        const std::string_view master_seq = aligned(master);
        const auto seq1 = ranges::to<std::vector<std::string_view>>(r1 | ranges::view::transform(aligned));
        const auto seq2 = ranges::to<std::vector<std::string_view>>(r2 | ranges::view::transform(aligned));

        for (auto pos : positions_with_differences(master_seq, ranges::view::all(seq1), ranges::view::all(seq2))) {
            const auto aa = master_seq[pos];
            fmt::format_to(out, "{:3d} {:>2c}", pos + 1, aa);
            const auto format_aa = [aa, pos, &out](std::string_view seq) { fmt::format_to(out, " {:>2c}", seq_aa(aa, seq, pos)); };
            ranges::for_each(seq1, format_aa);
            fmt::format_to(out, compare_report_text_split_space);
            ranges::for_each(seq2, format_aa);
            fmt::format_to(out, "\n");
        }

        return compare_report_text_header(master, std::forward<R1>(r1), std::forward<R2>(r2)) + fmt::to_string(out);
    }

} // namespace local

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subset& sequences, size_t split)
{
    if (split > 0)
        return local::compare_report_text(sequences[0], ranges::view::drop(sequences, 1) | ranges::view::take(split - 1), ranges::view::drop(sequences, static_cast<ssize_t>(split)));
    else
        return local::compare_report_text(sequences[0], ranges::view::drop(sequences, 1), ranges::empty_view<ref>{});

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subset& set1, const subset& set2)
{
    return local::compare_report_text(set1[0], set1 | ranges::view::drop(1), ranges::view::all(set2));

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(const subset& sequences, size_t split)
{
    return {};

} // acmacs::seqdb::v3::compare_report_html

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(const subset& set1, const subset& set2)
{
    return "<html></html>";

} // acmacs::seqdb::v3::compare_report_html

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
