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
        return positions_with_differences(master_sequence, std::forward<R>(r1), ranges::views::empty<R>());
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
        const auto seq1 = ranges::to<std::vector<std::string_view>>(r1 | ranges::views::transform(aligned));
        const auto seq2 = ranges::to<std::vector<std::string_view>>(r2 | ranges::views::transform(aligned));

        for (auto pos : positions_with_differences(master_seq, ranges::views::all(seq1), ranges::views::all(seq2))) {
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
        return local::compare_report_text(sequences[0], ranges::views::drop(sequences, 1) | ranges::views::take(split - 1), ranges::views::drop(sequences, static_cast<ssize_t>(split)));
    else
        return local::compare_report_text(sequences[0], ranges::views::drop(sequences, 1), ranges::empty_view<ref>{});

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subset& set1, const subset& set2)
{
    return local::compare_report_text(set1[0], set1 | ranges::views::drop(1), ranges::views::all(set2));

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

namespace local
{
    extern const char* sReportHtml;

    template <typename R1, typename R2> std::string compare_report_html(std::string_view title, const acmacs::seqdb::v3::ref& master, R1&& r1, R2&& r2)
    {
        const auto aligned = [](const auto& ref) -> std::string_view { return ref.seq().aa_aligned(); };
        const std::string_view master_seq = aligned(master);
        const auto seq1 = ranges::to<std::vector<std::string_view>>(r1 | ranges::views::transform(aligned));
        const auto seq2 = ranges::to<std::vector<std::string_view>>(r2 | ranges::views::transform(aligned));
        const auto positions = positions_with_differences(master_seq, ranges::views::all(seq1), ranges::views::all(seq2));

        fmt::memory_buffer rows;
        fmt::format_to(rows, "<tr><td class='fasta-aa-name'></td>", master.seq_id());
        for (auto pos : positions)
            fmt::format_to(rows, "<td class='pos-no'>{}</td>", pos + 1);
        fmt::format_to(rows, "</tr>\n");

        fmt::format_to(rows, "<tr><td class='fasta-aa-name'>{}</td>", master.seq_id());
        for (auto pos : positions)
            fmt::format_to(rows, "<td class='aa' title='{pos} {aa} {seq_id}' a{aa}>{aa}</td>", fmt::arg("aa", master_seq[pos]), fmt::arg("pos", pos + 1), fmt::arg("seq_id", master.seq_id()));
        fmt::format_to(rows, "</tr>\n");

        const auto make_row = [&](const auto& ref) {
            fmt::format_to(rows, "<tr><td class='fasta-aa-name'>{}</td>", ref.seq_id());
            const auto seq = aligned(ref);
            for (auto pos : positions)
                fmt::format_to(rows, "<td class='aa' title='{} {} {}' a{}>{}</td>", pos + 1, seq[pos], ref.seq_id(), seq[pos], seq_aa(master_seq[pos], seq, pos));
            fmt::format_to(rows, "</tr>\n");
        };

        ranges::for_each(r1, make_row);
        fmt::format_to(rows, "<tr><td class='fasta-aa-name'>--</td></tr>\n");
        ranges::for_each(r2, make_row);

        const auto data = fmt::format("<table class='fasta-aa-entries'>{}</table>\n", fmt::to_string(rows));
        return fmt::format(sReportHtml, fmt::arg("title", title), fmt::arg("data", data));
    }
}

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(std::string_view title, const subset& sequences, size_t split)
{
    if (split > 0)
        return local::compare_report_html(title, sequences[0], ranges::views::drop(sequences, 1) | ranges::views::take(split - 1), ranges::views::drop(sequences, static_cast<ssize_t>(split)));
    else
        return local::compare_report_html(title, sequences[0], ranges::views::drop(sequences, 1), ranges::empty_view<ref>{});

} // acmacs::seqdb::v3::compare_report_html

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(std::string_view title, const subset& set1, const subset& set2)
{
    return local::compare_report_html(title, set1[0], set1 | ranges::views::drop(1), ranges::views::all(set2));

} // acmacs::seqdb::v3::compare_report_html

// ----------------------------------------------------------------------

const char* local::sReportHtml = R"(<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8" />
    <style>
     table.fasta-aa-entries {{ white-space: nowrap; }}
     table.fasta-aa-entries td.fasta-aa-name.fasta-aa-master {{ color: magenta; }}
     table.fasta-aa-entries .row-even {{ background-color: #F0F0F0; }}
     table.fasta-aa-entries table.fasta-aa-sequence {{ border-spacing: 0; font-family: Menlo, monospace; }}
     table.fasta-aa-entries td.fasta-aa-name {{ min-width: 30em; }}
     table.fasta-aa-entries td.pos-no {{ text-align: center; width: 2em; }}
     table.fasta-aa-entries td.aa {{ text-align: center; }}
     [aA] {{ color: blue; }}
     [aC] {{ color: salmon; }}
     [aD] {{ color: magenta; }}
     [aE] {{ color: magenta; }}
     [aF] {{ color: blue; }}
     [aG] {{ color: #FF8C00; }}
     [aH] {{ color: #008B8B; }}
     [aI] {{ color: blue; }}
     [aK] {{ color: red; }}
     [aL] {{ color: blue; }}
     [aM] {{ color: blue; }}
     [aN] {{ color: #228B22; }}
     [aP] {{ color: goldenrod; }}
     [aQ] {{ color: #228B22; }}
     [aR] {{ color: red; }}
     [aS] {{ color: #228B22; }}
     [aT] {{ color: #228B22; }}
     [aV] {{ color: blue; }}
     [aW] {{ color: blue; }}
     [aY] {{ color: #008B8B; }}
     [aX] {{ color: grey; }}
     p.table-title {{ font-weight: bold; margin: 3em 0 0 2em; }}
    </style>
    <title>{title}</title>
  </head>
  <body>
    <h2>{title}</h2>
{data}
  </body>
</html>
)";

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
