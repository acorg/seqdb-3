#include "acmacs-base/enumerate.hh"
#include "acmacs-base/range.hh"
#include "acmacs-base/range-v3.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

namespace local
{
    class sequence_no_found : public std::runtime_error { public: sequence_no_found() : std::runtime_error{"sequence not found in seqdb"} {} };

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

    template <typename Same, typename Range, typename... Ranges> inline bool not_all_same(Same same_f, Range&& rng, Ranges&&... rngs)
    {
        if (!ranges::all_of(rng, same_f))
            return true;
        if constexpr (sizeof...(rngs) > 0)
            return not_all_same(same_f, std::forward<Ranges>(rngs)...);
        return false;
    }

    template <typename ... Ranges> inline std::vector<size_t> positions_with_differences(std::string_view master_sequence, Ranges&& ... rngs)
    {
        std::vector<size_t> result;
        for (size_t pos = 0; pos < master_sequence.size(); ++pos) {
            const auto aa_same = [aa = master_sequence[pos], pos](std::string_view seq) { return seq.size() <= pos || seq[pos] == aa; };
            if (not_all_same(aa_same, std::forward<Ranges>(rngs)...))
                result.push_back(pos);
        }
        return result;
    }

    // ----------------------------------------------------------------------

    inline std::string_view aligned(const auto& ref, enum acmacs::seqdb::compare cmp_nuc_aa)
    {
        if (!ref)
            throw sequence_no_found{};
        switch (cmp_nuc_aa) {
            case acmacs::seqdb::compare::nuc:
                return ref.nuc_aligned(acmacs::seqdb::get());
            case acmacs::seqdb::compare::aa:
                return ref.aa_aligned(acmacs::seqdb::get());
        }
    };

    constexpr const char* compare_report_text_split_space = "  ";

    template <typename R1, typename R2> std::string compare_report_text_header(const acmacs::seqdb::v3::ref& master, std::string_view title1, R1&& r1, std::string_view title2, R2&& r2)
    {
        fmt::memory_buffer out, col_labels;
        fmt::format_to(col_labels, "     0");
        if (!title1.empty())
            fmt::format_to(out, "{}\n", title1);
        fmt::format_to(out, " 0 {}\n", master.seq_id());
        size_t col_no = 1;
        for (const auto& ref : r1) {
            fmt::format_to(out, "{:2d} {}\n", col_no, ref.seq_id());
            fmt::format_to(col_labels, " {:2d}", col_no);
            ++col_no;
        }
        fmt::format_to(out, "\n");
        fmt::format_to(col_labels, compare_report_text_split_space);
        if (!title2.empty())
            fmt::format_to(out, "{}\n", title2);
        for (const auto& ref : r2) {
            fmt::format_to(out, "{:2d} {}\n", col_no, ref.seq_id());
            fmt::format_to(col_labels, " {:2d}", col_no);
            ++col_no;
        }
        fmt::format_to(out, "\n{}\n", fmt::to_string(col_labels));
        return fmt::to_string(out);
    }

    template <typename R1, typename R2> std::string compare_report_text(const acmacs::seqdb::v3::ref& master, std::string_view title1, enum acmacs::seqdb::compare cmp_nuc_aa, R1&& r1, std::string_view title2, R2&& r2)
    {
        if (!master)
            throw sequence_no_found{};

        fmt::memory_buffer out;

        const auto aligned = [cmp_nuc_aa](const auto& ref) -> std::string_view { return local::aligned(ref, cmp_nuc_aa); };
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

        return compare_report_text_header(master, title1, std::forward<R1>(r1), title2, std::forward<R2>(r2)) + fmt::to_string(out);
    }

} // namespace local

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subset& sequences, enum compare cmp_nuc_aa, size_t split)
{
    if (split > 0)
        return local::compare_report_text(sequences[0], "", cmp_nuc_aa, ranges::views::drop(sequences, 1) | ranges::views::take(split - 1), "", ranges::views::drop(sequences, static_cast<ssize_t>(split)));
    else
        return local::compare_report_text(sequences[0], "", cmp_nuc_aa, ranges::views::drop(sequences, 1), "", ranges::empty_view<ref>{});

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subsets_by_title_t& subsets, enum compare cmp_nuc_aa)
{
    if (subsets.size() == 1) {
        const auto& [title, subset] = *subsets.begin();
        return local::compare_report_text(subset[0], title, cmp_nuc_aa, ranges::views::drop(subset, 1), "", ranges::empty_view<ref>{});
    }
    else if (subsets.size() == 2) {
        const auto& [title1, subset1] = *subsets.begin();
        const auto& [title2, subset2] = *std::next(subsets.begin());
        return local::compare_report_text(subset1[0], title1, cmp_nuc_aa, ranges::views::drop(subset1, 1), title2, ranges::views::all(subset2));
    }
    else {
        throw std::runtime_error("too many or too few subsets (supported 1 or 2 at the time being)");
    }

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subset& set1, const subset& set2, enum compare cmp_nuc_aa)
{
    return local::compare_report_text(set1[0], "", cmp_nuc_aa, set1 | ranges::views::drop(1), "", ranges::views::all(set2));

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

namespace local
{
    extern const char* sReportHtml;

    template <typename R1, typename R2>
    std::string compare_report_html(std::string_view title, const acmacs::seqdb::v3::ref& master, std::string_view title1, enum acmacs::seqdb::compare cmp_nuc_aa, R1&& r1, std::string_view title2,
                                    R2&& r2)
    {
        if (!master)
            throw sequence_no_found{};

        const auto aligned = [cmp_nuc_aa](const auto& ref) -> std::string_view { return local::aligned(ref, cmp_nuc_aa); };
        const std::string_view master_seq = aligned(master);
        const auto seq1 = ranges::to<std::vector<std::string_view>>(r1 | ranges::views::transform(aligned));
        const auto seq2 = ranges::to<std::vector<std::string_view>>(r2 | ranges::views::transform(aligned));
        const auto positions = positions_with_differences(master_seq, ranges::views::all(seq1), ranges::views::all(seq2));

        fmt::memory_buffer rows;
        fmt::format_to(rows, "<tr><td class='fasta-aa-name'></td>", master.seq_id());
        for (auto pos : positions)
            fmt::format_to(rows, "<td class='pos-no'>{}</td>", pos + 1);
        fmt::format_to(rows, "</tr>\n");

        if (!title1.empty())
            fmt::format_to(rows, "<tr class='group-title'><td colspan=1000>{}</td></tr>\n", title1);

        fmt::format_to(rows, "<tr><td class='fasta-aa-name'>{}</td>", master.seq_id());
        const auto nuc_aa = cmp_nuc_aa == acmacs::seqdb::compare::aa ? 'a' : 'n';
        for (auto pos : positions)
            fmt::format_to(rows, "<td class='aa' title='{pos} {aa} {seq_id}' {nuc_aa}{aa}>{aa}</td>", fmt::arg("nuc_aa", nuc_aa), fmt::arg("aa", master_seq[pos]), fmt::arg("pos", pos + 1), fmt::arg("seq_id", master.seq_id()));
        fmt::format_to(rows, "</tr>\n");

        const auto make_row = [&](const auto& ref) {
            fmt::format_to(rows, "<tr><td class='fasta-aa-name'>{}</td>", ref.seq_id());
            const auto seq = aligned(ref);
            for (auto pos : positions)
                fmt::format_to(rows, "<td class='aa' title='{} {} {}' {}{}>{}</td>", pos + 1, seq[pos], ref.seq_id(), nuc_aa, seq[pos], seq_aa(master_seq[pos], seq, pos));
            fmt::format_to(rows, "</tr>\n");
        };

        ranges::for_each(r1, make_row);
        fmt::format_to(rows, "<tr class='separator'><td colspan=1000></td></tr>\n");
        if (!title2.empty())
            fmt::format_to(rows, "<tr class='group-title'><td colspan=1000>{}</td></tr>\n", title2);
        ranges::for_each(r2, make_row);

        const auto data = fmt::format("<table class='fasta-aa-entries'>{}</table>\n", fmt::to_string(rows));
        return fmt::format(sReportHtml, fmt::arg("title", title), fmt::arg("data", data));
    }
} // namespace local

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(std::string_view title, const subsets_by_title_t& subsets, enum compare cmp_nuc_aa)
{
    if (subsets.size() == 1) {
        const auto& [title1, subset1] = *subsets.begin();
        return local::compare_report_html(title, subset1[0], title1, cmp_nuc_aa, ranges::views::drop(subset1, 1), "", ranges::empty_view<ref>{});
    }
    else if (subsets.size() == 2) {
        const auto& [title1, subset1] = *subsets.begin();
        const auto& [title2, subset2] = *std::next(subsets.begin());
        return local::compare_report_html(title, subset1[0], title1, cmp_nuc_aa, ranges::views::drop(subset1, 1), title2, ranges::views::all(subset2));
    }
    else {
        throw std::runtime_error("too many or too few subsets (supported 1 or 2 at the time being)");
    }

} // acmacs::seqdb::v3::compare_report_html

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(std::string_view title, const subset& sequences, enum compare cmp_nuc_aa, size_t split)
{
    if (split > 0)
        return local::compare_report_html(title, sequences[0], "", cmp_nuc_aa, ranges::views::drop(sequences, 1) | ranges::views::take(split - 1), "", ranges::views::drop(sequences, static_cast<ssize_t>(split)));
    else
        return local::compare_report_html(title, sequences[0], "", cmp_nuc_aa, ranges::views::drop(sequences, 1), "", ranges::empty_view<ref>{});

} // acmacs::seqdb::v3::compare_report_html

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_html(std::string_view title, const subset& set1, const subset& set2, enum compare cmp_nuc_aa)
{
    return local::compare_report_html(title, set1[0], "", cmp_nuc_aa, set1 | ranges::views::drop(1), "", ranges::views::all(set2));

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
     table.fasta-aa-entries tr.separator td {{ height: 1em; border-bottom: 1px solid #A0A0A0; }}
     table.fasta-aa-entries tr.group-title td {{ color: #FF8000; }}

     /* Ana and Lily colors of 2020-05-13 Subject: [Lab] Colour scheme for colouring amino acids by property group */
     [aA] {{ color: #A30000; }} /* hydrophobic */
     [aC] {{ color: #FFFF29; }} /* special */
     [aD] {{ color: #4CBB17; }} /* negative */
     [aE] {{ color: #1A4F30; }} /* negative */
     [aF] {{ color: #681226; }} /* hydrophobic */
     [aG] {{ color: #F7CF08; }} /* special */
     [aH] {{ color: #7571BC; }} /* positive */
     [aI] {{ color: #FF5A3D; }} /* hydrophobic */
     [aK] {{ color: #440FAE; }} /* positive */
     [aL] {{ color: #C34822; }} /* hydrophobic */
     [aM] {{ color: #EB1919; }} /* hydrophobic */
     [aN] {{ color: #18E0F7; }} /* polar uncharged */
     [aP] {{ color: #CD8D00; }} /* special */
     [aQ] {{ color: #0080FF; }} /* polar uncharged */
     [aR] {{ color: #B2A1FF; }} /* positive */
     [aS] {{ color: #0F52BA; }} /* polar uncharged */
     [aT] {{ color: #000080; }} /* polar uncharged */
     [aV] {{ color: #F9BDED; }} /* hydrophobic */
     [aW] {{ color: #C32290; }} /* hydrophobic */
     [aY] {{ color: #FF2976; }} /* hydrophobic */

     [aX] {{ color: grey; }}
     p.table-title {{ font-weight: bold; margin: 3em 0 0 2em; }}

     [nA] {{ color: #0080FF; }}
     [nC] {{ color: #4CBB17; }}
     [nG] {{ color: #EB1919; }}
     [nT] {{ color: #CD8D00; }}

    </style>
    <title>{title}</title>
  </head>
  <body>
    <h2>{title}</h2>
{data}
    <br>
    <br>
    <br>
    <table class="all-letters">
      <tr>
        <td aA>A</td>
        <td aC>C</td>
        <td aD>D</td>
        <td aE>E</td>
        <td aF>F</td>
        <td aG>G</td>
        <td aH>H</td>
        <td aI>I</td>
        <td aK>K</td>
        <td aL>L</td>
        <td aM>M</td>
        <td aN>N</td>
        <td aP>P</td>
        <td aQ>Q</td>
        <td aR>R</td>
        <td aS>S</td>
        <td aT>T</td>
        <td aV>V</td>
        <td aW>W</td>
        <td aY>Y</td>
        <td aX>X</td>
      </tr>
      <tr>
        <td nA>A</td>
        <td nC>C</td>
        <td nG>G</td>
        <td nT>T</td>
      </tr>
    </table>
  </body>
</html>
)";

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
