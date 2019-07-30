#include "acmacs-base/enumerate.hh"
#include "acmacs-base/range.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::compare_report_text(const subset& sequences)
{
    fmt::memory_buffer out;
    for (auto [no, ref] : acmacs::enumerate(sequences))
        fmt::format_to(out, "{:2d} {}\n", no, ref.seq_id());
    fmt::format_to(out, "\n   ");
    for (auto no : acmacs::range(sequences.size()))
        fmt::format_to(out, " {:>2d}", no);
    fmt::format_to(out, "\n");

    std::vector<std::string_view> seqs(sequences.size());
    std::transform(std::begin(sequences), std::end(sequences), std::begin(seqs), [](const auto& ref) { return ref.seq().aa_aligned(); });

    // std::for_each(std::next(std::begin(seqs)), std::end(seqs), [](std::string_view seq) { fmt::format_to(out, "{}\n\n", seq.size()); });

    for (size_t pos = 0; pos < seqs[0].size(); ++pos) {
        const auto aa = seqs[0][pos];
        if (!std::all_of(std::next(std::begin(seqs)), std::end(seqs), [aa, pos](std::string_view seq) {
                return seq.size() <= pos || seq[pos] == aa;
            })) {
            fmt::format_to(out, "{:3d} {:>2c}", pos + 1, aa);
            std::for_each(std::next(std::begin(seqs)), std::end(seqs), [aa, pos, &out](std::string_view seq) {
                char seq_aa = '-';
                if (seq.size() <= pos)
                    seq_aa = '-';
                else if (seq[pos] == aa)
                    seq_aa = '.';
                else
                    seq_aa = seq[pos];
                fmt::format_to(out, " {:>2c}", seq_aa);
            });
            fmt::format_to(out, "\n");
        }
    }

    return fmt::to_string(out);

} // acmacs::seqdb::v3::compare_report_text

// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
