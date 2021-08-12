#include "acmacs-base/string.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/range.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/color-amino-acid.hh"
#include "acmacs-base/string-compare.hh"
#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/color-amino-acid.hh"
#include "acmacs-base/html.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

namespace local
{
    class sequence_no_found : public std::runtime_error
    {
      public:
        sequence_no_found() : std::runtime_error{"sequence not found in seqdb"} {}
    };

} // namespace local

// ----------------------------------------------------------------------

std::vector<acmacs::seqdb::v3::pos0_t> acmacs::seqdb::v3::subset_to_compare_base_t::positions_to_report() const
{
    std::vector<pos0_t> positions;
    for (size_t pos{0}; pos < counters.size(); ++pos) {
        // AD_DEBUG("pos:{:3d} {:2d} {}", pos + 1, counters[pos].size(), counters[pos].report_sorted_max_first());
        if (counters[pos].size() > 1)
            positions.push_back(pos0_t{pos});
    }
    return positions;

} // acmacs::seqdb::v3::subset_to_compare_base_t::positions_to_report

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset_to_compare_base_t::most_frequent(const std::vector<pos0_t>& positions) const
{
    std::string result(positions.size(), '/');
    for (size_t pp{0}; pp < positions.size(); ++pp)
        result[pp] = counters[*positions[pp]].max().first;
    return result;

} // acmacs::seqdb::v3::subset_to_compare_base_t::most_frequent

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset_to_compare_base_t::format_summary(const std::vector<pos0_t>& positions, std::string_view prefix, size_t name_width, size_t column_width, const std::string* most_frequent) const
{
    const auto num_rows{max_counter_size()};
    fmt::memory_buffer output;
    for (size_t row_no{0}; row_no < num_rows; ++row_no) {
        if (row_no == 0)
            fmt::format_to_mb(output, "{}{:{}s}", prefix, name, name_width);
        else
            fmt::format_to_mb(output, "{}{:{}c}", prefix, ' ', name_width);
        for (size_t pp{0}; pp < positions.size(); ++pp) {
            const auto pos{positions[pp]};
            if (row_no == 0 && most_frequent) {
                if (const auto aa{counters[*pos].sorted()[row_no]}; aa != (*most_frequent)[pp])
                    fmt::format_to_mb(output, "{:^{}c}", aa, column_width);
                else
                    fmt::format_to_mb(output, "{:^{}c}", '.', column_width);
            }
            else if (const auto aa{counters[*pos].sorted()}; row_no < aa.size())
                fmt::format_to_mb(output, "{:^{}c}", aa[row_no], column_width);
            else
                fmt::format_to_mb(output, "{:^{}c}", ' ', column_width);
        }
        fmt::format_to_mb(output, "\n");
    }
    return fmt::to_string(output);

} // acmacs::seqdb::v3::subset_to_compare_base_t::format_summary

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset_to_compare_base_t::format_summary(const std::vector<pos0_t>& positions, std::string_view prefix, size_t name_width, size_t column_width, const std::string* most_frequent, double threshold) const
{
    fmt::memory_buffer output;
    fmt::format_to_mb(output, "{}{:{}s}", prefix, name, name_width);
    for (size_t pp{0}; pp < positions.size(); ++pp) {
        const auto pos{positions[pp]};
        const auto aa_pairs{counters[*pos].pairs(counter_t::sorted::yes)};
        const auto total{static_cast<double>(std::accumulate(std::begin(aa_pairs), std::end(aa_pairs), 0ul, [](size_t sum, const auto& ap) { return sum + ap.second; }))};
        std::string aas{aa_pairs.front().first, ' ', ' '};
        size_t offset{1};
        for (auto aapp{std::next(std::begin(aa_pairs))}; aapp != std::end(aa_pairs); ++aapp) {
            if ((static_cast<double>(aapp->second) / total) > threshold)
                aas[offset++] = aapp->first;
        }
        if (most_frequent)
            ::string::replace_in_place(aas, (*most_frequent)[pp], '.');
        fmt::format_to_mb(output, "{:^{}s}", aas, column_width);
    }
    fmt::format_to_mb(output, "\n");
    return fmt::to_string(output);

} // acmacs::seqdb::v3::subset_to_compare_base_t::format_summary

// ======================================================================

void acmacs::seqdb::v3::subset_to_compare_t::make_counters(enum compare cmp_nuc_aa)
{
    for (const auto& ref : subset) {
        const auto seq{aligned(ref, cmp_nuc_aa)};
        if (counters.size() < *seq.size())
            counters.resize(*seq.size());
        for (pos0_t pos{0}; pos < seq.size(); ++pos)
            counters[*pos].count(seq.at(pos));
    }

} // acmacs::seqdb::v3::subset_to_compare_t::make_counters

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset_to_compare_t::format_seq_ids(size_t indent) const
{
    fmt::memory_buffer output;
    fmt::format_to_mb(output, "{:{}c}{}\n", ' ', indent, name);
    for (const auto& ref : subset) {
        fmt::format_to_mb(output, "{:{}c}{}\n", ' ', indent + 2, ref.seq_id());
        // fmt::format_to_mb(output, "{:{}c}{}\n", ' ', indent + 2, local::aligned(ref, compare::aa));
    }
    return fmt::to_string(output);

} // acmacs::seqdb::v3::subset_to_compare_t::format_seq_ids

// ----------------------------------------------------------------------

acmacs::seqdb::sequence_aligned_t acmacs::seqdb::v3::subset_to_compare_t::aligned(const acmacs::seqdb::ref& ref, enum acmacs::seqdb::compare cmp_nuc_aa)
{
    if (!ref)
        throw local::sequence_no_found{};
    switch (cmp_nuc_aa) {
        case acmacs::seqdb::compare::nuc:
            return acmacs::seqdb::sequence_aligned_t{ref.nuc_aligned(acmacs::seqdb::get())};
        case acmacs::seqdb::compare::aa:
            return acmacs::seqdb::sequence_aligned_t{ref.aa_aligned(acmacs::seqdb::get())};
    }
    AD_ERROR("unreachable code");
    throw std::runtime_error{"unreachable code"}; // hey g++9

} // acmacs::seqdb::v3::subset_to_compare_t::aligned

// ======================================================================

void acmacs::seqdb::v3::subset_to_compare_selected_t::make_counters(enum compare cmp_nuc_aa)
{
    for (const auto [ag_no, ag] : selected) {
        acmacs::seqdb::sequence_aligned_t seq{cmp_nuc_aa == acmacs::seqdb::compare::nuc ? ag->sequence_nuc() : ag->sequence_aa()};
        if (counters.size() < *seq.size())
            counters.resize(*seq.size());
        for (pos0_t pos{0}; pos < seq.size(); ++pos)
            counters[*pos].count(seq.at(pos));
    }

} // acmacs::seqdb::v3::subset_to_compare_selected_t::make_counters

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset_to_compare_selected_t::format_seq_ids(size_t indent) const
{
    fmt::memory_buffer output;
    fmt::format_to_mb(output, "{:{}c}{}\n", ' ', indent, name);
    for (const auto [ag_no, ag] : selected)
        fmt::format_to_mb(output, "{:{}c}{}\n", ' ', indent + 2, ag->name_full());
    return fmt::to_string(output);

} // acmacs::seqdb::v3::subset_to_compare_selected_tt::format_seq_ids

// ----------------------------------------------------------------------

acmacs::seqdb::sequence_aligned_t acmacs::seqdb::v3::subset_to_compare_selected_t::aligned(const acmacs::chart::Antigen& antigen, enum acmacs::seqdb::compare cmp_nuc_aa)
{
    switch (cmp_nuc_aa) {
        case acmacs::seqdb::compare::nuc:
            return sequence_aligned_t{antigen.sequence_nuc()};
        case acmacs::seqdb::compare::aa:
            return sequence_aligned_t{antigen.sequence_aa()};
    }
    throw std::runtime_error{"unreachable code"}; // hey g++9

} // acmacs::seqdb::v3::subset_to_compare_selected_t::aligned

// ======================================================================

void acmacs::seqdb::v3::update_common(sequence_aligned_t& target, const subset& source, enum compare cmp_nuc_aa)
{
    for (const auto& ref : source) {
        const auto seq{subset_to_compare_t::aligned(ref, cmp_nuc_aa)};
        const auto end{target.size()};
        if (end < seq.size()) {
            target.resize(seq.size());
            std::copy(std::next(seq->begin(), static_cast<ssize_t>(*end)), seq->end(), std::next(target.get().begin(), static_cast<ssize_t>(*end)));
        }
        for (acmacs::seqdb::pos0_t pos{0}; pos < end; ++pos) {
            if (seq.at(pos) != target.at(pos))
                target.set(pos, ' ');
        }
    }

} // acmacs::seqdb::v3::update_common

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::detail::generate_html(std::string_view html_filename, std::string_view data_filename_name, std::string_view data_var_name)
{
    using namespace std::string_view_literals;

    const auto templates_dir{fmt::format("{}/share/templates/seqdb-3", acmacs::acmacsd_root())};
    acmacs::html::Generator html;
    html.title("Compare sequences"sv);
    html.add_css(acmacs::amino_acid_nucleotide_color_css());
    html.add_css(static_cast<std::string>(acmacs::file::read(fmt::format("{}/compare-sequences.css", templates_dir))));
    html.add_script_link(data_filename_name);
    html.add_script(static_cast<std::string>(acmacs::file::read(fmt::format("{}/compare-sequences.js", templates_dir))));
    html.add_script(fmt::format(R"(document.addEventListener("DOMContentLoaded", function() {{ compare_sequences({}); }});)", data_var_name));
    html.add_to_body(static_cast<std::string>(acmacs::file::read(fmt::format("{}/compare-sequences.body.html", templates_dir))));
    acmacs::file::write(html_filename, html.generate());

} // acmacs::seqdb::v3::detail::generate_html

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
