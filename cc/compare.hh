#pragma once

#include "acmacs-base/counter.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    enum class compare { aa, nuc };

    struct subset_to_compare_t
    {
        using counter_t = acmacs::CounterCharSome<' ', '['>;
        using counters_t = std::vector<counter_t>;

        std::string name;
        acmacs::seqdb::subset subset;
        counters_t counters;

        subset_to_compare_t(std::string_view a_name) : name{a_name} {}
        subset_to_compare_t(std::string_view a_name, acmacs::seqdb::subset&& a_subset) : name{a_name}, subset{std::move(a_subset)} {}
        void make_counters(enum compare cmp_nuc_aa);
        size_t max_counter_size() const { return std::accumulate(std::begin(counters), std::end(counters), size_t{0}, [](size_t max, const auto& counter) { return std::max(max, counter.size()); }); }
        std::vector<pos0_t> positions_to_report() const;
        std::string most_frequent(const std::vector<pos0_t>& positions) const;
        std::string format_summary(const std::vector<pos0_t>& positions, std::string_view prefix, size_t name_width, size_t column_width, const std::string* most_frequent = nullptr) const;
        std::string format_summary(const std::vector<pos0_t>& positions, std::string_view prefix, size_t name_width, size_t column_width, const std::string* most_frequent, double threshold) const;
        std::string format_seq_ids(size_t indent) const;
    };

    struct subsets_to_compare_t
    {
        std::vector<subset_to_compare_t> subsets;
        enum compare cmp_nuc_aa;

        subsets_to_compare_t(enum compare nuc_aa) : cmp_nuc_aa{nuc_aa} {}
        void make_counters() { std::for_each(std::begin(subsets), std::end(subsets), [this](auto& ss) { ss.make_counters(cmp_nuc_aa); }); }
        std::vector<pos0_t> positions_to_report() const;

        size_t max_name() const { return std::accumulate(std::begin(subsets), std::end(subsets), 0ul, [](size_t max, const auto& ss) { return std::max(max, ss.name.size()); }); }
        std::string format_summary(size_t indent = 2, size_t column_width = 5, double threshold = -1.0) const;
        std::string format_seq_ids(size_t indent) const;
        std::string format_json(size_t indent) const;
    };

    void compare_sequences_generate_html(std::string_view html_filename, const subsets_to_compare_t& data);

    // ----------------------------------------------------------------------

    void update_common(sequence_aligned_t& target, const subset& source, enum compare cmp_nuc_aa);
    // sequence_aligned_t find_common(const subsets_to_compare_t& subsets, enum compare cmp_nuc_aa);

    template <typename... Subset> sequence_aligned_t find_common(const Subset&... subsets)
    {
        sequence_aligned_t target;
        for (const auto& ss : {(subsets, ...)}) {
            update_common(target, ss);
        }
        return target;
    }

} // namespace acmacs::seqdb::inline v3

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
