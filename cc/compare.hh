#pragma once

#include "acmacs-base/counter.hh"
#include "acmacs-base/to-json.hh"
#include "acmacs-base/read-file.hh"
#include "seqdb-3/seqdb.hh"
#include "acmacs-chart-2/selected-antigens-sera.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    enum class compare { aa, nuc };

    struct subset_to_compare_base_t
    {
        using counter_t = acmacs::CounterCharSome<' ', '['>;
        using counters_t = std::vector<counter_t>;

        std::string name;
        counters_t counters;

        subset_to_compare_base_t(std::string_view a_name) : name{a_name} {}
        subset_to_compare_base_t(const subset_to_compare_base_t&) = default;
        virtual ~subset_to_compare_base_t() = default;
        subset_to_compare_base_t& operator=(const subset_to_compare_base_t&) = default;

        virtual void make_counters(enum compare cmp_nuc_aa) = 0;
        size_t max_counter_size() const
        {
            return std::accumulate(std::begin(counters), std::end(counters), size_t{0}, [](size_t max, const auto& counter) { return std::max(max, counter.size()); });
        }
        std::vector<pos0_t> positions_to_report() const;
        std::string most_frequent(const std::vector<pos0_t>& positions) const;
        std::string format_summary(const std::vector<pos0_t>& positions, std::string_view prefix, size_t name_width, size_t column_width, const std::string* most_frequent = nullptr) const;
        std::string format_summary(const std::vector<pos0_t>& positions, std::string_view prefix, size_t name_width, size_t column_width, const std::string* most_frequent, double threshold) const;
    };

    struct subset_to_compare_t : public subset_to_compare_base_t
    {
        acmacs::seqdb::subset subset;

        using subset_to_compare_base_t::subset_to_compare_base_t;
        subset_to_compare_t(std::string_view a_name, acmacs::seqdb::subset&& a_subset) : subset_to_compare_base_t{a_name}, subset{std::move(a_subset)} {}

        void make_counters(enum compare cmp_nuc_aa) override;
        std::string format_seq_ids(size_t indent) const;

        auto begin() const { return subset.begin(); }
        auto end() const { return subset.end(); }
        bool empty() const { return subset.empty(); }

        static acmacs::seqdb::sequence_aligned_t aligned(const acmacs::seqdb::ref& ref, enum acmacs::seqdb::compare cmp_nuc_aa);
        static std::string seq_id(const acmacs::seqdb::ref& ref) { return *ref.seq_id(); }
    };

    struct subset_to_compare_selected_t : public subset_to_compare_base_t
    {
        acmacs::chart::SelectedAntigensModify selected;

        using subset_to_compare_base_t::subset_to_compare_base_t;
        subset_to_compare_selected_t(std::string_view a_name, acmacs::chart::SelectedAntigensModify&& a_selected) : subset_to_compare_base_t{a_name}, selected{std::move(a_selected)} {}
        subset_to_compare_selected_t(std::string_view a_name, const acmacs::chart::SelectedAntigensModify& a_selected) : subset_to_compare_base_t{a_name}, selected{a_selected} {}

        void make_counters(enum compare cmp_nuc_aa) override;
        std::string format_seq_ids(size_t indent) const;

        auto begin() const { return selected.begin(); }
        auto end() const { return selected.end(); }
        bool empty() const { return selected.empty(); }

        static acmacs::seqdb::sequence_aligned_t aligned(const acmacs::chart::Antigen& antigen, enum acmacs::seqdb::compare cmp_nuc_aa);
        static inline acmacs::seqdb::sequence_aligned_t aligned(const std::pair<size_t, std::shared_ptr<acmacs::chart::AntigenModify>>& antigen, enum acmacs::seqdb::compare cmp_nuc_aa)
        {
            return aligned(*antigen.second, cmp_nuc_aa);
        }
        static std::string seq_id(const std::pair<size_t, std::shared_ptr<acmacs::chart::AntigenModify>>& ref) { return ref.second->name_full(); }
    };

    // ----------------------------------------------------------------------

    template <typename Subset> struct subsets_to_compare_t
    {
        std::vector<Subset> subsets;
        enum compare cmp_nuc_aa;

        subsets_to_compare_t(enum compare nuc_aa) : cmp_nuc_aa{nuc_aa} {}
        void make_counters()
        {
            std::for_each(std::begin(subsets), std::end(subsets), [this](auto& ss) { ss.make_counters(cmp_nuc_aa); });
        }

        std::vector<pos0_t> positions_to_report() const
        {
            subset_to_compare_base_t::counters_t merged_counters; // (subsets.front().counters.size());
            for (const auto& ssc : subsets) {
                if (!ssc.empty()) {
                    if (merged_counters.size() < ssc.counters.size())
                        merged_counters.resize(ssc.counters.size());
                    for (size_t pos{0}; pos < ssc.counters.size(); ++pos)
                        merged_counters[pos] = merge_CounterCharSome(merged_counters[pos], ssc.counters[pos]);
                }
                else
                    AD_WARNING("subset empty: {}", ssc.name);
            }

            std::vector<pos0_t> positions;
            for (size_t pos{0}; pos < merged_counters.size(); ++pos) {
                if (merged_counters[pos].size() > 1)
                    positions.push_back(pos0_t{pos});
            }
            return positions;
        }

        size_t max_name() const
        {
            return std::accumulate(std::begin(subsets), std::end(subsets), 0ul, [](size_t max, const auto& ss) { return std::max(max, ss.name.size()); });
        }

        std::string format_summary(size_t indent = 2, size_t column_width = 5, double threshold = -1.0) const
        {
            const std::string prefix(indent, ' ');
            const auto name_width{max_name()};
            const auto positions{positions_to_report()};
            fmt::memory_buffer output;
            const auto most_frequent{subsets.front().most_frequent(positions)};
            // fmt::format_to_mb(output, "{}\n\n", most_frequent);
            fmt::format_to_mb(output, "{}{:{}c}", prefix, ' ', name_width);
            for (const auto pos : positions)
                fmt::format_to_mb(output, "{:^{}d}", pos, column_width);
            fmt::format_to_mb(output, "\n");
            bool first_group{true};
            for (const auto& ssc : subsets) {
                if (threshold > 0.0)
                    fmt::format_to_mb(output, "{}", ssc.format_summary(positions, prefix, name_width, column_width, first_group ? nullptr : &most_frequent, threshold));
                else
                    fmt::format_to_mb(output, "{}", ssc.format_summary(positions, prefix, name_width, column_width, first_group ? nullptr : &most_frequent));
                first_group = false;
            }
            return fmt::to_string(output);
        }

        std::string format_seq_ids(size_t indent) const
        {
            fmt::memory_buffer output;
            for (const auto& ssc : subsets)
                fmt::format_to_mb(output, "{}\n", ssc.format_seq_ids(indent));
            return fmt::to_string(output);
        }

        std::string format_json(size_t indent) const
        {
            using namespace to_json;
            using namespace std::string_view_literals;

            const auto positions{positions_to_report()};

            const auto make_group_pos = [&positions](const subset_to_compare_base_t& group) {
                const auto make_aa_counter = [](const auto& aap) { return object{key_val{"a"sv, std::string(1, aap.first)}, key_val{"c", aap.second}}; };
                object result;
                for (const auto pos : positions) {
                    const auto aa_pairs{group.counters[*pos].pairs(subset_to_compare_base_t::counter_t::sorted::yes)};
                    result << key_val{fmt::format("{}", pos), array{std::begin(aa_pairs), std::end(aa_pairs), make_aa_counter, array::compact_output::yes}};
                }
                return result;
            };

            const auto make_group = [this, make_group_pos](const Subset& group) {
                const auto make_sequence = [this](const auto& ref) {
                    return object{key_val{"id"sv, Subset::seq_id(ref)}, key_val{"seq"sv, fmt::format("{}", Subset::aligned(ref, cmp_nuc_aa))}};
                };
                return object{
                    key_val{"name"sv, group.name},                                      //
                    key_val{"pos1"sv, make_group_pos(group)},                           //
                    key_val{"seq"sv, array{group.begin(), group.end(), make_sequence}}, //
                };
            };

            const object data{
                key_val{"pos1"sv, array{std::begin(positions), std::end(positions), [](pos0_t pos) { return *pos + 1; }, array::compact_output::yes}}, //
                key_val{"groups"sv, array{subsets.begin(), subsets.end(), make_group}}                                                                 //
            };
            return fmt::format(fmt::runtime(fmt::format("{{:{}}}", indent)), data);
        }
    };

    namespace detail
    {
        void generate_html(std::string_view html_filename, std::string_view data_filename_name, std::string_view data_var_name);
    }

    template <typename Subset> void compare_sequences_generate_html(std::string_view html_filename, const subsets_to_compare_t<Subset>& data)
    {
        using namespace std::string_view_literals;

        const auto prefix{html_filename.substr(0, html_filename.size() - 5)};
        const auto data_filename{fmt::format("{}.data.js", prefix)};
        const auto data_var_name{fmt::format("compare_sequences_{}", ::string::replace(prefix, "/"sv, "_"sv, "-"sv, "_"sv))};
        acmacs::file::write(data_filename, fmt::format("const {} =\n{}", data_var_name, data.format_json(2)));

        std::string data_filename_name{data_filename};
        if (const auto pos = std::string_view{data_filename}.find_last_of('/'); pos != std::string_view::npos)
            data_filename_name.erase(0, pos + 1);

        detail::generate_html(html_filename, data_filename_name, data_var_name);

        // const auto templates_dir{fmt::format("{}/share/templates/seqdb-3", acmacs::acmacsd_root())};
        // acmacs::html::Generator html;
        // html.title("Compare sequences"sv);
        // html.add_css(acmacs::amino_acid_nucleotide_color_css());
        // html.add_css(static_cast<std::string>(acmacs::file::read(fmt::format("{}/compare-sequences.css", templates_dir))));
        // html.add_script_link(data_filename_name);
        // html.add_script(static_cast<std::string>(acmacs::file::read(fmt::format("{}/compare-sequences.js", templates_dir))));
        // html.add_script(fmt::format(R"(document.addEventListener("DOMContentLoaded", function() {{ compare_sequences({}); }});)", data_var_name));
        // html.add_to_body(static_cast<std::string>(acmacs::file::read(fmt::format("{}/compare-sequences.body.html", templates_dir))));
        // acmacs::file::write(html_filename, html.generate());
    }

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
