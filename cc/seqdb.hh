#pragma once

#include <map>
#include <mutex>

#include "acmacs-base/log.hh"
#include "acmacs-base/string-join.hh"
#include "acmacs-base/uppercase.hh"
#include "acmacs-base/flat-map.hh"
#include "acmacs-chart-2/chart-modify.hh"
#include "seqdb-3/aa-at-pos.hh"
#include "seqdb-3/seq-id.hh"
#include "seqdb-3/sequence-issues.hh"

// ----------------------------------------------------------------------

namespace acmacs::seqdb::inline v3
{
    struct SeqdbSeq;
    struct SeqdbEntry;
    struct ref;
    class subset;

    enum class even_if_already_popuplated { no, yes };

    using seq_id_index_t = map_with_duplicating_keys_t<seq_id_t, ref>; // duplicating seq_ids without hash present (for backward compatibility)
    using hi_name_index_t = map_with_unique_keys_t<std::string_view, ref>;
    using lab_id_index_t = map_with_duplicating_keys_t<std::string, ref>;
    using hash_index_t = map_with_duplicating_keys_t<std::string_view, ref>;

    class Seqdb
    {
      public:
        static const Seqdb& get();
        bool empty() const { return entries_.empty(); }

        subset all() const;
        subset select_by_seq_id(std::string_view seq_id) const;
        subset select_by_seq_id(const std::vector<std::string_view>& seq_ids) const;
        subset select_by_name(std::string_view name) const;
        subset select_by_name_hash(std::string_view name, std::string_view hash) const;
        subset select_by_name(const std::vector<std::string_view>& names) const;
        subset select_by_accession_number(const std::vector<std::string_view>& accession_numbers) const;
        subset select_by_regex(std::string_view re) const;
        subset select_by_lab_ids(const chart::LabIds& lab_ids) const;
        subset select_slaves() const;
        ref find_hi_name(std::string_view full_name) const;
        const seq_id_index_t& seq_id_index() const;
        const hi_name_index_t& hi_name_index() const;
        const lab_id_index_t& lab_id_index() const;
        const hash_index_t& hash_index() const;

        // returned subset contains elements for each antigen, i.e. it may contain empty ref's
        template <typename AgSr> subset match(const AgSr& antigens_sera, std::string_view aChartVirusType = {}) const;

        using aas_indexes_t = std::map<std::string, std::vector<size_t>>;
        aas_indexes_t aa_at_pos1_for_antigens(const acmacs::chart::Antigens& aAntigens, const std::vector<size_t>& aPositions1) const;

        using clade_t = std::string_view;
        using clades_t = std::vector<clade_t>;
        enum class clades_for_name_inclusive { no /* only common clades for matching sequences */, yes /* all possible clades */ };
        clades_t clades_for_name(std::string_view name, clades_for_name_inclusive inclusive = clades_for_name_inclusive::no) const;

        void populate(acmacs::chart::ChartModify& chart) const;

        // returns subset where each entry corresponds to the entry in seq_ids
        subset find_by_seq_ids(const std::vector<std::string_view>& seq_ids) const;

        // returns json with data for ace-view/2018 sequences_of_chart command
        std::string sequences_of_chart_for_ace_view_1(const acmacs::chart::Chart& chart) const;
        // returns sequences in the fasta format
        std::string sequences_of_chart_as_fasta(const acmacs::chart::Chart& chart) const;

        void find_slaves() const;

      private:
        std::string json_text_;
        std::vector<SeqdbEntry> entries_;
        mutable seq_id_index_t seq_id_index_;
        mutable hi_name_index_t hi_name_index_;
        mutable lab_id_index_t lab_id_index_;
        mutable hash_index_t hash_index_;
        mutable std::mutex index_access_; // acmacs-api is multi-threaded app
        mutable bool slaves_found_{false};

        Seqdb(std::string_view filename);
        // Seqdb(std::string&& source);

        void select_by_name(std::string_view name, subset& subs) const;

        // supports seq_ids before 2020-03-12, e.g. with _d1 etc. suffixes and without has suffix
        using seq_id_iter = typename seq_id_index_t::const_iterator;
        std::pair<seq_id_iter, seq_id_iter> find_seq_id(std::string_view seq_id) const;
    };

    void setup(std::string_view filename);
    inline const Seqdb& get() { return Seqdb::get(); }
    void populate(acmacs::chart::ChartModify& chart, even_if_already_popuplated eiap = even_if_already_popuplated::no);

    extern template subset Seqdb::match(const acmacs::chart::Antigens&, std::string_view) const;
    extern template subset Seqdb::match(const acmacs::chart::AntigensModify&, std::string_view) const;
    extern template subset Seqdb::match(const acmacs::chart::Sera&, std::string_view) const;
    extern template subset Seqdb::match(const acmacs::chart::SeraModify&, std::string_view) const;

    // ----------------------------------------------------------------------

    struct export_options
    {
        enum class format { fasta_aa, fasta_nuc };
        enum class aligned { no, yes };
        enum class most_common_length { no, yes };

        format e_format{format::fasta_nuc};
        size_t e_wrap_at{0};
        aligned e_aligned{aligned::yes};
        most_common_length e_most_common_length{most_common_length::no};
        std::string e_name_format{"{seq_id}"};
        size_t e_length{0};     // truncate/extend all sequences to this length
        size_t e_deletion_report_threshold{4}; // if sequence has this or more deletions, report the name.
        bool e_report_deletions_at_the_end{false};

        export_options& fasta(bool nucs)
        {
            e_format = nucs ? format::fasta_nuc : format::fasta_aa;
            return *this;
        }
        export_options& wrap(size_t wrap_at)
        {
            e_wrap_at = wrap_at;
            return *this;
        }
        export_options& no_wrap()
        {
            e_wrap_at = 0;
            return *this;
        }
        export_options& aligned(aligned a_aligned = aligned::yes)
        {
            e_aligned = a_aligned;
            return *this;
        }
        export_options& most_common_length(most_common_length a_most_common_length = most_common_length::yes)
        {
            e_most_common_length = a_most_common_length;
            return *this;
        }
        export_options& length(size_t len)
        {
            e_length = len;
            return *this;
        }
        export_options& name_format(std::string_view name_format)
        {
            e_name_format = name_format;
            return *this;
        }
        export_options& deletion_report_threshold(size_t deletion_report_threshold)
        {
            e_deletion_report_threshold = deletion_report_threshold;
            return *this;
        }
        export_options& deletion_report_threshold(std::string_view subtype)
        {
            if (subtype == "B" && e_deletion_report_threshold < 9) // do not report 3-del mutatants
                e_deletion_report_threshold = 9;
            return *this;
        }
        export_options& report_deletions_at_the_end(bool report)
        {
            e_report_deletions_at_the_end = report;
            return *this;
        }
    };

    // ----------------------------------------------------------------------

    struct SeqdbSeq
    {
        using lab_ids_t = std::vector<std::string_view>;
        using labs_t = std::vector<std::pair<std::string_view, lab_ids_t>>;

        struct gisaid_data_t
        {
            std::vector<std::string_view> isolate_ids; // gisaid accession numbers
            std::vector<std::string_view> sample_ids_by_sample_provider; // ncbi accession numbers
        };

        struct master_ref_t
        {
            std::string_view name;
            std::string_view hash;

            constexpr bool operator==(const master_ref_t& rhs) const { return name == rhs.name && hash == rhs.hash; }
            constexpr bool operator!=(const master_ref_t& rhs) const { return !operator==(rhs); }
            // std::string_view annotations;
            // std::string_view reassortant;
            // std::string_view passage;
            // constexpr bool operator==(const master_ref_t& rhs) const { return name == rhs.name && annotations == rhs.annotations && reassortant == rhs.reassortant && passage == rhs.passage; }
        };

        // sequence either contains nucs, amino_acids, clades or reference master sequence with the same nucs
        master_ref_t master; // for slave only
        sequence_with_alignment_ref_t amino_acids; // for master only
        sequence_with_alignment_ref_t nucs; // for master only
        std::string_view annotations;
        std::vector<std::string_view> reassortants;
        std::vector<std::string_view> passages;
        std::vector<std::string_view> clades; // for master only
        std::vector<std::string_view> hi_names;
        std::string_view hash;
        sequence::issues_t issues;
        labs_t lab_ids;
        gisaid_data_t gisaid;
        mutable std::unique_ptr<std::vector<ref>> slaves_; // for master only, list of slaves pointing to this master

        bool has_lab(std::string_view lab) const
        {
            return std::any_of(std::begin(lab_ids), std::end(lab_ids), [lab](const auto& en) { return en.first == lab; });
        }
        bool has_reassortant(std::string_view reassortant) const { return std::find(std::begin(reassortants), std::end(reassortants), reassortant) != std::end(reassortants); }
        bool matches(const amino_acid_at_pos1_eq_list_t& aa_at_pos1_eq) const { return acmacs::seqdb::matches(acmacs::seqdb::aligned(amino_acids), aa_at_pos1_eq); }
        bool matches(const amino_acid_at_pos1_list_t& aa_at_pos1) const { return acmacs::seqdb::matches(acmacs::seqdb::aligned(amino_acids), aa_at_pos1); }
        bool matches(const nucleotide_at_pos1_eq_list_t& nuc_at_pos1_eq) const { return acmacs::seqdb::matches(acmacs::seqdb::aligned(nucs), nuc_at_pos1_eq); }
        bool matches(const nucleotide_at_pos1_list_t& nuc_at_pos1) const { return acmacs::seqdb::matches(acmacs::seqdb::aligned(nucs), nuc_at_pos1); }

        constexpr bool matches_without_name(const master_ref_t& other_reference) const
        {
            return hash == other_reference.hash;
            // return annotations == other_reference.annotations &&
            //        ((other_reference.reassortant.empty() && reassortants.empty()) || std::find(std::begin(reassortants), std::end(reassortants), other_reference.reassortant) != std::end(reassortants)) &&
            //        ((other_reference.passage.empty() && passages.empty()) || (!passages.empty() && passages.front() == other_reference.passage)); // the first passage must match
        }

        // must not be used for slaves
        bool has_clade_master(std::string_view clade) const
        {
            if (!is_master())
                throw std::runtime_error(fmt::format("SeqdbSeq::has_clade_master is used for seq with the reference to {}, hi_names: {}", master.name, hi_names));
            return std::find(std::begin(clades), std::end(clades), clade) != std::end(clades);
        }

        constexpr sequence_aligned_ref_t aa_aligned_master(size_t length = std::string_view::npos) const { return acmacs::seqdb::aligned(amino_acids, length); }
        constexpr sequence_aligned_ref_t nuc_aligned_master(size_t length = std::string_view::npos) const { return acmacs::seqdb::aligned(nucs, length); }

        size_t aa_aligned_length_master() const { return acmacs::seqdb::aligned_length(amino_acids); }
        size_t nuc_aligned_length_master() const { return acmacs::seqdb::aligned_length(nucs); }

        char aa_at_pos_master(pos0_t pos0) const { return acmacs::seqdb::at_pos(amino_acids, pos0); }
        char aa_at_pos_master(pos1_t pos1) const { return acmacs::seqdb::at_pos(amino_acids, pos1); }

        std::string_view lab() const { return lab_ids.empty() ? std::string_view{} : lab_ids.front().first; }
        std::string_view lab_id() const { return (lab_ids.empty() || lab_ids.front().second.empty()) ? std::string_view{} : lab_ids.front().second.front(); }
        std::string_view passage() const { return passages.empty() ? std::string_view{} : passages.front(); }
        std::vector<std::string> designations(bool just_first = false) const;
        std::string designation() const { return designations(true).front(); }

        bool is_master() const { return master.name.empty(); }
        // bool is_slave() const { return !is_master(); }
        const SeqdbSeq& with_sequence(const Seqdb& seqdb) const { return is_master() ? *this : find_master(seqdb); }
        const SeqdbSeq& find_master(const Seqdb& seqdb) const;
        void add_slave(const ref& slave) const;
        const std::vector<ref>& slaves() const;
    };

    struct SeqdbEntry
    {
        std::string_view name;
        std::string_view continent;
        std::string_view country;
        std::vector<std::string_view> dates;
        std::string_view lineage;
        std::string_view virus_type;
        std::vector<SeqdbSeq> seqs;

        std::string host() const;
        bool date_within(std::string_view start, std::string_view end) const { return !dates.empty() && (start.empty() || dates.front() >= start) && (end.empty() || dates.front() < end); }
        std::string_view date() const;
        bool has_date(std::string_view date) const { return std::find(std::begin(dates), std::end(dates), date) != std::end(dates); }
        std::string location() const;
    };

    // ----------------------------------------------------------------------

    struct ref
    {
        const SeqdbEntry* entry;
        size_t seq_index;
        size_t group_no{0};         // for group_by_hamming_distance
        size_t hamming_distance{0}; // for group_by_hamming_distance and nuc_hamming_distance_to_base, printed using {hamming_distance}
        bool marked_for_removal{false}; // py-seqdb remove_nuc_duplicates_by_aligned_truncated

        ref() : entry{nullptr}, seq_index{static_cast<size_t>(-1)} {}
        ref(const SeqdbEntry* a_entry, size_t a_index) : entry{a_entry}, seq_index{a_index} {}
        ref(const SeqdbEntry& a_entry, size_t a_index) : entry{&a_entry}, seq_index{a_index} {}

        constexpr bool operator==(const ref& rhs) const { return entry == rhs.entry && seq_index == rhs.seq_index; }
        constexpr bool operator!=(const ref& rhs) const { return !operator==(rhs); }
        constexpr explicit operator bool() const { return entry != nullptr; }
        constexpr bool empty() const { return entry == nullptr; }

        const SeqdbSeq& seq() const { return entry->seqs[seq_index]; }
        const SeqdbSeq& seq_with_sequence(const Seqdb& seqdb) const { return seq().with_sequence(seqdb); }
        bool is_master() const { return seq().is_master(); }
        bool is_hi_matched() const { return !seq().hi_names.empty(); }
        seq_id_t seq_id() const;
        std::string full_name() const { return acmacs::string::join(acmacs::string::join_space, entry->name, acmacs::string::join(acmacs::string::join_space, seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()); }
        std::string full_name_with_date() const
        {
            return fmt::format("{} [{}]", acmacs::string::join(acmacs::string::join_space, entry->name, acmacs::string::join(acmacs::string::join_space, seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()),
                               entry->date());
        }
        std::string hi_name_or_full_name() const
        {
            if (seq().hi_names.empty())
                return full_name();
            else
                return std::string{seq().hi_names.front()};
        }
        bool has_issues() const { return !seq().issues.none(); }
        bool has_lab(std::string_view lab) const { return seq().has_lab(lab); }
        bool has_clade(const Seqdb& seqdb, std::string_view clade) const { return seq_with_sequence(seqdb).has_clade_master(clade); }
        bool has_hi_names() const { return !seq().hi_names.empty(); }
        bool matches(const amino_acid_at_pos1_eq_list_t& aa_at_pos1) const { return seq().matches(aa_at_pos1); }
        bool matches(const amino_acid_at_pos1_list_t& aa_at_pos1) const { return seq().matches(aa_at_pos1); }
        constexpr bool matches(const SeqdbSeq::master_ref_t& master) const { return entry->name == master.name && seq().matches_without_name(master); }

        sequence_aligned_ref_t aa_aligned(const Seqdb& seqdb, size_t length = std::string_view::npos) const { return seq_with_sequence(seqdb).aa_aligned_master(length); }
        sequence_aligned_ref_t nuc_aligned(const Seqdb& seqdb, size_t length = std::string_view::npos) const { return seq_with_sequence(seqdb).nuc_aligned_master(length); }
        size_t aa_aligned_length(const Seqdb& seqdb) const { return seq_with_sequence(seqdb).aa_aligned_length_master(); }
        size_t nuc_aligned_length(const Seqdb& seqdb) const { return seq_with_sequence(seqdb).nuc_aligned_length_master(); }
        char aa_at_pos(const Seqdb& seqdb, pos0_t pos0) const { return seq_with_sequence(seqdb).aa_at_pos_master(pos0); }
        char aa_at_pos(const Seqdb& seqdb, pos1_t pos1) const { return seq_with_sequence(seqdb).aa_at_pos_master(pos1); }
    };

    class subset
    {
      public:
        using refs_t = std::vector<ref>;
        using iterator = typename refs_t::iterator;
        using const_iterator = typename refs_t::const_iterator;

        enum class sorting { none, name_asc, name_desc, date_asc, date_desc };
        enum class master_only { no, yes };

        subset() = default;
        subset(iterator first, iterator last) : refs_{first, last} {} // copy of the part of another subset
        subset(const_iterator first, const_iterator last) : refs_{first, last} {} // copy of the part of another subset

        auto empty() const { return refs_.empty(); }
        auto size() const { return refs_.size(); }
        const_iterator begin() const { return refs_.begin(); }
        const_iterator end() const { return refs_.end(); }
        iterator begin() { return refs_.begin(); }
        iterator end() { return refs_.end(); }
        const auto& operator[](size_t index) const { return refs_.at(index); }
        const auto& front() const { return refs_.front(); }

        subset& multiple_dates(bool do_filter = true);
        subset& subtype(const acmacs::uppercase& virus_type);
        subset& lineage(const acmacs::uppercase& lineage);
        subset& lab(const acmacs::uppercase& lab);
        subset& whocc_lab(bool do_filter = true);
        subset& host(const acmacs::uppercase& host);
        subset& dates(std::string_view start, std::string_view end);
        subset& continent(const acmacs::uppercase& continent);
        subset& country(const acmacs::uppercase& country);
        subset& with_issues(bool keep_with_issues);
        subset& clade(const Seqdb& seqdb, const acmacs::uppercase& clade);
        subset& recent(size_t recent, master_only master);
        subset& recent_matched(const std::vector<size_t>& recent_matched, master_only master);
        subset& random(size_t random);
        subset& subset_every_month(double fraction);
        subset& group_by_hamming_distance(const Seqdb& seqdb, size_t dist_threshold, size_t output_size);  // Eu's algorithm 2019-07-23
        subset& subset_by_hamming_distance_random(const Seqdb& seqdb, bool do_subset, size_t output_size); // davipatti algorithm 2019-07-23
        subset& remove_nuc_duplicates(bool do_remove, bool keep_hi_matched);
        subset& remove_empty(const Seqdb& seqdb, bool nuc);
        subset& remove_with_front_back_deletions(const Seqdb& seqdb, bool remove, size_t nuc_length);
        subset& remove_with_deletions(const Seqdb& seqdb, bool remove, size_t threshold); // remove if number of deletions >= threshold
        subset& remove_marked(); // ref::marked_for_removal
        subset& with_hi_name(bool with_hi_name);
        subset& aa_at_pos(const Seqdb& seqdb, const amino_acid_at_pos1_eq_list_t& aa_at_pos);
        subset& nuc_at_pos(const Seqdb& seqdb, const nucleotide_at_pos1_eq_list_t& nuc_at_pos);
        subset& min_aa_length(const Seqdb& seqdb, size_t length);
        subset& min_nuc_length(const Seqdb& seqdb, size_t length);
        subset& names_matching_regex(const std::vector<std::string_view>& regex_list);
        subset& names_matching_regex(std::string_view re) { return names_matching_regex(std::vector<std::string_view>{re}); }
        subset& exclude(const std::vector<std::string_view>& seq_ids);
        subset& keep_master_only();
        subset& prepend(std::string_view seq_id, const Seqdb& seqdb);
        subset& prepend(const std::vector<std::string_view>& seq_ids, const Seqdb& seqdb);
        // subset& prepend_single_matching(std::string_view re, const Seqdb& seqdb);
        subset& nuc_hamming_distance_mean(size_t threshold, size_t size_threshold = 1000);
        subset& nuc_hamming_distance_to(size_t threshold, std::string_view seq_id);
        subset& nuc_hamming_distance_to_base(size_t threshold, bool do_filter = true);
        subset& sort(sorting srt);
        std::pair<size_t, std::string> export_sequences(const Seqdb& seqdb, const export_options& options) const; // returns {num_sequences, fasta_data}
        subset& export_sequences(std::string_view filename, const Seqdb& seqdb, const export_options& options) const;
        subset& export_json_sequences(std::string_view filename, const Seqdb& seqdb, const export_options& options);
        subset& print(const Seqdb& seqdb, std::string_view name_format, std::string_view header = {}, bool do_print = true) const;
        subset& report_stat(const Seqdb& seqdb, bool do_report = true);
        subset& report_stat_month_region(bool do_report = true);
        subset& report_aa_at(const Seqdb& seqdb, const pos1_list_t& pos1_list);
        subset& report_hamming_distance(bool do_report);
        subset& report_hamming_bins(const Seqdb& seqdb, size_t bin_size);

        subset& append(const ref& seq)
        {
            if (std::find_if(std::begin(refs_), std::end(refs_), [&seq](const auto& en) { return en.entry == seq.entry && en.seq_index == seq.seq_index; }) == std::end(refs_))
                refs_.push_back(seq);
            return *this;
        }

        subset& append(const subset& another) {
            for (const auto& en : another)
                append(en);
            return *this;
        }

        // returns new subset, this subset is not modified
        enum class matched_only { no, yes };
        subset filter_by_indexes(const acmacs::chart::PointIndexList& indexes, enum matched_only matched_only = matched_only::yes) const;

        void sort_by_nuc_aligned_truncated(const Seqdb& seqdb, size_t truncate_at)
        {
            std::sort(std::begin(refs_), std::end(refs_), [&seqdb, truncate_at](const auto& e1, const auto& e2) { return e1.nuc_aligned(seqdb, truncate_at) < e2.nuc_aligned(seqdb, truncate_at); });
        }

      private:
        refs_t refs_;
        using ref_indexes = std::vector<size_t>;

        subset(size_t size) : refs_(size) {}

        void sort_by_name_asc()
        {
            std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.seq_id() < e2.seq_id(); });
        }
        void sort_by_name_desc()
        {
            std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.seq_id() > e2.seq_id(); });
        }
        void sort_by_date_recent_first()
        {
            std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.entry->date() > e2.entry->date(); });
        }
        void sort_by_date_oldest_first()
        {
            std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.entry->date() < e2.entry->date(); });
        }
        void sort_by_hamming_distance()
        {
            std::sort(std::begin(refs_), std::end(refs_), [](const auto& e1, const auto& e2) { return e1.hamming_distance < e2.hamming_distance; });
        }

        refs_t::iterator most_recent_with_hi_name();
        void remove(ref_indexes& to_remove);
        void keep(ref_indexes& to_keep);

        void resize(size_t size) { refs_.resize(size); }

        struct collected_entry_t
        {
            std::string seq_id;
            std::string sequence;
        };
        using collected_t = std::vector<collected_entry_t>;

        collected_t export_collect(const Seqdb& seqdb, const export_options& options) const;
        std::string export_fasta(const collected_t& entries, const export_options& options) const;
        std::string export_json(const collected_t& entries, const export_options& options) const;
        std::string make_name(const Seqdb& seqdb, std::string_view name_format, const ref& entry) const;

        friend class Seqdb;
    };

    // ----------------------------------------------------------------------

    void remove_nuc_duplicates(subset::refs_t& refs, bool keep_hi_matched = false);

} // namespace acmacs::seqdb::inlinev3

// ----------------------------------------------------------------------

template <> struct fmt::formatter<acmacs::seqdb::v3::ref> : fmt::formatter<acmacs::fmt_helper::default_formatter> {
    template <typename FormatCtx> auto format(const acmacs::seqdb::v3::ref& rf, FormatCtx& ctx)
    {
        return format_to(ctx.out(), "{}", rf.seq_id());
    }
};

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
