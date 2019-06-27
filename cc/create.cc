#include <memory>

#include "acmacs-base/to-json.hh"
#include "acmacs-base/read-file.hh"
#include "seqdb-3/create.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

struct filter_base
{
    virtual ~filter_base() = default;
    virtual bool good(const acmacs::seqdb::sequence_t& seq) const = 0;
};

struct filter_all_aligned : public filter_base
{
    bool good(const acmacs::seqdb::sequence_t& seq) const override { return seq.aligned(); }
};

struct filter_h1_h3_b_aligned : public filter_all_aligned
{
    bool good(const acmacs::seqdb::sequence_t& seq) const override
    {
        return filter_all_aligned::good(seq) && (seq.type_subtype() == acmacs::virus::type_subtype_t{"B"} || seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H1N1)"} ||
                                                 seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H3N2)"});
    }
};

struct filter_whocc_aligned : public filter_h1_h3_b_aligned
{
    bool good(const acmacs::seqdb::sequence_t& seq) const override
    {
        return filter_h1_h3_b_aligned::good(seq) && (seq.lab() == "CDC" || seq.lab() == "Crick" || seq.lab() == "NIID" || seq.lab() == "VIDRL");
    }
};

static void generate(std::string_view filename, const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const filter_base& filter);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::create(std::string_view prefix, std::vector<fasta::scan_result_t>& sequences)
{
    acmacs::seqdb::fasta::sort_by_name(sequences);
#pragma omp parallel sections
    {
        generate(fmt::format("{}/seqdb-all.json.xz", prefix), sequences, filter_all_aligned{});
#pragma omp section
        generate(fmt::format("{}/seqdb-h1-h3-b.json.xz", prefix), sequences, filter_h1_h3_b_aligned{});
#pragma omp section
        generate(fmt::format("{}/seqdb.json.xz", prefix), sequences, filter_whocc_aligned{});
    }

} // acmacs::seqdb::create

// ----------------------------------------------------------------------

void generate(std::string_view filename, const std::vector<acmacs::seqdb::fasta::scan_result_t>& sequences, const filter_base& filter)
{
    to_json::array seqdb_data;
    to_json::object entry;
    to_json::array entry_seqs;
    std::vector<std::string> dates;
    std::string previous;

    const auto flush = [&]() {
        if (!entry.empty()) {
            if (!dates.empty()) {
                std::sort(std::begin(dates), std::end(dates));
                const auto end = std::unique(std::begin(dates), std::end(dates));
                entry << to_json::key_val("d", to_json::array(std::begin(dates), end, to_json::json::compact_output::yes));
            }
            entry << to_json::key_val("s", entry_seqs);
            seqdb_data << std::move(entry);
        }
        dates.clear();
        entry_seqs = to_json::array{};
    };

    const auto make_seq_name = [](const auto& seq) {
        if (seq.annotations().empty())
            return std::string{seq.name()};
        else
            return fmt::format("{} {}", *seq.name(), seq.annotations());
    };

    size_t num_sequences = 0;
    for (const auto& en : sequences) {
        const auto& seq = en.sequence;
        const auto name = make_seq_name(seq);
        if (filter.good(seq)) { //  && seq.type_subtype() == acmacs::virus::type_subtype_t{"B"}) {
            if (name != previous) {
                flush();
                entry = to_json::object(to_json::key_val("N", name), to_json::key_val("v", *seq.type_subtype()));
                if (!seq.lineage().empty())
                    entry << to_json::key_val("l", *seq.lineage());
                if (!seq.country().empty())
                    entry << to_json::key_val("c", seq.country());
                if (!seq.continent().empty())
                    entry << to_json::key_val("C", seq.continent());
                previous = name;
            }

            {
                to_json::object entry_seq;
                if (!seq.reassortant().empty())
                    entry_seq << to_json::key_val("r", to_json::array(*seq.reassortant(), to_json::json::compact_output::yes));
                if (!seq.passage().empty())
                    entry_seq << to_json::key_val("p", to_json::array(*seq.passage(), to_json::json::compact_output::yes));
                if (!seq.hi_names().empty())
                    entry_seq << to_json::key_val("h", to_json::array(seq.hi_names().begin(), seq.hi_names().end(), to_json::json::compact_output::yes));
                if (!seq.aa().empty())
                    entry_seq << to_json::key_val("a", seq.aa_format_not_aligned());
                if (!seq.nuc().empty())
                    entry_seq << to_json::key_val("n", seq.nuc_format_not_aligned());
                if (seq.shift_aa() != acmacs::seqdb::sequence_t::shift_t{0})
                    entry_seq << to_json::key_val("s", - static_cast<ssize_t>(*seq.shift_aa()));
                if (seq.shift_nuc() != acmacs::seqdb::sequence_t::shift_t{0})
                    entry_seq << to_json::key_val("t", - static_cast<ssize_t>(*seq.shift_nuc()));
                if (!seq.clades().empty())
                    entry_seq << to_json::key_val("c", to_json::array(seq.clades().begin(), seq.clades().end(), [](const auto& clade) { return *clade; }, to_json::json::compact_output::yes));
                if (!seq.lab().empty()) {
                    if (!seq.lab_id().empty())
                        entry_seq << to_json::key_val("l", to_json::object(to_json::key_val(seq.lab(), to_json::array{seq.lab_id()}), to_json::json::compact_output::yes));
                    else
                        entry_seq << to_json::key_val("l", to_json::object(to_json::key_val(seq.lab(), to_json::array{}), to_json::json::compact_output::yes));
                }
                // "g": "gene: HA|NA", // HA if omitted
                entry_seqs << std::move(entry_seq);
            }

            std::copy(std::begin(seq.dates()), std::end(seq.dates()), std::back_inserter(dates));
            ++num_sequences;
            // if (num_sequences > 5)
            //     break;
        }
    }
    flush();

    const auto js = to_json::object(to_json::key_val("_", "-*- js-indent-level: 1 -*-"), to_json::key_val("  version", "sequence-database-v2"), to_json::key_val("  date", current_date_time()),
                              to_json::key_val("data", std::move(seqdb_data)));
    acmacs::file::write(filename, fmt::format("{:1}\n", js));
    fmt::print("INFO: {} sequences written to {}\n", num_sequences, filename);

} // generate

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
