#include <memory>

#include "acmacs-base/to-json.hh"
#include "acmacs-base/read-file.hh"
#include "seqdb-3/create.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

struct filter_base
{
    virtual ~filter_base() = default;
    virtual bool good(const acmacs::seqdb::scan::sequence_t& seq) const = 0;
};

struct filter_all_aligned : public filter_base
{
    bool good(const acmacs::seqdb::scan::sequence_t& seq) const override { return seq.good(); }
};

struct filter_h1_h3_b_aligned : public filter_all_aligned
{
    bool good(const acmacs::seqdb::scan::sequence_t& seq) const override
    {
        // CDC sometimes puts H3N0 into gisaid and there is no HI match
        // We need to have them because they may be referenced by slave sequence entries
        return filter_all_aligned::good(seq) && (seq.type_subtype() == acmacs::virus::type_subtype_t{"B"} || seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H1N1)"} ||
                                                 seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H1)"} || seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H3N2)"} ||
                                                 seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H3)"});
    }
};

struct filter_whocc_aligned : public filter_h1_h3_b_aligned
{
    bool good(const acmacs::seqdb::scan::sequence_t& seq) const override
    {
        return filter_h1_h3_b_aligned::good(seq) && seq.lab_in({"CDC", "CRICK", "NIID", "VIDRL"}); // must be all uppercase because lab and lab_id are acmacs::uppercase
    }
};

static void generate(std::string_view filename, const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& sequences, const filter_base& filter);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::create(std::string_view prefix, std::vector<scan::fasta::scan_result_t>& sequences, create_dbs cdb)
{
    acmacs::seqdb::scan::fasta::sort_by_name(sequences);
#pragma omp parallel sections
    {
        if (cdb == create_dbs::all)
            generate(fmt::format("{}/seqdb-all.json.xz", prefix), sequences, filter_all_aligned{});
// #pragma omp section
//         if (cdb == create_dbs::all)
//             generate(fmt::format("{}/seqdb-h1-h3-b.json.xz", prefix), sequences, filter_h1_h3_b_aligned{});
#pragma omp section
        generate(fmt::format("{}/seqdb.json.xz", prefix), sequences, filter_h1_h3_b_aligned{}); // filter_whocc_aligned{});
    }

} // acmacs::seqdb::create

// ----------------------------------------------------------------------

void generate(std::string_view filename, const std::vector<acmacs::seqdb::scan::fasta::scan_result_t>& sequences, const filter_base& filter)
{
    to_json::array seqdb_data;
    to_json::object entry;
    to_json::array entry_seqs;
    acmacs::flat_set_t<std::string> dates;
    std::string previous;

    const auto flush = [&]() {
        if (!entry.empty()) {
            if (!dates.empty()) {
                if (dates.size() > 1 && std::any_of(std::begin(dates), std::end(dates), acmacs::seqdb::scan::not_empty_month_or_day))
                    dates.erase_if(acmacs::seqdb::scan::empty_month_or_day);
                dates.sort();
                entry << to_json::key_val("d", to_json::array(std::begin(dates), std::end(dates), to_json::json::compact_output::yes));
            }
            entry << to_json::key_val("s", entry_seqs);
            seqdb_data << std::move(entry);
        }
        dates.clear();
        entry_seqs = to_json::array{};
    };

    size_t num_sequences = 0;
    for (const auto& en : sequences) {
        const auto& seq = en.sequence;
        const std::string_view name = seq.name().get();
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
                if (!seq.annotations().empty())
                    entry_seq << to_json::key_val("A", seq.annotations());
                if (!seq.reassortant().empty())
                    entry_seq << to_json::key_val("r", to_json::array(*seq.reassortant(), to_json::json::compact_output::yes));
                if (!seq.passages().empty())
                    entry_seq << to_json::key_val("p", to_json::array(
                                                           seq.passages().begin(), seq.passages().end(), [](const auto& passage) { return *passage; }, to_json::json::compact_output::yes));
                if (en.reference) {
                    to_json::object reference = to_json::object{to_json::key_val("N", *en.reference->name), to_json::key_val("H", en.reference->hash)};
                    // if (!en.reference->annotations.empty())
                    //     reference << to_json::key_val("A", en.reference->annotations);
                    // if (!en.reference->reassortant.empty())
                    //     reference << to_json::key_val("r", *en.reference->reassortant);
                    // if (!en.reference->passage.empty())
                    //     reference << to_json::key_val("p", *en.reference->passage);
                    reference.make_compact();
                    entry_seq << to_json::key_val("R", std::move(reference));
                }
                else {
                    if (!seq.hash().empty())
                        entry_seq << to_json::key_val("H", seq.hash());
                    if (!seq.aa().empty())
                        entry_seq << to_json::key_val("a", seq.aa_format_not_aligned());
                    if (!seq.nuc().empty())
                        entry_seq << to_json::key_val("n", seq.nuc_format_not_aligned());
                    if (seq.shift_aa() != acmacs::seqdb::scan::sequence_t::shift_t{0})
                        entry_seq << to_json::key_val("s", -static_cast<ssize_t>(*seq.shift_aa()));
                    if (seq.shift_nuc() != acmacs::seqdb::scan::sequence_t::shift_t{0})
                        entry_seq << to_json::key_val("t", -static_cast<ssize_t>(*seq.shift_nuc()));
                    if (!seq.clades().empty())
                        entry_seq << to_json::key_val("c", to_json::array(
                                                               seq.clades().begin(), seq.clades().end(), [](const auto& clade) { return *clade; }, to_json::json::compact_output::yes));
                }
                if (!seq.hi_names().empty())
                    entry_seq << to_json::key_val("h", to_json::array(seq.hi_names().begin(), seq.hi_names().end(), to_json::json::compact_output::yes));
                if (!seq.lab_ids().empty()) {
                    to_json::object lab_ids;
                    for (const auto& [lab, ids] : seq.lab_ids()) {
                        lab_ids << to_json::key_val(
                            lab, to_json::array(
                                     std::begin(ids), std::end(ids), [](const auto& lab_id) { return *lab_id; }, to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    }
                    lab_ids.make_compact();
                    entry_seq << to_json::key_val("l", std::move(lab_ids));
                }
                // "g": "gene: HA|NA", // HA if omitted

                {
                    to_json::object gisaid_data;
                    if (!seq.isolate_id().empty())
                        gisaid_data << to_json::key_val("i",
                                                        to_json::array(seq.isolate_id().begin(), seq.isolate_id().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!seq.submitters().empty())
                        gisaid_data << to_json::key_val("S",
                                                        to_json::array(seq.submitters().begin(), seq.submitters().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!seq.sample_id_by_sample_provider().empty())
                        gisaid_data << to_json::key_val("s", to_json::array(seq.sample_id_by_sample_provider().begin(), seq.sample_id_by_sample_provider().end(), to_json::json::compact_output::yes,
                                                                            to_json::json::escape_double_quotes::yes));
                    if (!seq.gisaid_last_modified().empty())
                        gisaid_data << to_json::key_val(
                            "m", to_json::array(seq.gisaid_last_modified().begin(), seq.gisaid_last_modified().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!seq.originating_lab().empty())
                        gisaid_data << to_json::key_val(
                            "o", to_json::array(seq.originating_lab().begin(), seq.originating_lab().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!seq.gisaid_segment_number().empty())
                        gisaid_data << to_json::key_val(
                            "n", to_json::array(seq.gisaid_segment_number().begin(), seq.gisaid_segment_number().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!seq.gisaid_identifier().empty())
                        gisaid_data << to_json::key_val(
                            "t", to_json::array(seq.gisaid_identifier().begin(), seq.gisaid_identifier().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!seq.gisaid_dna_accession_no().empty())
                        gisaid_data << to_json::key_val("D", to_json::array(seq.gisaid_dna_accession_no().begin(), seq.gisaid_dna_accession_no().end(), to_json::json::compact_output::yes,
                                                                            to_json::json::escape_double_quotes::yes));
                    if (!seq.gisaid_dna_insdc().empty())
                        gisaid_data << to_json::key_val(
                            "d", to_json::array(seq.gisaid_dna_insdc().begin(), seq.gisaid_dna_insdc().end(), to_json::json::compact_output::yes, to_json::json::escape_double_quotes::yes));
                    if (!gisaid_data.empty())
                        entry_seq << to_json::key_val("G", std::move(gisaid_data));
                }

                entry_seqs << std::move(entry_seq);
            }

            dates.merge_from(seq.dates());
            ++num_sequences;
            // if (num_sequences > 5)
            //     break;
        }
    }
    flush();

    const auto js = to_json::object(to_json::key_val("_", "-*- js-indent-level: 1 -*-"), to_json::key_val("  version", "sequence-database-v3"), to_json::key_val("  date", date::current_date_time()),
                                    to_json::key_val("data", std::move(seqdb_data)));
    acmacs::file::write(filename, fmt::format("{:1}\n", js));
    fmt::print("INFO: {} sequences written to {}\n", num_sequences, filename);

} // generate

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
