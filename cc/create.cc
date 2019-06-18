#include <memory>

#include "acmacs-base/to-json.hh"
#include "seqdb-3/create.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

namespace
{
    struct filter_base
    {
        virtual ~filter_base() = default;
        virtual bool good(const acmacs::seqdb::sequence_t& seq) const  = 0;
    };

    struct filter_all_aligned : public filter_base
    {
        bool good(const acmacs::seqdb::sequence_t& seq) const override { return seq.aligned(); }
    };

    struct filter_h1_h3_b_aligned : public filter_base
    {
        bool good(const acmacs::seqdb::sequence_t& seq) const override
        {
            return seq.aligned() && (seq.type_subtype() == acmacs::virus::type_subtype_t{"B"} || seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H1N1)"} || seq.type_subtype() == acmacs::virus::type_subtype_t{"A(H3N2)"});
        }
    };

    struct filter_whocc_aligned : public filter_base
    {
        bool good(const acmacs::seqdb::sequence_t& seq) const override { return seq.aligned(); }
    };

    inline std::unique_ptr<filter_base> make_filter(acmacs::seqdb::create_filter cfil)
    {
        switch (cfil) {
          case acmacs::seqdb::create_filter::all_aligned:
              return std::make_unique<filter_all_aligned>();
          case acmacs::seqdb::create_filter::h1_h3_b_aligned:
              return std::make_unique<filter_h1_h3_b_aligned>();
          case acmacs::seqdb::create_filter::whocc_aligned:
              return std::make_unique<filter_whocc_aligned>();
        }
    }
}

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::create(std::string_view filename, std::vector<fasta::scan_result_t>& sequences, create_filter cfilter)
{
    std::sort(std::begin(sequences), std::end(sequences), [](const auto& e1, const auto& e2) { return e1.sequence.name() < e2.sequence.name(); });

    to_json::array seqdb_data;
    to_json::object entry;
    to_json::array entry_seq;
    std::vector<std::string> dates;
    acmacs::virus::virus_name_t previous;

    const auto flush = [&]() {
        if (!entry.empty()) {
            if (!dates.empty()) {
                std::sort(std::begin(dates), std::end(dates));
                const auto end = std::unique(std::begin(dates), std::end(dates));
                entry << to_json::key_val("d", to_json::array(std::begin(dates), end));
            }
            entry << to_json::key_val("s", entry_seq);
            seqdb_data << std::move(entry);
        }
        dates.clear();
        entry_seq = to_json::array{};
    };

    auto filter = make_filter(cfilter);
    size_t num_sequences = 0;

    for (const auto& en : sequences) {
        const auto& seq = en.sequence;
        if (filter->good(seq)) {
            if (seq.name() == previous) {
            }
            else {
                flush();
                entry = to_json::object(to_json::key_val("N", *seq.name()),
                                        to_json::key_val("v", *seq.type_subtype()));
                if (!seq.lineage().empty())
                    entry << to_json::key_val("l", *seq.lineage());
                // "C": "continent", "c": "country",
            }
            entry_seq << to_json::object(to_json::key_val("a", seq.aa()));
            if (seq.date())
                dates.push_back(seq.date().display());
            ++num_sequences;
            if (num_sequences > 20)
                break;
        }
    }
    flush();

    auto js = to_json::object(to_json::key_val("_", "-*- js-indent-level: 1 -*-"),
                              to_json::key_val("  version", "sequence-database-v2"),
                              to_json::key_val("  date", current_date_time()),
                              to_json::key_val("data", std::move(seqdb_data)));
    fmt::print("{:1}\n", js);

} // acmacs::seqdb::create

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
