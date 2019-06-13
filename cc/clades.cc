#include "seqdb-3/clades.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

namespace local
{
    void detect_B_lineage(acmacs::seqdb::v3::sequence_t& sequence);
}

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::detect_lineages_clades(std::vector<fasta::scan_result_t>& sequences)
{
#pragma omp parallel for default(shared) schedule(static, 256)
    for (size_t e_no = 0; e_no < sequences.size(); ++e_no) {
        if (auto& entry = sequences[e_no]; fasta::is_aligned(entry)) {
            const auto& type_subtype = entry.sequence.type_subtype();
            if (*type_subtype == "B") {
                local::detect_B_lineage(entry.sequence);
            }
        }
    }
        // detect B lineage and VIC deletion mutants, adjust deletions
        // detect clades

} // acmacs::seqdb::v3::detect_lineages_clades

// ----------------------------------------------------------------------

void local::detect_B_lineage(acmacs::seqdb::v3::sequence_t& sequence)
{
    auto& deletions = sequence.deletions();
    if (deletions.empty()) {
        if (sequence.lineage() != "VICTORIA")
            fmt::print(stderr, "WARNING: lineage {} and no deletions\n", sequence.lineage());
    }
    else {
        // fmt::print(stderr, "{} :: {}\n", sequence.full_name(), format(sequence.deletions()));
    }

} // local::detect_B_lineage

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
