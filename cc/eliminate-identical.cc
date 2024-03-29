#include "seqdb-3/eliminate-identical.hh"
#include "seqdb-3/scan-fasta.hh"

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::eliminate_identical(std::vector<fasta::scan_result_t>& sequences)
{
    std::sort(std::begin(sequences), std::end(sequences), [](const auto& e1, const auto& e2) -> bool {
        if (e1.fasta.type_subtype == e2.fasta.type_subtype) {
            if (e1.sequence.nuc_shift() == e2.sequence.nuc_shift()) {
                if (const auto n1 = e1.sequence.nuc(), n2 = e2.sequence.nuc(); n1 == n2)
                    return (e1.sequence.good() ? 1 : 0) > (e2.sequence.good() ? 1 : 0); // good first, to avoid making bad (and excluded from seqdb) as a master
                else
                    return n1 < n2;
            }
            else
                return e1.sequence.nuc_shift() < e2.sequence.nuc_shift();
        }
        else
            return e1.fasta.type_subtype < e2.fasta.type_subtype;
    });

    size_t dups{0};
    for (auto seq = std::next(std::begin(sequences)), master = std::begin(sequences); seq != std::end(sequences); ++seq) {
        if (master->sequence.aligned() && seq->fasta.type_subtype == master->fasta.type_subtype && !seq->sequence.nuc().empty() && seq->sequence.nuc_shift() == master->sequence.nuc_shift() && seq->sequence.nuc() == master->sequence.nuc()) {
            if (!master->sequence.good() && seq->sequence.good())
                AD_WARNING("Master with issues ({}) for good {}", master->sequence.name(), seq->sequence.name());
            // fmt::print("DEBUG: identical: \"{} {}\"   \"{} {}\"\n", seq->fasta.name, seq->fasta.passage, master->fasta.name, master->fasta.passage);
            seq->reference = fasta::master_ref_t{master->sequence.name(), std::string{master->sequence.hash()}}; // std::string{master->sequence.annotations()}, master->sequence.reassortant(), master->sequence.passage()};
            ++dups;
        }
        else
            master = seq;
    }
    fmt::print(stderr, "INFO: entries with identical sequences: {}\n", dups);

} // acmacs::seqdb::v3::scan::eliminate_identical

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
