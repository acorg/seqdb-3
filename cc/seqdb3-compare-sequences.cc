#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/quicklook.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  db{*this, "db"};
    option<str>  html{*this, "html", desc{"generate html"}};
    option<bool> nuc{*this, "nuc", desc{"compare nucleotide sequences"}};
    option<bool> open{*this, "open", desc{"open html"}};

    argument<str_array> seq_ids{*this, arg_name{"seq-id or :T:<group-title>"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        if (opt.seq_ids->size() < 2)
            throw std::runtime_error("too few seq ids: nothing to compare");

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        const auto nuc_aa{opt.nuc ? acmacs::seqdb::compare::nuc : acmacs::seqdb::compare::aa};
        acmacs::seqdb::v3::subsets_to_compare_t subsets_to_compare{nuc_aa};

        acmacs::seqdb::subset* comparing_subset{nullptr};
        for (const auto& seq_id : *opt.seq_ids) {
            if (seq_id.size() >= 3 && seq_id.substr(0, 3) == ":T:") {
                comparing_subset = &subsets_to_compare.subsets.emplace_back(seq_id.substr(3)).subset;
            }
            else if (comparing_subset) {
                if (const auto selected = seqdb.select_by_seq_id(seq_id); !selected.empty())
                    comparing_subset->append(selected);
                else
                    AD_WARNING("No sequences found by seq_id: {}", seq_id);
            }
            else
                throw std::runtime_error{fmt::format("The first argument must be title (e.g. :T:name), found: \"{}\"", seq_id)};
        }
        subsets_to_compare.make_counters();

        if (opt.html) {
            acmacs::seqdb::compare_sequences_generate_html(opt.html, subsets_to_compare);
            acmacs::open_or_quicklook(opt.open && opt.html != "-" && opt.html != "=", false, opt.html);
        }
        else
            fmt::print("{}\n", subsets_to_compare.format_summary());

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
