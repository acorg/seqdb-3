#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/acmacsd.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  name{*this, 'n', "name", dflt{""}};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        seqdb::Seqdb seqdb(acmacs::acmacsd_root() + "/data/seqdb.json.xz");

        seqdb::Seqdb::refs_t refs;
        if (!!opt.name)
            refs = seqdb.select_by_name(*opt.name);

        size_t selected = 0;
        for (const auto& [entry, seq_no] : refs) {
            fmt::print("{} <{}> [{}] {} {} {}\n", entry->name, entry->lineage, entry->dates, entry->seqs[seq_no].reassortants, entry->seqs[seq_no].passages, entry->seqs[seq_no].clades);
            ++selected;
        }
        fmt::print("INFO: selected: {}\n", selected);

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
