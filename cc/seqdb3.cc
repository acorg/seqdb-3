#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/acmacsd.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  subtype{*this, "flu", dflt{""}};
    option<str>  host{*this, "host", dflt{""}};
    option<str>  lab{*this, "lab", dflt{""}};
    option<str>  lineage{*this, "lineage", dflt{""}};
    option<str>  start_date{*this, "start-date", dflt{""}};
    option<str>  end_date{*this, "end-date", dflt{""}};
    option<str>  continent{*this, "continent", dflt{""}};
    option<str>  country{*this, "country", dflt{""}};
    option<str>  clade{*this, "clade", dflt{""}};
    option<str>  db{*this, "db", dflt{""}};

//   aa at pos, not aa at pos
//   random N
//   recent N
//   with-hi-name
//   name matches regex (multiple regex possible, export all matching)
//   Base sequence to export together with other sequences (regex -> to select just one sequence)
//   HAMMING_DISTANCE_THRESHOLD - relative to base seq

    option<str>  name{*this, 'n', "name", dflt{""}};
};

// export
//   fasta (phylip?)
//   aligned
//   aa/nuc
//   with deletions inserted
//   wrapped
//   Truncate or extend with - all sequences to make them all of the same length, most common among original sequences.
//   sort by date or name
//   name format: {seq_id} {hi_name_or_seq_name_with_passage} {name} {date} {lab_id} {passage} {lab}
//   name encode
// select
//   subtype
//   lineage
//   date range
//   lab
//   continent/country
//   host (h1 swine)
//   clade
//   aa at pos, not aa at pos
//   random N
//   recent N
//   with-hi-name
//   name matches regex (multiple regex possible, export all matching)
//   Base sequence to export together with other sequences (regex -> to select just one sequence)
//   HAMMING_DISTANCE_THRESHOLD - relative to base seq

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        std::string db_filename{acmacs::acmacsd_root() + "/data/seqdb.json.xz"};
        if (!opt.db->empty())
            db_filename = opt.db;
        seqdb::Seqdb seqdb(db_filename);

        const auto init = [&] {
            if (!opt.name->empty())
                return seqdb.select_by_name(*opt.name);
            else
                return seqdb.all();
        };

        const auto subset = init();
        for (const auto& ref : subset)
            fmt::print("{} {}\n", ref.seq_id(), ref.entry->dates);

        // size_t selected = 0;
        // for (const auto& [entry, seq_no] : refs) {
        //     fmt::print("{} <{}> [{}] {} {} {}\n", entry->name, entry->lineage, entry->dates, entry->seqs[seq_no].reassortants, entry->seqs[seq_no].passages, entry->seqs[seq_no].clades);
        //     ++selected;
        // }
        // fmt::print("INFO: selected: {}\n", selected);

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
