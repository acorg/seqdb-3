#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/counter.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db"};

    argument<str_array> seqids{*this, arg_name{"seq-id or - to read them from stdin"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        acmacs::seqdb::subset subset;
        std::vector<std::string_view> seq_ids;
        if (opt.seqids->size() == 1 && opt.seqids->at(0) == "-")
            subset = seqdb.select_by_seq_id(acmacs::string::split(acmacs::file::read_stdin()));
        else
            subset = seqdb.select_by_seq_id(opt.seqids);

        std::vector<acmacs::CounterChar> aa_at_pos;
        for (const auto& ref : subset) {
            const auto aa = ref.aa_aligned(seqdb);
            if (aa_at_pos.size() < *aa.size())
                aa_at_pos.resize(*aa.size());
            for (acmacs::seqdb::pos0_t pos{0}; pos < aa.size(); ++pos)
                aa_at_pos[*pos].count(aa.at(pos));
        }

        for (size_t pos{0}; pos < aa_at_pos.size(); ++pos) {
            if (aa_at_pos[pos].size() < 2)
                continue;       // just one AA at this pos in all sequences
            fmt::print("{:3d} {}\n", pos + 1, aa_at_pos[pos].report_sorted_max_first("  {value}:{counter}"));
        }

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
