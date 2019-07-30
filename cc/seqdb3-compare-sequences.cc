#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/enumerate.hh"
#include "acmacs-base/range.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

static void report_text(const acmacs::seqdb::subset& subset);

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};

    argument<str_array> seq_ids{*this, arg_name{"seq-id"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        if (opt.seq_ids->size() < 2)
            throw std::runtime_error("too few seq ids: nothing to compare");

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();
        const auto subset = seqdb.find_by_seq_ids(*opt.seq_ids);
        report_text(subset);

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------

void report_text(const acmacs::seqdb::subset& subset)
{
    for (auto [no, ref] : acmacs::enumerate(subset))
        fmt::print("{:2d} {}\n", no, ref.seq_id());
    fmt::print("\n   ");
    for (auto no : acmacs::range(subset.size()))
        fmt::print(" {:>2d}", no);
    fmt::print("\n");

    std::vector<std::string_view> seqs(subset.size());
    std::transform(std::begin(subset), std::end(subset), std::begin(seqs), [](const auto& ref) { return ref.seq().aa_aligned(); });

    // std::for_each(std::next(std::begin(seqs)), std::end(seqs), [](std::string_view seq) { fmt::print("{}\n\n", seq.size()); });

    for (size_t pos = 0; pos < seqs[0].size(); ++pos) {
        const auto aa = seqs[0][pos];
        if (!std::all_of(std::next(std::begin(seqs)), std::end(seqs), [aa, pos](std::string_view seq) {
                return seq.size() <= pos || seq[pos] == aa;
            })) {
            fmt::print("{:3d} {:>2c}", pos + 1, aa);
            std::for_each(std::next(std::begin(seqs)), std::end(seqs), [aa, pos](std::string_view seq) {
                char seq_aa = '-';
                if (seq.size() <= pos)
                    seq_aa = '-';
                else if (seq[pos] == aa)
                    seq_aa = '.';
                else
                    seq_aa = seq[pos];
                fmt::print(" {:>2c}", seq_aa);
            });
            fmt::print("\n");
        }
    }

} // report_text

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
