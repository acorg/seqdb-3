#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/read-file.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db"};
    option<str> html{*this, "html", desc{"generate html"}};

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
        if (opt.html)
            acmacs::file::write(opt.html, acmacs::seqdb::compare_report_html("", subset));
        else
            fmt::print("{}\n", acmacs::seqdb::compare_report_text(subset));

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
