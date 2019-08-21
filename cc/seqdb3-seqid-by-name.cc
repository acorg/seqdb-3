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

    argument<str_array> names{*this, arg_name{"name"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();
        for (const auto& name : *opt.names) {
            if (auto subset = seqdb.select_by_name(name); !subset.empty())
                subset.print("{seq_id}");
            else
                fmt::print(stderr, "ERROR: \"{}\" not found\n", name);
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
