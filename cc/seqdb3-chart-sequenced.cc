#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};

    argument<str> chart_name{*this, arg_name{"chart_name"}, mandatory};

};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        auto chart = acmacs::chart::import_from_file(opt.chart_name);
        auto antigens = chart->antigens();
        size_t matched_by_name = 0, matched_by_name_reassortant = 0;
        for (auto [ag_no, antigen] : acmacs::enumerate(*antigens)) {
            if (const auto subset = seqdb.select_by_name(antigen->name()); !subset.empty()) {
                ++matched_by_name;
                fmt::print("{} {}\n", ag_no, antigen->full_name());
                bool this_reassortant_match = false;
                for (const auto& entry : subset) {
                    char match_symbol = ' ';
                    if (antigen->reassortant().empty() && entry.seq().reassortants.empty()) {
                        this_reassortant_match = true;
                        match_symbol = '+';
                    }
                    else if (entry.seq().has_reassortant(*antigen->reassortant())) {
                        this_reassortant_match = true;
                        match_symbol = 'R';
                    }
                    fmt::print("  {} {} {}\n", match_symbol, entry.full_name(), entry.seq().hi_names);
                }
                if (this_reassortant_match)
                    ++matched_by_name_reassortant;
            }
        }
        fmt::print(stderr, "INFO: matched_by_name {}\n", matched_by_name);
        fmt::print(stderr, "INFO: matched_by_name_reassortant {}\n", matched_by_name_reassortant);

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
