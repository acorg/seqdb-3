#include "acmacs-base/argv.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/factory-export.hh"
#include "acmacs-chart-2/chart-modify.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};
    option<str_array> verbose{*this, 'v', "verbose", desc{"comma separated list (or multiple switches) of enablers"}};

    argument<str_array> chart_name{*this, arg_name{"chart_name"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        acmacs::log::enable(opt.verbose);

        acmacs::seqdb::setup(opt.db);
        for (const auto& chart_name : *opt.chart_name) {
            acmacs::chart::ChartModify chart{acmacs::chart::import_from_file(chart_name)};
            const auto [matched_antigens, matched_sera] = acmacs::seqdb::get().populate(chart);
            AD_PRINT(FMT_STRING("{}\n  antigens: {:5d} (of {:5d})\n  sera:     {:5d} (of {:5d})"), chart_name, matched_antigens, chart.number_of_antigens(), matched_sera, chart.number_of_sera());
            acmacs::chart::export_factory(chart, chart_name, opt.program_name());
        }
        return 0;
    }
    catch (std::exception& err) {
        AD_ERROR("{}", err);
        return 1;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
