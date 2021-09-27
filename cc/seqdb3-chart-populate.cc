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
    option<bool> no_export{*this, 'n', "no-export"};
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
            AD_PRINT("{}", chart_name);
            const auto [matched_antigens, matched_sera] = acmacs::seqdb::get().populate(chart);

            AD_PRINT("AG matched: {} (of {})", matched_antigens.size(), chart.number_of_antigens());
            if (matched_antigens.size() < chart.number_of_antigens()) {
                AD_PRINT("AG NOT matched: {}", chart.number_of_antigens() - matched_antigens.size());
                for (const auto ag_no : range_from_0_to(chart.number_of_antigens())) {
                    if (!matched_antigens.contains(ag_no))
                        AD_PRINT("  {:5d} {}", ag_no, chart.antigens()->at(ag_no)->name_full());
                }
            }

            AD_PRINT("SR matched: {} (of {})", matched_sera.size(), chart.number_of_sera());
            if (matched_sera.size() < chart.number_of_sera()) {
                AD_PRINT("SR NOT matched: {}", chart.number_of_sera() - matched_sera.size());
                for (const auto sr_no : range_from_0_to(chart.number_of_sera())) {
                    if (!matched_sera.contains(sr_no))
                        AD_PRINT("  {:5d} {}", sr_no, chart.sera()->at(sr_no)->name_full());
                }
            }

            // for (const auto ag_no : matched_antigens)
            //     AD_PRINT("    {:5d} {} {}", ag_no, chart.antigens()->at(ag_no)->name_full(), chart.antigens()->at(ag_no)->clades());
            // AD_PRINT("  sera: {:5d} (of {:5d})", matched_sera.size(), chart.number_of_sera());
            // for (const auto sr_no : matched_sera)
            //     AD_PRINT("    {:5d} {} {}", sr_no, chart.sera()->at(sr_no)->name_full(), chart.sera()->at(sr_no)->clades());

            if (!opt.no_export)
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
