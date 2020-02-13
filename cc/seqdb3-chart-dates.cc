#include "acmacs-base/argv.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/factory-export.hh"
#include "acmacs-chart-2/chart-modify.hh"
#include "seqdb-3/seqdb.hh"

// puts dates from seqdb into chart, if absent in chart, report conflicts

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};
    option<bool> verbose{*this, 'v', "verbose", desc{"\n\nProgram puts dates from seqdb into chart, if absent in chart, reports conflicts.\n"}};

    argument<str> source_chart{*this, arg_name{"source-chart"}, mandatory};
    argument<str> output_chart{*this, arg_name{"output-chart"}};

};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();
        acmacs::chart::ChartModify chart{acmacs::chart::import_from_file(opt.source_chart)};
        auto antigens = chart.antigens_modify();
        const auto subset = seqdb.match(*antigens);
        for (const auto ag_no : acmacs::range(antigens->size())) {
            auto& antigen = antigens->at(ag_no);
            const auto& ref = subset[ag_no];
            const auto sequenced = !ref.empty();
            if (opt.verbose) {
                fmt::print("AG {:4d} {} [{}]", ag_no, antigen.full_name(), antigen.date());
                if (sequenced)
                    fmt::print("  {} [{}]\n", ref.seq_id(), ref.entry->dates);
                else
                    fmt::print(" *not sequenced*\n");
            }
            if (sequenced) {
                if (antigen.date().empty()) {
                    antigen.date(ref.entry->date());
                }
                else if (!ref.entry->has_date(antigen.date())) {
                    fmt::print(stderr, "WARNING: AG {} {} has different dates: table: {} seqdb: {}\n", ag_no, antigen.full_name(), antigen.date(), ref.entry->dates);
                }
            }
        }

        if (opt.output_chart)
            acmacs::chart::export_factory(chart, opt.output_chart, opt.program_name());
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
