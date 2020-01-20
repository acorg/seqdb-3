#include "acmacs-base/argv.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};
    option<str> clade{*this, "clade", desc{"report antigens/sera of that clade only"}};
    option<bool> indexes_only{*this, "indexes-only"};

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
        auto sera = chart->sera();
        auto layout = chart->number_of_projections() > 0 ? chart->projection(0)->layout() : std::shared_ptr<acmacs::Layout>{};

        const auto print = [&](bool is_ag, auto ag_no, auto antigen, const auto& clades) {
            if (!opt.indexes_only) {
                fmt::print("{} {:4d} {}{}  ::", is_ag ? "AG" : "SR", ag_no, antigen->full_name(), (is_ag && layout && !layout->point_has_coordinates(ag_no)) ? " <not-shown-on-map>" : "");
                for (const auto& clade : clades) {
                    // if (opt.gly || (clade != "GLY" && clade != "NO-GLY"))
                    fmt::print(" {}", clade);
                }
                fmt::print("\n");
            }
        };

        const auto show = [&](const auto& ag_sr, bool is_ag) {
            std::vector<size_t> indexes;
            for (auto [ag_no, antigen] : acmacs::enumerate(ag_sr)) {
                const auto clades = seqdb.clades_for_name(antigen->name());
                if (opt.clade.has_value()) {
                    if (std::find(std::begin(clades), std::end(clades), *opt.clade) != std::end(clades)) {
                        indexes.push_back(ag_no);
                        print(is_ag, ag_no, antigen, clades);
                    }
                }
                else {
                    print(is_ag, ag_no, antigen, clades);
                }
            }
            if (!indexes.empty()) {
                fmt::print("{} ({}) {}", is_ag ? "AG" : "SR", indexes.size(), indexes.front());
                for (auto ind = std::next(indexes.begin()); ind != indexes.end(); ++ind)
                    fmt::print(",{}", *ind);
                fmt::print("\n");
            }
        };

        show(*antigens, true);
        show(*sera, false);
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
