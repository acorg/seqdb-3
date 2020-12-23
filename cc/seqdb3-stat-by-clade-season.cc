#include <set>

#include "acmacs-base/argv.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/date.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db"};

    // select
    option<str>       subtype{*this, "flu", desc{"B, A(H1N1), H1, A(H3N2), H3"}};
    option<str>       lineage{*this, "lineage"};
    option<str>       start_date{*this, "start-date"};
    option<str>       end_date{*this, "end-date"};
    option<str>       clades{*this, "clades", desc{"comma separated list"}};

};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        const auto date_within = [&opt](std::string_view date) {
            return !date.empty() && (opt.start_date.empty() || date >= *opt.start_date) && (opt.end_date.empty() || date < *opt.end_date);
        };

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        auto ss = seqdb.all()
                .subtype(opt.subtype)
                .lineage(opt.lineage)
                .host("HUMAN");


        using data_key = std::pair<std::string, std::string>;
        std::map<data_key, acmacs::Counter<std::string>> data; // {continent, season} -> clade -> count
        std::set<std::string> all_clades;
        std::set<std::string> continents{"all"};

        for (const auto& ref : ss) {
            if (const auto date = date::from_string(ref.entry->date(), date::allow_incomplete::yes, date::throw_on_error::no); date::month_ok(date) && date_within(ref.entry->date()) && !ref.entry->continent.empty()) {
                std::string season;
                if (const auto month = date::get_month(date); month >= 4 && month < 10)
                    season = fmt::format("{}-04", date::get_year(date));
                else if (month >= 10)
                    season = fmt::format("{}-10", date::get_year(date));
                else
                    season = fmt::format("{}-10", date::get_year(date) - 1);
                const std::string continent{ref.entry->continent};
                continents.insert(continent);
                for (const auto& clade : ref.seq_with_sequence(seqdb).clades) {
                    all_clades.insert(std::string{clade});
                    data[data_key{continent, season}].count(clade);
                    data[data_key{"all", season}].count(clade);
                }
            }
        }

        fmt::print("Clades: {}\nContinents: {}\n", all_clades, continents);
        if (opt.clades.empty()) {
            for (const auto& [key, counter] : data) {
                fmt::print("{} {} {}\n", key.first, key.second, counter.report_sorted_max_first("  {first}: {second}"));
            }
        }
        else {
            const auto clades = acmacs::string::split(*opt.clades);
            for (const auto& continent : continents) {
                fmt::print("{}\n           ", continent);
                for (size_t i = 0; i < 2; ++i) {
                    for (const auto& clade : clades)
                        fmt::print("{:^12s}", clade);
                    fmt::print("       ");
                }
                fmt::print("\n");
                for (const auto& [key, counter] : data) {
                    if (key.first == continent) {
                        fmt::print("{}  ", key.second);
                        size_t sum{0};
                        for (const auto& clade : clades) {
                            const auto val = counter[std::string{clade}];
                            fmt::print("  {:7d}    ", val);
                            sum += val;
                        }
                        fmt::print("      ");
                        for (const auto& clade : clades) {
                            if (const auto val = counter[std::string{clade}]; val > 0)
                                fmt::print("  {:5.1f}%    ", static_cast<double>(val) / static_cast<double>(sum) * 100.0);
                            else
                                fmt::print("            ");
                        }
                        fmt::print("\n");
                    }
                }
                fmt::print("\n\n");
            }
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
