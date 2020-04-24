#include "acmacs-base/argv.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/string.hh"
// #include "acmacs-base/date.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db"};

    // option<str>       subtype{*this, "flu", desc{"B, A(H1N1), H1, A(H3N2), H3"}};
    // option<str>       lineage{*this, "lineage"};
    // option<str>       start_date{*this, "start-date"};
    // option<str>       end_date{*this, "end-date"};
    // option<str>       clade{*this, "clade"};

};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        using namespace std::string_view_literals;
        using pp  = std::pair<std::string_view, std::string_view>;

        std::map<std::string_view, std::map<std::string_view, acmacs::CounterChar>> table; // clade -> year -> aa -> count
        std::set<char> all_aa;
        const std::array years{pp{"2016"sv, "2017"sv}, pp{"2017"sv, "2018"sv}, pp{"2018"sv, "2019"sv}, pp{"2019"sv, "2020"sv}, pp{"2020"sv, "2021"sv}};
        const std::array clades{"3C.2A"sv, "3C.2A1"sv, "3C.2A1A"sv, "3C.2A1B"sv, "3C.2A2"sv, "3C.2A3"sv, "3C.2A4"sv, "3C.3A"sv};

        for (const auto& [start, end] : years) {
            for (const auto clade : clades) {
                auto ss = seqdb.all()
                        .subtype("H3")
                        .host("HUMAN")
                        .dates(start, end)
                        .clade(seqdb, clade)
                        ;
                for (const auto& ref : ss) {
                    const auto aa = ref.aa_at_pos(seqdb, acmacs::seqdb::pos1_t{142});
                    table[clade][start].count(aa);
                    all_aa.insert(aa);
                }
            }
        }

        fmt::print("{}\n\n", all_aa);

        constexpr const int aa_sep_width{2};
        constexpr const int aa_width{4};
        constexpr const int aa_with_sep_width{aa_width + aa_sep_width * 2};
        constexpr const int year_sep_width{6};
        constexpr const int year_width{4};
        constexpr const int clade_width{8};
        const int year_sep_prefix_width{static_cast<int>(all_aa.size() * aa_with_sep_width + year_sep_width) / 2};

        fmt::print("{:{}s}", "", clade_width - 3);
        for (const auto& [year, end] : years)
            fmt::print("{:{}s}{}{:{}s}|", "", year_sep_prefix_width, year, "", year_sep_prefix_width);
        fmt::print("\n");

        fmt::print("{:{}s}", "", clade_width);
        for (const auto& y : years) {
            fmt::print("{:{}s}", "", aa_with_sep_width / 2);
            for (const auto aa : all_aa)
                fmt::print("{:{}c}", aa, aa_with_sep_width);
            fmt::print("{:{}s}|{:{}s}", "", year_sep_width / 2, "", year_sep_width / 2);
        }
        fmt::print("\n");

        for (const auto& clade : clades) {
            fmt::print("{:{}s}", ::string::lower(clade.substr(3)), clade_width);
            for (const auto& [start, end] : years) {
                for (const auto aa : all_aa) {
                    if (const auto count{table[clade][start][aa]}; count > 0)
                        fmt::print("{:{}s}{:{}d}{:{}s}", "", aa_sep_width, table[clade][start][aa], aa_width, "", aa_sep_width);
                    else
                        fmt::print("{:{}s}{:{}s}{:{}s}", "", aa_sep_width, "", aa_width, "", aa_sep_width);
                }
                fmt::print("{:{}s}|{:{}s}", "", year_sep_width + 1, "", year_sep_width / 2);
            }
            fmt::print("\n");
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
