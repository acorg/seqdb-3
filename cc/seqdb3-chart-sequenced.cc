#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};
    option<str> aa_at_pos{*this, "aa-at-pos", desc{"comma separated list to filter: 162N,74R,!167X"}};

    argument<str> chart_name{*this, arg_name{"chart_name"}, mandatory};

};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        std::vector<acmacs::seqdb::subset::amino_acid_at_pos0_t> aa_at_pos;
        if (!opt.aa_at_pos->empty()) {
            const auto fields = acmacs::string::split(*opt.aa_at_pos, ",");
            aa_at_pos = ranges::to<decltype(aa_at_pos)>(fields | ranges::views::transform([](const auto& source) -> acmacs::seqdb::subset::amino_acid_at_pos0_t {
                            if (source.size() >= 2 && source.size() <= 4 && std::isdigit(source.front()) && std::isalpha(source.back()))
                                return {string::from_chars<size_t>(source.substr(0, source.size() - 1)) - 1, source.back(), true};
                            else if (source.size() >= 3 && source.size() <= 5 && source.front() == '!' && std::isdigit(source[1]) && std::isalpha(source.back()))
                                return {string::from_chars<size_t>(source.substr(1, source.size() - 2)) - 1, source.back(), false};
                            else
                                throw std::runtime_error{fmt::format("--aa-at: cannot parse entry: {}", source)};
            }));
        }

        auto chart = acmacs::chart::import_from_file(opt.chart_name);
        auto antigens = chart->antigens();
        auto sera = chart->sera();
        size_t matched_by_name = 0, matched_by_name_reassortant = 0;
        for (auto [ag_no, antigen] : acmacs::enumerate(*antigens)) {
            if (const auto subset = seqdb.select_by_name(antigen->name()).aa_at_pos(aa_at_pos); !subset.empty()) {
                ++matched_by_name;
                std::string matching_sera;
                if (const auto serum_indexes = sera->find_by_name(antigen->name()); !serum_indexes->empty())
                    matching_sera = fmt::format("    SR:{}", serum_indexes);
                fmt::print("{} {}{}\n", ag_no, antigen->full_name(), matching_sera);
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
