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

    argument<str> chart_name{*this, arg_name{"chart.ace"}, mandatory};
    argument<str> fasta_name{*this, arg_name{"output.fasta"}, mandatory};

};

// static acmacs::seqdb::v3::ref* match_by_lab_id(const acmacs::seqdb::Seqdb& seqdb, const acmacs::chart::Antigen& antigen);

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        auto chart = acmacs::chart::import_from_file(opt.chart_name);
        auto antigens = chart->antigens();
        auto subset = seqdb.match(*antigens, chart->info()->virus_type(acmacs::chart::Info::Compute::Yes));
        size_t matched { 0 };
        for (const auto& ref : subset) {
            if (!ref.empty())
                ++matched;
        }
        AD_INFO("matched: {}", matched);

        // auto sera = chart->sera();
        // size_t matched_by_name = 0, matched_by_name_reassortant = 0;
        // for (auto [ag_no, antigen] : acmacs::enumerate(*antigens)) {
        //     auto* ref = match_by_lab_id(seqdb, antigen);

        //     if (const auto subset = seqdb.select_by_name(antigen->name()); !subset.empty()) {
        //         ++matched_by_name;
        //         std::string matching_sera;
        //         if (const auto serum_indexes = sera->find_by_name(antigen->name()); !serum_indexes->empty())
        //             matching_sera = fmt::format("    SR:{}", serum_indexes);
        //         fmt::print("{} {}{}\n", ag_no, antigen->format("{name_full}"), matching_sera);
        //         bool this_reassortant_match = false;
        //         for (const auto& entry : subset) {
        //             char match_symbol = ' ';
        //             if (antigen->reassortant().empty() && entry.seq().reassortants.empty()) {
        //                 this_reassortant_match = true;
        //                 match_symbol = '+';
        //             }
        //             else if (entry.seq().has_reassortant(*antigen->reassortant())) {
        //                 this_reassortant_match = true;
        //                 match_symbol = 'R';
        //             }
        //             fmt::print("  {} {} {}\n", match_symbol, entry.full_name(), entry.seq().hi_names);
        //         }
        //         if (this_reassortant_match)
        //             ++matched_by_name_reassortant;
        //     }
        // }
        // fmt::print(stderr, "INFO: matched_by_name {}\n", matched_by_name);
        // fmt::print(stderr, "INFO: matched_by_name_reassortant {}\n", matched_by_name_reassortant);

        return 0;
    }
    catch (std::exception& err) {
        AD_ERROR("{}", err);
        return 1;
    }
}

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::ref* match_by_lab_id(const acmacs::seqdb::Seqdb& seqdb, const acmacs::chart::Antigen& antigen)
// {
// seqdb.select_by_name
// } // match_by_lab_id

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
