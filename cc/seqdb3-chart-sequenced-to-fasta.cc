#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db", dflt{""}};
    option<bool> nuc{*this, "nuc"};

    argument<str> chart_name{*this, arg_name{"chart.ace"}, mandatory};
    argument<str> fasta_name{*this, arg_name{"output.fasta"}, mandatory};

};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        auto chart = acmacs::chart::import_from_file(opt.chart_name);
        auto antigens = chart->antigens();
        auto subset = seqdb.match(*antigens, chart->info()->virus_type(acmacs::chart::Info::Compute::Yes));
        fmt::memory_buffer out;
        size_t matched{0};
        for (size_t ag_no{0}; ag_no < antigens->size(); ++ag_no) {
            if (!subset[ag_no].empty()) {
                fmt::format_to(out, ">{}\n{}\n", antigens->at(ag_no)->name_full(), *opt.nuc ? subset[ag_no].nuc_aligned(seqdb) : subset[ag_no].aa_aligned(seqdb));
                ++matched;
            }
        }
        AD_INFO("matched: {}", matched);
        if (matched)
            acmacs::file::write(opt.fasta_name, fmt::to_string(out));

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
