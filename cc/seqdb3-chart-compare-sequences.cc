#include "acmacs-base/argv.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/quicklook.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-chart-2/factory-import.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  db{*this, "db"};
    option<str>  html{*this, "html", desc{"generate html"}};
    option<str>  json{*this, 'j', "json", desc{"generate json"}};
    option<bool> nuc{*this, "nuc", desc{"compare nucleotide sequences"}};
    option<bool> open{*this, "open", desc{"open html"}};

    argument<str> chart{*this, arg_name{"chart"}, mandatory};
    argument<str_array> groups{*this, arg_name{"group-name,antigen-index,..."}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        acmacs::seqdb::setup(opt.db);
        auto chart = acmacs::chart::import_from_file(opt.chart);
        auto matched_seqdb = acmacs::seqdb::get().match(*chart->antigens(), chart->info()->virus_type(acmacs::chart::Info::Compute::Yes));

        const auto nuc_aa{opt.nuc ? acmacs::seqdb::compare::nuc : acmacs::seqdb::compare::aa};
        acmacs::seqdb::subsets_to_compare_t<acmacs::seqdb::subset_to_compare_t> subsets_to_compare{nuc_aa};
        for (const auto& group_desc : opt.groups) {
            const auto fields{acmacs::string::split(group_desc)};
            auto& subset = subsets_to_compare.subsets.emplace_back(fields[0]).subset;
            for (auto indp{std::next(std::begin(fields))}; indp != std::end(fields); ++indp) {
                if (const auto ind{acmacs::string::from_chars<size_t>(*indp)}; matched_seqdb[ind])
                    subset.append(matched_seqdb[ind]);
            }
        }
        subsets_to_compare.make_counters();

        if (opt.html) {
            acmacs::seqdb::compare_sequences_generate_html(opt.html, subsets_to_compare);
            acmacs::open_or_quicklook(opt.open && opt.html != "-" && opt.html != "=", false, opt.html);
        }
        if (opt.json)
            acmacs::file::write(opt.json, subsets_to_compare.format_json(2));

        fmt::print("{}\n\n{}\n\n", subsets_to_compare.format_seq_ids(0), subsets_to_compare.format_summary(0, 5, 0.2));

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------
