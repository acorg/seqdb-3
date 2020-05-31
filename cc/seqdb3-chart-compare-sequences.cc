#include "acmacs-base/argv.hh"
#include "acmacs-base/string-split.hh"
// #include "acmacs-base/read-file.hh"
#include "acmacs-base/quicklook.hh"
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

        acmacs::seqdb::v3::subsets_by_title_t subsets_by_title;
        for (const auto& group_desc : opt.groups) {
            const auto fields{acmacs::string::split(group_desc)};
            if (subsets_by_title.find(fields[0]) != subsets_by_title.end())
                    throw std::runtime_error(fmt::format("Group tag already exists: \"{}\"", fields[0]));
            auto& subset = subsets_by_title.emplace(fields[0], acmacs::seqdb::subset{}).first->second;
            for (auto indp{std::next(std::begin(fields))}; indp != std::end(fields); ++indp) {
                if (const auto ind{acmacs::string::from_chars<size_t>(*indp)}; matched_seqdb[ind])
                    subset.append(matched_seqdb[ind]);
            }
        }
        fmt::print("{}\n", acmacs::seqdb::compare_report_text(subsets_by_title, opt.nuc ? acmacs::seqdb::compare::nuc : acmacs::seqdb::compare::aa));

        // std::string title;
        // std::vector<std::string_view> seq_ids;
        // for (const auto& seq_id : *opt.seq_ids) {
        //     if (seq_id.size() >= 3 && seq_id.substr(0, 3) == ":T:") {
        //         if (!seq_ids.empty()) {
        //             subsets_by_title[title] = seqdb.find_by_seq_ids(seq_ids);
        //             seq_ids.clear();
        //         }
        //         title = seq_id.substr(3);
        //         if (subsets_by_title.find(title) != subsets_by_title.end())
        //             throw std::runtime_error(fmt::format("Group tag already exists: \"{}\"", title));
        //     }
        //     else if (seq_id.empty() || seq_id.front() == ':') {
        //         throw std::runtime_error(fmt::format("Unrecognized command entry: \"{}\"", seq_id));
        //     }
        //     else
        //         seq_ids.push_back(seq_id);
        // }
        // if (!seq_ids.empty()) {
        //     subsets_by_title[title] = seqdb.find_by_seq_ids(seq_ids);
        //     bool not_found_present = false;
        //     for (size_t i = 0; i < seq_ids.size(); ++i) {
        //         if (!subsets_by_title[title][i]) {
        //             AD_ERROR("{} not found", seq_ids[i]);
        //             not_found_present = true;
        //         }
        //     }
        //     if (not_found_present)
        //         throw std::runtime_error{"few sequences not found"};
        // }

        // if (opt.html) {
        //     acmacs::file::write(opt.html, acmacs::seqdb::compare_report_html("", subsets_by_title, opt.nuc ? acmacs::seqdb::compare::nuc : acmacs::seqdb::compare::aa));
        //     acmacs::open_or_quicklook(opt.open && opt.html != "-" && opt.html != "=", false, opt.html);
        // }
        // else
        //     fmt::print("{}\n", acmacs::seqdb::compare_report_text(subsets_by_title, opt.nuc ? acmacs::seqdb::compare::nuc : acmacs::seqdb::compare::aa));

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
