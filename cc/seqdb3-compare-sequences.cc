#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/quicklook.hh"
#include "seqdb-3/compare.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>  db{*this, "db"};
    option<str>  html{*this, "html", desc{"generate html"}};
    option<bool> open{*this, "open", desc{"open html"}};

    argument<str_array> seq_ids{*this, arg_name{"seq-id or :T:<group-title>"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        if (opt.seq_ids->size() < 2)
            throw std::runtime_error("too few seq ids: nothing to compare");

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        acmacs::seqdb::v3::subsets_by_title_t subsets_by_title;

        std::string title;
        std::vector<std::string_view> seq_ids;
        for (const auto& seq_id : *opt.seq_ids) {
            if (seq_id.size() >= 3 && seq_id.substr(0, 3) == ":T:") {
                if (!seq_ids.empty()) {
                    subsets_by_title[title] = seqdb.find_by_seq_ids(seq_ids);
                    seq_ids.clear();
                }
                title = seq_id.substr(3);
                if (subsets_by_title.find(title) != subsets_by_title.end())
                    throw std::runtime_error(fmt::format("Group tag already exists: \"{}\"", title));
            }
            else if (seq_id.empty() || seq_id.front() == ':') {
                throw std::runtime_error(fmt::format("Unrecognized command entry: \"{}\"", seq_id));
            }
            else
                seq_ids.push_back(seq_id);
        }
        if (!seq_ids.empty())
            subsets_by_title[title] = seqdb.find_by_seq_ids(seq_ids);

        if (opt.html) {
            acmacs::file::write(opt.html, acmacs::seqdb::compare_report_html("", subsets_by_title));
            acmacs::open_or_quicklook(opt.open && opt.html != "-" && opt.html != "=", false, opt.html);
        }
        else
            fmt::print("{}\n", acmacs::seqdb::compare_report_text(subsets_by_title));

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
