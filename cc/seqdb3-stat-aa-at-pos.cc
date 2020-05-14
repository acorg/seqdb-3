#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/counter.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db"};
    option<double> second_counter_threshold{*this, "threshold", dflt{0.05}, "fraction of the second frequent AA at pos must be bigger that this value to report"};

    argument<str_array> seqids{*this, arg_name{"seq-id or - to read them from stdin"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);
        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        acmacs::seqdb::subset subset;
        std::vector<std::string_view> seq_ids;
        if (opt.seqids->size() == 1 && opt.seqids->at(0) == "-")
            subset = seqdb.select_by_seq_id(acmacs::string::split(acmacs::file::read_stdin()));
        else
            subset = seqdb.select_by_seq_id(opt.seqids);

        std::vector<acmacs::CounterChar> aa_at_pos;
        for (const auto& ref : subset) {
            const auto aa = ref.aa_aligned(seqdb);
            if (aa_at_pos.size() < *aa.size())
                aa_at_pos.resize(*aa.size());
            for (acmacs::seqdb::pos0_t pos{0}; pos < aa.size(); ++pos)
                aa_at_pos[*pos].count(aa.at(pos));
        }

        const auto min_second_to_report = static_cast<size_t>(static_cast<double>(subset.size()) * *opt.second_counter_threshold);
        // AD_DEBUG("subset:{} threshold:{} min-second:{}", subset.size(), static_cast<double>(subset.size()) * *opt.second_counter_threshold, min_second_to_report);
        for (size_t pos{0}; pos < aa_at_pos.size(); ++pos) {
            const auto data = aa_at_pos[pos].pairs(acmacs::CounterChar::sorted::yes);
            if (data.size() > 1 && data[1].first != 'X' && data[1].second > min_second_to_report)
                fmt::print("{:3d} {}\n", pos + 1, aa_at_pos[pos].report_sorted_max_first("  {value}:{counter}"));
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
