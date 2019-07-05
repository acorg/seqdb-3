#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/date.hh"
#include "acmacs-base/acmacsd.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/range-v3.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str>       db{*this, "db", dflt{""}};

    option<str>       subtype{*this, "flu", dflt{""}, desc{"B, A(H1N1), H1, A(H3N2), H3"}};
    option<str>       host{*this, "host", dflt{""}};
    option<str>       lab{*this, "lab", dflt{""}};
    option<str>       lineage{*this, "lineage", dflt{""}};
    option<str>       start_date{*this, "start-date", dflt{""}};
    option<str>       end_date{*this, "end-date", dflt{""}};
    option<str>       continent{*this, "continent", dflt{""}};
    option<str>       country{*this, "country", dflt{""}};
    option<str>       clade{*this, "clade", dflt{""}};
    option<str>       aa_at_pos{*this, "aa-at-pos", dflt{""}, desc{"comma separated list: 162N,74R,!167X"}};
    option<size_t>    recent{*this, "recent", dflt{0UL}};
    option<size_t>    random{*this, "random", dflt{0UL}};
    option<bool>      with_hi_name{*this, "with-hi-name"};
    option<str_array> name_regex{*this, "re", desc{"filter names by regex, multiple regex possible, all matching listed"}};

//   Base sequence to export together with other sequences (regex -> to select just one sequence)
//   HAMMING_DISTANCE_THRESHOLD - relative to base seq

    option<bool>  multiple_dates{*this, "multiple-dates"};

    option<str>  name{*this, 'n', "name", dflt{""}};
};

// export
//   fasta (phylip?)
//   aligned
//   aa/nuc
//   with deletions inserted
//   wrapped
//   Truncate or extend with - all sequences to make them all of the same length, most common among original sequences.
//   sort by date or name
//   name format: {seq_id} {hi_name_or_seq_name_with_passage} {name} {date} {lab_id} {passage} {lab}
//   name encode
// select
//   Base sequence to export together with other sequences (regex -> to select just one sequence)
//   HAMMING_DISTANCE_THRESHOLD - relative to base seq

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        std::string db_filename{acmacs::acmacsd_root() + "/data/seqdb.json.xz"};
        if (!opt.db->empty())
            db_filename = opt.db;
        seqdb::Seqdb seqdb(db_filename);

        const auto init = [&] {
            if (!opt.name->empty())
                return seqdb.select_by_name(*opt.name);
            else
                return seqdb.all();
        };

        const auto fix_date = [](std::string_view source) -> std::string {
            return source.empty() ? std::string{} : date::display(date::from_string(source, date::allow_incomplete::yes), date::allow_incomplete::yes);
        };

        const auto fix_country = [](std::string_view source) -> acmacs::uppercase {
            if (source == "USA" || source == "US")
                return "UNITED STATES OF AMERICA";
            if (source == "UK" || source == "GB" || source == "GREAT BRITAIN")
                return "UNITED KINGDOM";
            return source;
        };

        std::vector<seqdb::subset::amino_acid_at_pos0_t> aa_at_pos;
        if (!opt.aa_at_pos->empty()) {
            const auto fields = acmacs::string::split(*opt.aa_at_pos, ",");
            aa_at_pos = fields | ranges::view::transform([](const auto& source) -> seqdb::subset::amino_acid_at_pos0_t {
                            if (source.size() >= 2 && source.size() <= 4 && std::isdigit(source.front()) && std::isalpha(source.back()))
                                return {string::from_chars<size_t>(source.substr(0, source.size() - 1)) - 1, source.back(), true};
                            else if (source.size() >= 3 && source.size() <= 5 && source.front() == '!' && std::isdigit(source[1]) && std::isalpha(source.back()))
                                return {string::from_chars<size_t>(source.substr(1, source.size() - 2)) - 1, source.back(), false};
                            else
                                throw std::runtime_error{fmt::format("--aa-at: cannot parse entry: {}", source)};
                        });
        }

        init()
            .subtype(acmacs::uppercase{*opt.subtype})
            .lineage(acmacs::uppercase{*opt.lineage})
            .lab(acmacs::uppercase{*opt.lab})
            .host(acmacs::uppercase{*opt.host})
            .dates(fix_date(opt.start_date), fix_date(opt.end_date))
            .continent(acmacs::uppercase{*opt.continent})
            .country(fix_country(opt.country))
            .clade(acmacs::uppercase{*opt.clade})
            .aa_at_pos(aa_at_pos)
            .multiple_dates(opt.multiple_dates)
            .with_hi_name(opt.with_hi_name)
            .names_matching_regex(opt.name_regex)
            .recent(opt.recent)
            .random(opt.random)
            .print();

        return 0;
    }
    catch (std::exception& err) {
        fmt::print(stderr, "ERROR: {}\n", err);
        return 1;
    }
}

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
