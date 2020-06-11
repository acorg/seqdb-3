#include "acmacs-base/argv.hh"
#include "acmacs-base/fmt.hh"
#include "acmacs-base/date.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-whocc-data/labs.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/log.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    option<str> db{*this, "db"};

    // select
    option<str_array> seq_id{*this, "seq-id", desc{"initially filter by seq-id, all matching"}};
    option<str>       seq_id_from{*this, "seq-id-from", desc{"read list of seq ids from a file (one per line) and initially select them all"}};
    option<str>       name{*this, 'n', "name", desc{"initially filter by name (name only, full string equality)"}};
    option<str>       names_from{*this, "names-from", desc{"read names from a file (one per line)\n                                       and initially select them all (name only, full string equality)"}};
    option<str>       accession_numbers_from{*this, "accession-numbers-from", desc{"read accession numbers (gisaid and/or ncbi) names from a file (one per line)\n                                       and initially select them all (full string equality)"}};
    option<str>       subtype{*this, "flu", desc{"B, A(H1N1), H1, A(H3N2), H3"}};
    option<str>       host{*this, "host"};
    option<str>       lab{*this, "lab"};
    option<bool>      whocc_lab{*this, "whocc-lab", desc{"only 4 WHOCC labs"}};
    option<str>       lineage{*this, "lineage"};
    option<str>       start_date{*this, "start-date"};
    option<str>       end_date{*this, "end-date"};
    option<str>       continent{*this, "continent", desc{"africa antarctica asia australia-oceania central-america europe middle-east north-america russia south-america"}};
    option<str>       country{*this, "country"};
    option<str>       clade{*this, "clade"};
    option<str>       aa_at_pos{*this, "aa-at-pos", desc{"comma separated list: 162N,74R,!167X"}};
    option<str>       nuc_at_pos{*this, "nuc-at-pos", desc{"comma separated list: 618C"}};
    option<size_t>    recent{*this, "recent", dflt{0UL}};
    option<str>       recent_matched{*this, "recent-matched", desc{"num1,num2 - select num1 most recent,\n                                       then add num2 older which are also matched against hidb"}};
    option<size_t>    random{*this, "random", dflt{0UL}};
    option<bool>      with_hi_name{*this, "with-hi-name", desc{"matched against hidb"}};
    option<str_array> name_regex{*this, "re", desc{"filter names by regex, multiple regex possible, all matching listed"}};
    option<str_array> prepend{*this, "prepend", desc{"prepend with seq by seq-id, multiple possible, always included"}};
    option<str>       base_seq_id{*this, "base-seq-id", desc{"single base sequence (outgroup), always included"}};
    option<size_t>    nuc_hamming_distance_threshold{*this, "nuc-hamming-distance-threshold", dflt{140UL}, desc{"Select only sequences having hamming distance to the base sequence less than threshold."}};
    option<bool>      multiple_dates{*this, "multiple-dates"};
    option<str>       sort_by{*this, "sort", dflt{"none"}, desc{"none, name, -name, date, -date"}};
    option<str>       name_format{*this, 'f', "name-format", desc{"{{seq_id}} {{full_name}} {{hi_name_or_full_name}} {{hi_names}} {{hi_name}} {{lineage}} {{name}}\n                                       {{date}} {{dates}} {{lab_id}} {{passage}} {{clades}} {{lab}} {{country}} {{continent}} {{group_no}}\n                                       {{hamming_distance}} {{nuc_length}} {{aa_length}} {{gisaid_accession_numbers}} {{ncbi_accession_numbers}}"}};
    option<size_t>    group_by_hamming_distance{*this, "group-by-hamming", dflt{0ul}, desc{"Group sequences by hamming distance."}};
    option<bool>      subset_by_hamming_distance_random{*this, "subset-by-hamming-random", desc{"Subset using davipatti algorithm 2019-07-23."}};
    option<bool>      remove_nuc_duplicates{*this, "remove-nuc-duplicates", desc{""}};
    option<bool>      keep_all_hi_matched{*this, "keep-all-hi", desc{"do NOT remove HI matched when removing duplicates (--remove-nuc-duplicates)"}};
    option<size_t>    output_size{*this, "output-size", dflt{4000ul}, desc{"Number of sequences to use from grouped by hamming distance."}};
    option<size_t>    minimum_aa_length{*this, "minimum-aa-length", dflt{0ul}, desc{"Select only sequences having min number of AAs in alignment."}};
    option<size_t>    minimum_nuc_length{*this, "minimum-nuc-length", dflt{0ul}, desc{"Select only sequences having min number of nucs in alignment."}};

    // print
    option<bool>      print{*this, 'p', "print", desc{"force printing selected sequences"}};
    option<bool>      report_hamming_distance{*this, "report-hamming", desc{"Report hamming distance from base for all strains."}};
    option<str>       report_aa_at{*this, "report-aa-at", desc{"comma separated list: 142,144."}};
    option<bool>      no_stat{*this, "no-stat"};

    // export
    option<str>       fasta{*this, "fasta", desc{"export to fasta, - for stdout"}};
    option<bool>      wrap{*this, "wrap"};
    option<bool>      nucs{*this, "nucs", desc{"export nucleotide sequences instead of amino acid"}};
    option<bool>      not_aligned{*this, "not-aligned", desc{"do not align for exporting"}};
    option<bool>      most_common_length{*this, "most-common-length", desc{"truncate or extend with - all sequences to make them all of the same length,\n                                       most common among original sequences"}};
    option<size_t>    length{*this, "length", dflt{0ul}, desc{"truncate or extend with - all sequences to make them all of the same length,\n                                       0 - do not truncate/extend"}};

    option<str_array> verbose{*this, 'v', "verbose", desc{"comma separated list (or multiple switches) of enablers"}};
};

int main(int argc, char* const argv[])
{
    using namespace std::string_view_literals;

    try {
        acmacs::log::register_enabler_acmacs_base();
        acmacs::log::register_enabler("seq"sv, acmacs::log::sequences);
        acmacs::log::register_enabler("fasta"sv, acmacs::log::fasta);
        Options opt(argc, argv);
        acmacs::log::enable(opt.verbose);

        acmacs::seqdb::setup(opt.db);
        const auto& seqdb = acmacs::seqdb::get();

        const auto init = [&] {
            if (!opt.seq_id->empty())
                return seqdb.select_by_seq_id(*opt.seq_id);
            else if (opt.seq_id_from)
                return seqdb.select_by_seq_id(acmacs::string::split(static_cast<std::string>(acmacs::file::read(opt.seq_id_from)), "\n", acmacs::string::Split::StripRemoveEmpty));
            else if (opt.name)
                return seqdb.select_by_name(*opt.name);
            else if (opt.names_from)
                return seqdb.select_by_name(acmacs::string::split(static_cast<std::string>(acmacs::file::read(opt.names_from)), "\n", acmacs::string::Split::StripRemoveEmpty));
            else if (opt.accession_numbers_from)
                return seqdb.select_by_accession_number(acmacs::string::split(static_cast<std::string>(acmacs::file::read(opt.accession_numbers_from)), "\n", acmacs::string::Split::StripRemoveEmpty));
            else
                return seqdb.all();
        };

        const auto fix_date = [](std::string_view source) -> std::string {
            return source.empty() ? std::string{} : date::display(date::from_string(source, date::allow_incomplete::yes), date::allow_incomplete::yes);
        };

        const auto fix_country = [](const acmacs::uppercase& source) -> acmacs::uppercase {
            if (*source == "USA" || *source == "US")
                return "UNITED STATES OF AMERICA";
            if (*source == "UK" || *source == "GB" || *source == "GREAT BRITAIN")
                return "UNITED KINGDOM";
            return source;
        };

        const auto sorting_order = [](const acmacs::lowercase& desc) -> acmacs::seqdb::subset::sorting {
            if (desc == acmacs::lowercase{"none"})
                return acmacs::seqdb::subset::sorting::none;
            if (desc == acmacs::lowercase{"name"})
                return acmacs::seqdb::subset::sorting::name_asc;
            if (desc == acmacs::lowercase{"-name"})
                return acmacs::seqdb::subset::sorting::name_desc;
            if (desc == acmacs::lowercase{"date"})
                return acmacs::seqdb::subset::sorting::date_asc;
            if (desc == acmacs::lowercase{"-date"})
                return acmacs::seqdb::subset::sorting::date_desc;
            fmt::print(stderr, "WARNING: unrecognized soriting: {}\n", desc);
            return acmacs::seqdb::subset::sorting::name_asc;
        };

        acmacs::seqdb::amino_acid_at_pos1_eq_list_t aa_at_pos;
        if (!opt.aa_at_pos->empty())
            aa_at_pos = acmacs::seqdb::extract_aa_at_pos1_eq_list(*opt.aa_at_pos);

        acmacs::seqdb::nucleotide_at_pos1_eq_list_t nuc_at_pos;
        if (!opt.nuc_at_pos->empty())
            nuc_at_pos = acmacs::seqdb::extract_nuc_at_pos1_eq_list(*opt.nuc_at_pos);

        acmacs::seqdb::pos1_list_t aa_at_pos_report;
        if (!opt.report_aa_at->empty())
            aa_at_pos_report = acmacs::seqdb::extract_pos1_list(*opt.report_aa_at);

        if (opt.name_format->empty()) {
            if (opt.fasta->empty())
                opt.name_format.add("\"{full_name}\" {lineage} {dates} {country} {clades} \"{lab}\" {seq_id}");
            else
                opt.name_format.add("{seq_id}");
        }

        init()
            .subtype(acmacs::uppercase{*opt.subtype})
            .lineage(acmacs::uppercase{*opt.lineage})
            .lab(acmacs::whocc::lab_name_normalize(*opt.lab))
            .whocc_lab(opt.whocc_lab)
            .host(acmacs::uppercase{*opt.host})
            .dates(fix_date(opt.start_date), fix_date(opt.end_date))
            .continent(acmacs::uppercase{*opt.continent})
            .country(fix_country(acmacs::uppercase{*opt.country}))
            .clade(seqdb, acmacs::uppercase{*opt.clade})
            .aa_at_pos(seqdb, aa_at_pos)
            .nuc_at_pos(seqdb, nuc_at_pos)
            .min_aa_length(seqdb, opt.minimum_aa_length)
            .min_nuc_length(seqdb, opt.minimum_nuc_length)
            .multiple_dates(opt.multiple_dates)
            .with_hi_name(opt.with_hi_name)
            .names_matching_regex(opt.name_regex)
            .recent(opt.recent, opt.remove_nuc_duplicates ? acmacs::seqdb::subset::master_only::yes : acmacs::seqdb::subset::master_only::no)
            .recent_matched(acmacs::string::split_into_size_t(*opt.recent_matched, ","), opt.remove_nuc_duplicates ? acmacs::seqdb::subset::master_only::yes : acmacs::seqdb::subset::master_only::no)
            .random(opt.random)
            .group_by_hamming_distance(seqdb, opt.group_by_hamming_distance, opt.output_size)
            .subset_by_hamming_distance_random(seqdb, opt.subset_by_hamming_distance_random, opt.output_size)
            .remove_nuc_duplicates(opt.remove_nuc_duplicates, opt.keep_all_hi_matched)
            .remove_empty(seqdb, opt.nucs)
            .sort(sorting_order(opt.sort_by))
            .report_stat(seqdb, !opt.no_stat) // static_cast<bool>(opt.fasta))
            .report_aa_at(seqdb, aa_at_pos_report)
            .prepend(opt.prepend, seqdb)
            .prepend(opt.base_seq_id, seqdb)
            // .prepend_single_matching(opt.base_seq_regex, seqdb)
            .nuc_hamming_distance_to_base(opt.nuc_hamming_distance_threshold, !!opt.base_seq_id)
            .export_sequences(opt.fasta, seqdb,
                              acmacs::seqdb::export_options{}
                                  .fasta(opt.nucs)
                                  .wrap(opt.wrap ? 80 : 0)
                                  .aligned(opt.not_aligned ? acmacs::seqdb::export_options::aligned::no : acmacs::seqdb::export_options::aligned::yes)
                                  .most_common_length(opt.most_common_length ? acmacs::seqdb::export_options::most_common_length::yes : acmacs::seqdb::export_options::most_common_length::no)
                                  .length(opt.length)
                                  .name_format(opt.name_format))         // acmacs::seqdb::v3::subset::make_name
            .print(seqdb, opt.name_format, opt.print /* || opt.fasta */) // acmacs::seqdb::v3::subset::make_name
            .report_hamming_distance(opt.report_hamming_distance && !opt.base_seq_id->empty());

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
