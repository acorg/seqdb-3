#include <algorithm>
#include <random>
#include <regex>
#include <numeric>
#include <memory>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/acmacsd.hh"
#include "acmacs-virus/virus-name.hh"
#include "acmacs-chart-2/chart.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/seqdb-parse.hh"
#include "seqdb-3/hamming-distance.hh"

// ----------------------------------------------------------------------

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

static std::string sSeqdbFilename = acmacs::acmacsd_root() + "/data/seqdb.json.xz";

#pragma GCC diagnostic pop

void acmacs::seqdb::v3::setup(std::string_view filename)
{
    if (!filename.empty())
        sSeqdbFilename = filename;

} // acmacs::seqdb::v3::setup

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::Seqdb& acmacs::seqdb::v3::Seqdb::get()
{
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif
    static Seqdb sSeqdb(sSeqdbFilename);
#pragma GCC diagnostic pop

    return sSeqdb;

} // acmacs::seqdb::v3::get

// ----------------------------------------------------------------------

acmacs::seqdb::v3::Seqdb::Seqdb(std::string_view filename)
{
    try {
        json_text_ = static_cast<std::string>(acmacs::file::read(filename));
        parse(json_text_, entries_);
    }
    catch (std::exception&) {
        json_text_.clear();
        entries_.clear();
    }

} // acmacs::seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::Seqdb::Seqdb(std::string&& source)
//     : json_text_(std::move(source))
// {
//     parse(json_text_, entries_);

// } // acmacs::seqdb::v3::Seqdb::Seqdb

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::all() const
{
    subset ss;
    ss.refs_.reserve(entries_.size() * 2);
    for (const auto& entry : entries_) {
        for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no)
            ss.refs_.emplace_back(&entry, seq_no);
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::all

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_name(std::string_view name) const
{
    subset ss;
    if (const auto found = std::lower_bound(std::begin(entries_), std::end(entries_), name, [](const auto& entry, std::string_view nam) { return entry.name < nam; }); found != std::end(entries_) && found->name == name) {
        for (size_t seq_no = 0; seq_no < found->seqs.size(); ++seq_no)
            ss.refs_.emplace_back(&*found, seq_no);
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::select_by_regex(std::string_view re) const
{
    std::regex reg(std::begin(re), std::end(re), std::regex_constants::icase);
    subset ss;
    for (const auto& entry : entries_) {
        for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
            if (ref candidate{&entry, seq_no}; std::regex_search(candidate.full_name(), reg))
                ss.refs_.push_back(std::move(candidate));
        }
    }
    return ss;

} // acmacs::seqdb::v3::Seqdb::select_by_regex

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::seq_id_index_t& acmacs::seqdb::v3::Seqdb::seq_id_index() const
{
    if (seq_id_index_.empty()) {
        seq_id_index_.reserve(entries_.size() * 2);
        for (const auto& entry : entries_) {
            for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
                ref rf{&entry, seq_no};
                seq_id_index_.emplace(rf.seq_id(), std::move(rf));
            }
        }
        seq_id_index_.sort_by_key();
    }
    return seq_id_index_;

} // acmacs::seqdb::v3::Seqdb::seq_id_index

// ----------------------------------------------------------------------

const acmacs::seqdb::v3::hi_name_index_t& acmacs::seqdb::v3::Seqdb::hi_name_index() const
{
    if (hi_name_index_.empty()) {
        hi_name_index_.reserve(entries_.size() * 2);
        for (const auto& entry : entries_) {
            for (size_t seq_no = 0; seq_no < entry.seqs.size(); ++seq_no) {
                for (const auto& hi_name : entry.seqs[seq_no].hi_names) {
                    ref rf{&entry, seq_no};
                    hi_name_index_.emplace(hi_name, std::move(rf));
                }
            }
        }
        hi_name_index_.sort_by_key();
    }
    return hi_name_index_;

} // acmacs::seqdb::v3::Seqdb::hi_name_index

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::Seqdb::match(const acmacs::chart::Antigens& aAntigens, std::string_view /*aChartVirusType*/) const
{
    // check lineage?
    // check virus type

    subset result;
    const auto& hi_name_ind = hi_name_index();
    size_t matched = 0;
    for (auto antigen : aAntigens) {
        if (auto found_ref = hi_name_ind.find(antigen->full_name()); found_ref != hi_name_ind.end()) {
            result.refs_.push_back(std::move(found_ref->second));
            ++matched;
        }
        else if (antigen->passage().empty()) {
            bool found = false;
            for (const auto& selected : select_by_name(antigen->name())) {
                if (selected.seq().has_reassortant(*antigen->reassortant())) {
                    result.refs_.push_back(selected);
                    ++matched;
                    found = true;
                }
            }
            if (!found)
                result.refs_.emplace_back();
        }
        else
            result.refs_.emplace_back();
    }
    fmt::print("INFO: antigens from chart have sequences in seqdb: {}\n", matched);

    return result;

} // acmacs::seqdb::v3::Seqdb::match

// ----------------------------------------------------------------------

std::string_view acmacs::seqdb::v3::SeqdbEntry::host() const
{
    if (const auto ho = acmacs::virus::host(acmacs::virus::v2::virus_name_t{name}); !ho.empty())
        return ho;
    else
        return "HUMAN";

} // acmacs::seqdb::v3::SeqdbEntry::host

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::ref::seq_id() const
{
    const auto to_remove = [](char cc) {
        switch (cc) {
          case '(':
          case ')':
          case '[':
          case ']':
          case ':':
          case '\'':
          case ';':
              // not in seqdb 2019-07-01, but remove them just in case
          case '!':
          case '#':
          case '*':
          case '@':
          case '$':
              return true;
        }
        return false;
    };

    const auto to_replace_with_underscore = [](char cc) {
        switch (cc) {
          case ' ':
          case '&':
          case '=':
              return true;
        }
        return false;
    };

    const auto to_replace_with_slash = [](char cc) {
        switch (cc) {
          case ',':
          case '+':
              return true;
        }
        return false;
    };

    auto source = ::string::join(" ", {entry->name, seq().designation()});
    if (entry->seqs.size() > 1) {
        // there could be multiple seqs with the same designation, but seq_id must be unique, also garli does not like name duplicates
        std::vector<std::string> designations(entry->seqs.size());
        std::transform(std::begin(entry->seqs), std::end(entry->seqs), std::begin(designations), [](const auto& en) { return en.designation(); });
        if (std::count(std::begin(designations), std::end(designations), designations[seq_index]) > 1)
            source.append(fmt::format("_d{}", seq_index));
    }
    return source
            | ranges::view::remove_if(to_remove)
            | ranges::view::replace('?', 'x')
            | ranges::view::replace_if(to_replace_with_slash, '/') // usually in passages
            | ranges::view::replace_if(to_replace_with_underscore, '_')
            ;

} // acmacs::seqdb::v3::ref::seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::multiple_dates(bool do_filter)
{
    if (do_filter)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.entry->dates.size() < 2; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::multiple_dates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::subtype(const acmacs::uppercase& virus_type)
{
    if (!virus_type.empty()) {
        std::string_view vt = virus_type;
        if (vt == "H1")
            vt = "A(H1N1)";
        else if (vt == "H3")
            vt = "A(H3N2)";
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [vt](const auto& en) { return en.entry->virus_type != vt; }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::subtype

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::lineage(const acmacs::uppercase& lineage)
{
    if (!lineage.empty()) {
        std::string_view lin = lineage;
        switch (lin[0]) {
          case 'V': lin = "VICTORIA"; break;
          case 'Y': lin = "YAMAGATA"; break;
        }
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lin](const auto& en) { return en.entry->lineage != lin; }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::lineage

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::lab(const acmacs::uppercase& lab)
{
    if (!lab.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lab=static_cast<std::string_view>(lab)](const auto& en) { return !en.has_lab(lab); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::lab

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::host(const acmacs::uppercase& host)
{
    if (!host.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [host=static_cast<std::string_view>(host)](const auto& en) { return en.entry->host() != host; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::host

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::continent(const acmacs::uppercase& continent)
{
    if (!continent.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [continent=static_cast<std::string_view>(continent)](const auto& en) { return en.entry->continent != continent; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::continent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::country(const acmacs::uppercase& country)
{
    if (!country.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [country=static_cast<std::string_view>(country)](const auto& en) { return en.entry->country != country; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::country

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::clade(const acmacs::uppercase& clade)
{
    if (!clade.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [clade=static_cast<std::string_view>(clade)](const auto& en) { return !en.has_clade(clade); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::clade

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent(size_t recent)
{
    if (recent > 0 && refs_.size() > recent) {
        sort_by_date_recent_first();
        refs_.erase(std::next(std::begin(refs_), static_cast<ssize_t>(recent)), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::random(size_t random)
{
    if (random > 0 && refs_.size() > random) {
        std::mt19937 generator{std::random_device()()};
        std::uniform_int_distribution<size_t> distribution(0, refs_.size() - 1);
        for (size_t no = 0; no < random; ++no)
            refs_[distribution(generator)].to_be_removed = true;
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.to_be_removed; }), std::end(refs_));
        std::for_each(std::begin(refs_), std::end(refs_), [](auto& en) { en.to_be_removed = false; });
    }
    return *this;

} // acmacs::seqdb::v3::subset::random

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::group_by_hamming_distance(size_t dist_threshold)
{
    if (dist_threshold > 0) {
        auto group_master = most_recent_with_hi_name();
        auto group_master_aa_aligned = group_master->seq().aa_aligned();
        auto group_first = std::begin(refs_);
        acmacs::Counter<ssize_t> counter_group_size;
        for (size_t group_no = 1; group_first != std::end(refs_); ++group_no) {
            if (group_no != 1) {
                group_master = group_first;
                while (!group_master->has_hi_names() && group_master != std::end(refs_))
                    ++group_master;
                if (group_master == std::end(refs_))
                    group_master = group_first;
                group_master_aa_aligned = group_master->seq().aa_aligned();
            }
            // fmt::print("DEBUG: group {} master: {} {} rest size: {}\n", group_no, group_master->seq_id(), group_master->entry->date(), std::end(refs_) - group_first);
            for (auto refp = group_first; refp != std::end(refs_); ++refp)
                refp->hamming_distance = hamming_distance(group_master_aa_aligned, refp->seq().aa_aligned());
            std::sort(group_first, std::end(refs_), [](const auto& e1, const auto& e2) { return e1.hamming_distance < e2.hamming_distance; });
            for (auto refp = group_first; ; ++refp) {
                if (refp == std::end(refs_) || refp->hamming_distance >= dist_threshold) {
                    // fmt::print(stderr, "DEBUG: group {} size: {}\n", group_no, refp - group_first);
                    counter_group_size.count(refp - group_first);
                    group_first = refp;
                    break;
                }
                else
                    refp->group_no = group_no;
            }
        }
        fmt::print(stderr, "DEBUG: (group-size, num-groups): {}\n", counter_group_size);
        fmt::print(stderr, "DEBUG: total groups: {}\n", refs_.back().group_no);
    }
    return *this;

} // acmacs::seqdb::v3::subset::group_by_hamming_distance

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::refs_t::const_iterator acmacs::seqdb::v3::subset::most_recent_with_hi_name() const
{
    refs_t::const_iterator result = std::end(refs_);
    std::string_view date;
    for (auto refp = std::begin(refs_); refp != std::end(refs_); ++refp) {
        if (refp->has_hi_names() && refp->entry->date() > date) { // refs_[no].seq().reassortants.empty() &&
            result = refp;
            date = refp->entry->date();
        }
    }
    return result;

} // acmacs::seqdb::v3::subset::most_recent_with_hi_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::with_hi_name(bool with_hi_name)
{
    if (with_hi_name)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.has_hi_names(); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::with_hi_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::aa_at_pos(const std::vector<amino_acid_at_pos0_t>& aa_at_pos0)
{
    if (!aa_at_pos0.empty()) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&aa_at_pos0](const auto& en) {
                                       try {
                                           constexpr bool keep = false, remove = true;
                                           if (en.seq().amino_acids.empty())
                                               return remove;
                                           for (const auto& [pos0, aa, equal] : aa_at_pos0) {
                                               if ((en.seq().aa_at(pos0) == aa) != equal)
                                                   return remove;
                                           }
                                           return keep;
                                       }
                                       catch (std::exception& err) {
                                           throw std::runtime_error{fmt::format("{}, full_name: {}", err, en.full_name())};
                                       }
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::aa_at_pos

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::names_matching_regex(const std::vector<std::string_view>& regex_list)
{
    if (!regex_list.empty()) {
        std::vector<std::regex> re_list(regex_list.size());
        std::transform(std::begin(regex_list), std::end(regex_list), std::begin(re_list),
                       [](const auto& regex_s) { return std::regex(std::begin(regex_s), std::end(regex_s), std::regex_constants::icase); });
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&re_list](const auto& en) {
                                       return std::none_of(std::begin(re_list), std::end(re_list), [full_name = en.full_name()](const auto& re) { return std::regex_search(full_name, re); });
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::names_matching_regex

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::dates(std::string_view start, std::string_view end)
{
    if (!start.empty() || !end.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [start,end](const auto& en) { return !en.entry->date_within(start, end); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::dates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend_single_matching(std::string_view re, const Seqdb& seqdb)
{
    if (!re.empty()) {
        auto candidates = seqdb.select_by_regex(re);
        if (candidates.size() != 1)
            throw std::runtime_error{fmt::format("regular expression must select single sequence: \"{}\", selected: {}", re, candidates.size())};
        refs_.erase(std::remove(std::begin(refs_), std::end(refs_), candidates.front()), std::end(refs_)); // remove it, if selected earlier
        refs_.insert(std::begin(refs_), candidates.front());
    }
    return *this;

} // acmacs::seqdb::v3::subset::prepend_single_matching

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base(size_t threshold, bool do_filter)
{
    if (do_filter) {
        refs_.erase(std::remove_if(std::next(std::begin(refs_)), std::end(refs_),
                                   [threshold, base_seq = refs_.front().seq().aa_aligned()](const auto& en) { return hamming_distance(en.seq().aa_aligned(), base_seq) >= threshold; }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_hamming_distance_to_base

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::sort(sorting srt)
{
    switch (srt) {
        case sorting::none:
            break;
        case sorting::name_asc:
            sort_by_name_asc();
            break;
        case sorting::name_desc:
            sort_by_name_desc();
            break;
        case sorting::date_asc:
            sort_by_date_oldest_first();
            break;
        case sorting::date_desc:
            sort_by_date_recent_first();
            break;
    }
    return *this;

} // acmacs::seqdb::v3::subset::sort

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_stat(bool do_report)
{
    if (do_report) {
        if (!refs_.empty()) {
            size_t with_hi_names = 0;
            std::string_view min_date = refs_.front().entry->date(), max_date = min_date;
            for (const auto& ref : refs_) {
                const auto date = ref.entry->date();
                if (date < min_date)
                    min_date = date;
                else if (date > max_date)
                    max_date = date;
                if (!ref.seq().hi_names.empty())
                    ++with_hi_names;
            }
            fmt::print(stderr, "Sequences: {}\nDate range: {} - {}\nHiDb matches: {}\n", refs_.size(), min_date, max_date, with_hi_names);
        }
        else {
            fmt::print(stderr, "No sequences selected\n");
        }
    }

    return *this;

} // acmacs::seqdb::v3::subset::report_stat

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::export_sequences(std::string_view filename, const export_options& options)
{
    if (!filename.empty()) {
        auto to_export = export_collect(options);

        if (options.e_most_common_length) {
            const acmacs::Counter counter(to_export, [](const auto& en) { return en.second.size(); });
            ranges::for_each(to_export, [most_common_length = counter.max().first](auto& en) { en.second.resize(most_common_length, '-'); });
        }

        switch (options.e_format) {
            case export_options::format::fasta_aa:
            case export_options::format::fasta_nuc:
                acmacs::file::write(filename, export_fasta(to_export, options));
                break;
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::export_sequences

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::make_name(std::string_view name_format, const ref& entry) const
{
    return fmt::format(name_format,
                       fmt::arg("seq_id", entry.seq_id()),
                       fmt::arg("full_name", entry.full_name()),
                       fmt::arg("hi_name_or_full_name", entry.hi_name_or_full_name()),
                       fmt::arg("lineage", entry.entry->lineage),
                       fmt::arg("name", entry.entry->name),
                       fmt::arg("date", entry.entry->date()),
                       fmt::arg("dates", entry.entry->dates),
                       fmt::arg("lab_id", entry.seq().lab_id()),
                       fmt::arg("passage", entry.seq().passage()),
                       fmt::arg("clades", entry.seq().clades),
                       fmt::arg("lab", entry.seq().lab()),
                       fmt::arg("country", entry.entry->country),
                       fmt::arg("continent", entry.entry->continent));

} // acmacs::seqdb::v3::subset::make_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::collected_t acmacs::seqdb::v3::subset::export_collect(const export_options& options) const
{
    const auto get_seq = [&options](const auto& entry) -> std::string_view {
        if (options.e_format == export_options::format::fasta_aa) {
            if (options.e_aligned)
                return entry.seq().aa_aligned();
            else
                return entry.seq().amino_acids;
        }
        else {
            if (options.e_aligned)
                return entry.seq().nuc_aligned();
            else
                return entry.seq().nucs;
        }
    };

    collected_t result(refs_.size()); // {seq_id, sequence}
    std::transform(std::begin(refs_), std::end(refs_), std::begin(result),
                   [this, &options, &get_seq](const auto& en) -> collected_entry_t { return std::pair(make_name(options.e_name_format, en), std::string{get_seq(en)}); });
    return result;

} // acmacs::seqdb::v3::subset::export_collect

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::export_fasta(const collected_t& entries, const export_options& options)
{
    std::string output;
    const auto output_size =
        std::accumulate(std::begin(entries), std::end(entries), 0UL, [](size_t size, const auto& en) { return size + en.first.size() + en.second.size() + 2 + en.second.size() / 40; });
    output.reserve(output_size);
    for (const auto& en : entries) {
        output.append(1, '>');
        output.append(en.first);
        output.append(1, '\n');
        if (options.e_wrap_at == 0 || options.e_wrap_at >= en.second.size()) {
            output.append(en.second);
            output.append(1, '\n');
        }
        else {
            for (const auto chunk : en.second | ranges::view::chunk(options.e_wrap_at)) {
                output.append(chunk);
                output.append(1, '\n');
            }
        }
    }
    return output;

} // acmacs::seqdb::v3::subset::export_fasta

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::print(std::string_view name_format, bool do_print)
{
    if (do_print) {
        for (const auto& ref : *this)
            fmt::print("{}\n", make_name(name_format, ref));
    }
    return *this;

} // acmacs::seqdb::v3::subset::print

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
