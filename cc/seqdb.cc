#include <algorithm>
#include <random>
#include <regex>

#include "acmacs-base/read-file.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-virus/virus-name.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/seqdb-parse.hh"
#include "seqdb-3/hamming-distance.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::v3::Seqdb::Seqdb(const std::string& filename)
{
    json_text_ = static_cast<std::string>(acmacs::file::read(filename));
    parse(json_text_, entries_);

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

    const auto source = ::string::join(" ", {entry->name, ::string::join(" ", seq().reassortants), seq().passages.empty() ? std::string_view{} : seq().passages.front()});
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
            refs_[distribution(generator)].selected = true;
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.selected; }), std::end(refs_));
        std::for_each(std::begin(refs_), std::end(refs_), [](auto& en) { en.selected = false; });
    }
    return *this;

} // acmacs::seqdb::v3::subset::random

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::with_hi_name(bool with_hi_name)
{
    if (with_hi_name)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.seq().hi_names.empty(); }), std::end(refs_));
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

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::export_sequences(std::string_view filename, const export_options& options)
{
    if (!filename.empty()) {
        std::string output;
        output.reserve(refs_.size() * static_cast<size_t>(refs_.front().seq_id().size() * 1.5 + refs_.front().seq().nucs.size()));
        for (const auto& en : refs_) {
            switch (options.e_format) {
              case export_options::format::fasta_aa:
              case export_options::format::fasta_nuc:
                  export_fasta(en, options, output);
                  break;
            }
        }
        acmacs::file::write(filename, output);
    }
    return *this;

} // acmacs::seqdb::v3::subset::export_sequences

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::subset::export_fasta(const ref& entry, const export_options& options, std::string& output)
{
    const auto get_seq = [&entry, &options] {
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

    output.append(entry.seq_id());
    output.append(1, '\n');
    const auto seq = get_seq();
    if (options.e_wrap_at == 0 || options.e_wrap_at >= seq.size()) {
        output.append(seq);
        output.append(1, '\n');
    }
    else {
        for (size_t start = 0; start < seq.size(); start += options.e_wrap_at) {
            output.append(seq.substr(start, options.e_wrap_at));
            output.append(1, '\n');
        }
    }

} // acmacs::seqdb::v3::subset::export_fasta

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::print(print_options po, bool do_print)
{
    const auto make_details = [](const auto& ref) {
        fmt::print("{}{}{} {} {} {}\n", ref.full_name(), ref.entry->lineage.empty() ? "" : " L:", ref.entry->lineage.empty() ? std::string_view{} : ref.entry->lineage, ref.entry->dates,
                   ref.entry->country, ref.seq().clades);
    };

    const auto make_seq_id = [](const auto& ref) { fmt::print("{}\n", ref.seq_id()); };
    const auto make_passage = [](const auto& ref) { for (const auto& passage : ref.seq().passages) fmt::print("{}\n", passage); };

    if (do_print) {
        for (const auto& ref : *this) {
            switch (po) {
                case print_options::details:
                    make_details(ref);
                    break;
                case print_options::seq_id:
                    make_seq_id(ref);
                    break;
                case print_options::passage:
                    make_passage(ref);
                    break;
            }
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::print

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
