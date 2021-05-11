#include <cctype>
#include <map>
#include <tuple>
#include <array>
#include <algorithm>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/string-join.hh"
#include "acmacs-base/algorithm.hh"
#include "acmacs-base/hash.hh"
#include "acmacs-base/counter.hh"
#include "seqdb-3/scan-sequence.hh"
#include "seqdb-3/scan-align.hh"

// ----------------------------------------------------------------------

namespace local
{
    static std::string translate_nucleotides_to_amino_acids(std::string_view nucleotides, size_t offset);
}

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::format_aa(const std::vector<deletions_insertions_t::pos_num_t>& pos_num, std::string_view sequence, char deletion_symbol)
{
    fmt::memory_buffer out;
    pos0_t pos{0};
    for (const auto& en : pos_num) {
        fmt::format_to(out, "{}{}", sequence.substr(*pos, *(en.pos - pos)), std::string(en.num, deletion_symbol));
        pos = en.pos;
    }
    fmt::format_to(out, "{}", sequence.substr(*pos));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::format

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::format(const deletions_insertions_t& deletions)
{
    fmt::memory_buffer out;
    const auto frmt = [&out](const char* prefix, const auto& num_pos) {
        if (!num_pos.empty()) {
            fmt::format_to(out, "{}[{}](", prefix, num_pos.size());
            bool first = true;
            for (const auto& en : num_pos) {
                if (first)
                    first = false;
                else
                    fmt::format_to(out, " ");
                fmt::format_to(out, "{}:{}", *(en.pos + 1UL), en.num);
            }
            fmt::format_to(out, ")");
        }
    };

    frmt("DEL", deletions.deletions);
    frmt(" INS", deletions.insertions);
    fmt::format_to(out, "<pos-1-based>");
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::format

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::sequence_t::aa_format() const
{
    fmt::memory_buffer out;
    const auto aa = aa_aligned();
    pos0_t pos{0};
    for (const auto& en : deletions_.deletions) {
        fmt::format_to(out, "{}{:->{}s}", aa.substr(*pos, *(en.pos - pos)), "", en.num);
        pos = en.pos;
    }
    fmt::format_to(out, "{}", aa.substr(*pos));
    // fmt::print(stderr, "aa_format s:{}\n  {}\n  {}\n", pos, aa, fmt::to_string(out));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::sequence_t::aa_format

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::sequence_t::aa_format_not_aligned() const
{
    std::string_view aav{aa_};
    fmt::memory_buffer out;
    if (shift_aa_ > shift_t{0}) {
        fmt::format_to(out, "{}", aav.substr(0, *shift_aa_));
        aav.remove_prefix(*shift_aa_);
    }
    pos0_t pos{0};
    for (const auto& en : deletions_.deletions) {
        fmt::format_to(out, "{}{:->{}s}", aav.substr(*pos, *(en.pos - pos)), "", en.num);
        pos = en.pos;
    }
    fmt::format_to(out, "{}", aav.substr(*pos));
    // fmt::print(stderr, "aa_format_not_aligned s:{}\n  {}\n  {}\n", shift_aa_, aa_, fmt::to_string(out));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::sequence_t::aa_format_not_aligned

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::sequence_t::nuc_format() const
{
    fmt::memory_buffer out;
    const auto nuc = nuc_aligned();
    pos0_t pos{0};
    for (const auto& en : deletions_.deletions) {
        fmt::format_to(out, "{}{:->{}s}", nuc.substr(*pos, *(en.pos.aa_to_nuc() - pos)), "", en.num * 3);
        pos = en.pos.aa_to_nuc();
    }
    fmt::format_to(out, "{}", nuc.substr(*pos));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::sequence_t::nuc_format

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::sequence_t::nuc_format_not_aligned() const
{
    std::string_view nucv{nuc_};
    fmt::memory_buffer out;
    if (shift_nuc_ > shift_t{0}) {
        fmt::format_to(out, "{}", nucv.substr(0, *shift_nuc_));
        nucv.remove_prefix(*shift_nuc_);
    }
    pos0_t pos{0};
    for (const auto& en : deletions_.deletions) {
        fmt::format_to(out, "{}{:->{}s}", nucv.substr(*pos, *(en.pos.aa_to_nuc() - pos)), "", en.num * 3);
        pos = en.pos.aa_to_nuc();
    }
    fmt::format_to(out, "{}", nucv.substr(*pos));
    return fmt::to_string(out);

} // acmacs::seqdb::v3::scan::sequence_t::nuc_format_not_aligned

// ----------------------------------------------------------------------

// Some sequences from CNIC (and perhaps from other labs) have initial
// part of nucleotides with stop codons inside. To figure out correct
// translation we have first to translate with all possible offsets
// (0, 1, 2) and not stoppoing at stop codons, then try to align all
// of them. Most probably just one offset leads to finding correct
// align shift.
void acmacs::seqdb::v3::scan::sequence_t::translate()
{
    constexpr size_t MINIMUM_SEQUENCE_AA_LENGTH = 200; // throw away everything shorter, HA1 is kinda 318, need to have just HA1 sequences to be able to make HA1 trees

    if (!nuc_.empty()) {
        using translated_t = std::tuple<std::string, size_t>; // translated longest part, nuc prefix size at which translation started

        const auto transformation = [this](size_t offset) -> translated_t {
            const auto find_longest_part = [offset](std::string aa) -> translated_t {
                auto split_data = acmacs::string::split(aa, "*");
                const auto longest_part = std::max_element(std::begin(split_data), std::end(split_data), [](const auto& e1, const auto& e2) { return e1.size() < e2.size(); });
                return {std::string(*longest_part), offset + static_cast<size_t>(longest_part->data() - aa.data()) * 3};
            };
            return find_longest_part(local::translate_nucleotides_to_amino_acids(this->nuc_, offset));
        };

        std::array<translated_t, 3> translated;
        ranges::transform(ranges::views::iota(0UL, translated.size()), std::begin(translated), transformation);
        const auto longest_translated =
            std::max_element(std::begin(translated), std::end(translated), [](const auto& e1, const auto& e2) { return std::get<std::string>(e1).size() < std::get<std::string>(e2).size(); });
        if (std::get<std::string>(*longest_translated).size() >= MINIMUM_SEQUENCE_AA_LENGTH)
            std::tie(aa_, nuc_translation_offset_) = *longest_translated;
    }

    aa_trim_absent();

} // acmacs::seqdb::v3::scan::sequence_t::translate

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::sequence_t::date_simulated() const noexcept
{
    if (!dates_.empty())
        return dates_.front();
    else if (auto yr = acmacs::virus::year(name()); yr.has_value())
        return fmt::format("{}-01-01", *yr);
    else
        return "1800-01-01";

} // acmacs::seqdb::v3::scan::sequence_t::date_simulated

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::aa_trim_absent()
{
    if (!aa_.empty()) {
        // remove trailing X and - in aa
        if (const auto found = aa_.find_last_not_of("X-"); found != std::string::npos)
            aa_.erase(found + 1);
        else
            fmt::print(stderr, "WARNING: just X and - in AA sequence for {} ::: {}\n", full_name(), aa_);

        // remove leading X and -
        if (const auto found = aa_.find_first_not_of("X-"); found > 0 && found != std::string::npos) {
            // fmt::print(stderr, "leading X: {} ::: {}\n", full_name(), aa_);
            aa_.erase(0, found);
            nuc_translation_offset_ += found * 3;
        }
    }

} // acmacs::seqdb::v3::scan::sequence_t::aa_trim_absent

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::scan::sequence_t::full_name() const
{
    return acmacs::string::join(acmacs::string::join_space, *name(), *reassortant(), annotations(), passages_.empty() ? std::string{} : *passages_.front(), *lineage());

} // acmacs::seqdb::v3::scan::sequence_t::full_name

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::add_lab_id(const acmacs::uppercase& lab, const acmacs::uppercase& lab_id)
{
    if (lab.empty() && lab_id.empty())
        return;

    if (auto found = lab_ids_.find(lab); found == lab_ids_.end())
        lab_ids_.emplace(lab, std::set<acmacs::uppercase>{lab_id});
    else
        found->second.insert(lab_id);

} // acmacs::seqdb::v3::scan::sequence_t::add_lab_id

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::add_lab_id(const acmacs::uppercase& lab)
{
    if (auto found = lab_ids_.find(lab); found == lab_ids_.end())
        lab_ids_.emplace(lab, std::set<acmacs::uppercase>{});

} // acmacs::seqdb::v3::scan::sequence_t::add_lab_id

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::add_date(const std::string& date)
{
    if (!date.empty() && (dates_.empty() || !empty_month_or_day(date))) {
        dates_.add(date, flat_set_sort_afterwards::yes);
    }

} // acmacs::seqdb::v3::scan::sequence_t::add_date

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::scan::sequence_t::lab_in(std::initializer_list<std::string_view> labs) const
{
    return std::any_of(std::begin(labs), std::end(labs), [this](const auto lab) { return lab_ids_.find(lab) != lab_ids_.end(); });

} // acmacs::seqdb::v3::scan::sequence_t::lab_in

// ----------------------------------------------------------------------

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
#endif

static const std::map<std::string_view, char> CODON_TO_PROTEIN = {
    {"UGC", 'C'}, {"GTA", 'V'}, {"GTG", 'V'}, {"CCT", 'P'}, {"CUG", 'L'}, {"AGG", 'R'}, {"CTT", 'L'}, {"CUU", 'L'},
    {"CTG", 'L'}, {"GCU", 'A'}, {"CCG", 'P'}, {"AUG", 'M'}, {"GGC", 'G'}, {"UUA", 'L'}, {"GAG", 'E'}, {"UGG", 'W'},
    {"UUU", 'F'}, {"UUG", 'L'}, {"ACU", 'T'}, {"TTA", 'L'}, {"AAT", 'N'}, {"CGU", 'R'}, {"CCA", 'P'}, {"GCC", 'A'},
    {"GCG", 'A'}, {"TTG", 'L'}, {"CAT", 'H'}, {"AAC", 'N'}, {"GCA", 'A'}, {"GAU", 'D'}, {"UAU", 'Y'}, {"CAC", 'H'},
    {"AUA", 'I'}, {"GUC", 'V'}, {"TCG", 'S'}, {"GGG", 'G'}, {"AGC", 'S'}, {"CTA", 'L'}, {"GCT", 'A'}, {"CCC", 'P'},
    {"ACC", 'T'}, {"GAT", 'D'}, {"TCC", 'S'}, {"UAC", 'Y'}, {"CAU", 'H'}, {"UCG", 'S'}, {"CAA", 'Q'}, {"UCC", 'S'},
    {"AGU", 'S'}, {"TTT", 'F'}, {"ACA", 'T'}, {"ACG", 'T'}, {"CGC", 'R'}, {"TGT", 'C'}, {"CAG", 'Q'}, {"GUA", 'V'},
    {"GGU", 'G'}, {"AAG", 'K'}, {"AGA", 'R'}, {"ATA", 'I'}, {"TAT", 'Y'}, {"UCU", 'S'}, {"TCA", 'S'}, {"GAA", 'E'},
    {"AGT", 'S'}, {"TCT", 'S'}, {"ACT", 'T'}, {"CGA", 'R'}, {"GGT", 'G'}, {"TGC", 'C'}, {"UGU", 'C'}, {"CUC", 'L'},
    {"GAC", 'D'}, {"UUC", 'F'}, {"GTC", 'V'}, {"ATT", 'I'}, {"TAC", 'Y'}, {"CUA", 'L'}, {"TTC", 'F'}, {"GTT", 'V'},
    {"UCA", 'S'}, {"AUC", 'I'}, {"GGA", 'G'}, {"GUG", 'V'}, {"GUU", 'V'}, {"AUU", 'I'}, {"CGT", 'R'}, {"CCU", 'P'},
    {"ATG", 'M'}, {"AAA", 'K'}, {"TGG", 'W'}, {"CGG", 'R'}, {"AAU", 'N'}, {"CTC", 'L'}, {"ATC", 'I'},
    {"TAA", '*'}, {"UAA", '*'}, {"TAG", '*'}, {"UAG", '*'}, {"TGA", '*'}, {"UGA", '*'}, {"TAR", '*'}, {"TRA", '*'}, {"UAR", '*'}, {"URA", '*'},
};

#pragma GCC diagnostic pop

std::string local::translate_nucleotides_to_amino_acids(std::string_view nucleotides, size_t offset)
{
    using diff_t = decltype(CODON_TO_PROTEIN)::difference_type;

    std::string result((nucleotides.size() - offset) / 3 + 1, '-');
    auto result_p = result.begin();
    for (auto off = offset; off < (nucleotides.size() - 2); off += 3, ++result_p) {
        if (const auto it = CODON_TO_PROTEIN.find(std::string_view(nucleotides.data() + static_cast<diff_t>(off), 3)); it != CODON_TO_PROTEIN.end())
            *result_p = it->second;
        else
            *result_p = 'X';
    }
    result.resize(static_cast<size_t>(result_p - result.begin()));
    return result;

} // local::translate_nucleotides_to_amino_acids

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::import(std::string_view source)
{
    nuc_.resize(source.size(), '-');
    std::transform(std::begin(source), std::end(source), std::begin(nuc_), [](char c) { return std::toupper(c); });

    const auto freq = acmacs::CounterChar(nuc_).pairs(acmacs::CounterChar::sorted::yes);

    const auto most_freq_are_acgnt = [](const auto& frq) {
        const auto most_frequent_symbols_size = std::min(5UL, frq.size());
        std::string most_frequent_symbols(most_frequent_symbols_size, ' ');
        const auto mfs_end = acmacs::transform_if(frq.begin(), frq.begin() + static_cast<ssize_t>(most_frequent_symbols_size), most_frequent_symbols.begin(), [](const auto& e1) { return e1.second > 5; }, [](const auto& e2) { return e2.first; });
        most_frequent_symbols.resize(static_cast<size_t>(mfs_end - most_frequent_symbols.begin()));
        std::sort(std::begin(most_frequent_symbols), std::end(most_frequent_symbols));
        return std::string_view(most_frequent_symbols.data(), 4) == "ACGT" || most_frequent_symbols == "ACGNT" || most_frequent_symbols == "-ACGT";
    };

    if (freq.size() > 1 && /* freq.size() < 12 && */ ((freq[0].second > (nuc_.size() / 4) && freq[1].second > (nuc_.size() / 5)) || most_freq_are_acgnt(freq))) {
        // looks like nuc
        hash_ = acmacs::hash(nuc_);
    }
    else {
        // if (freq.size() < 12)
        //    fmt::print(stderr, "nuc freq: {} {}\n", nuc_.size(), acmacs::to_string(freq));
        aa_ = nuc_;
        nuc_.clear();
        hash_ = acmacs::hash(aa_);
    }

} // acmacs::seqdb::v3::scan::sequence_t::import

// ----------------------------------------------------------------------

acmacs::seqdb::v3::scan::sequence_t acmacs::seqdb::v3::scan::sequence_t::from_aligned_aa(const acmacs::virus::name_t& name, std::string_view source)
{
    sequence_t result;
    result.name_ = name;
    result.aa_ = source;
    result.shift_aa_ = shift_t{0};
    return result;

} // acmacs::seqdb::v3::scan::sequence_t::from_aligned_aa

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::set_shift(int shift_aa, std::optional<acmacs::virus::type_subtype_t> type_subtype)
{
    if (shift_aa < 0) {
        aa_.insert(0, static_cast<size_t>(-shift_aa), 'X');
        shift_aa_ = shift_t{0};
        nuc_.insert(0, static_cast<size_t>(-shift_aa) * 3, '-');
        shift_nuc_ = shift_t{nuc_translation_offset_};
    }
    else {
        shift_aa_ = shift_t{shift_aa};
        shift_nuc_ = shift_t{nuc_translation_offset_ + shift_aa * 3};
    }
    if (type_subtype.has_value()) {
        type_subtype_ = *type_subtype;
        acmacs::virus::set_type_subtype(name_, type_subtype_);
    }
    remove_issue_not_aligned();

} // acmacs::seqdb::v3::scan::sequence_t::set_shift

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::merge_from(const sequence_t& source)
{
    if (!source.country_.empty()) {
        if (country_.empty())
            country_ = source.country_;
        else if (country_ != source.country_)
            fmt::print(stderr, "WARNING: sequence_t::merge_from: {}: countries differ: \"{}\" vs. \"{}\"\n", name(), country_, source.country_);
    }

    if (!source.continent_.empty()) {
        if (continent_.empty())
            continent_ = source.continent_;
        else if (continent_ != source.continent_)
            fmt::print(stderr, "WARNING: sequence_t::merge_from: {}: continents differ: \"{}\" vs. \"{}\"\n", name(), continent_, source.continent_);
    }

    dates_.merge_from(source.dates_);
    passages_.merge_from(source.passages_, flat_set_sort_afterwards::yes);
    hi_names_.merge_from(source.hi_names_);
    isolate_id_.merge_from(source.isolate_id_);
    submitters_.merge_from(source.submitters_);
    gisaid_last_modified_.merge_from(source.gisaid_last_modified_);
    originating_lab_.merge_from(source.originating_lab_);
    gisaid_segment_.merge_from(source.gisaid_segment_);
    gisaid_segment_number_.merge_from(source.gisaid_segment_number_);
    gisaid_identifier_.merge_from(source.gisaid_identifier_);
    gisaid_dna_accession_no_.merge_from(source.gisaid_dna_accession_no_);
    gisaid_dna_insdc_.merge_from(source.gisaid_dna_insdc_);

    for (const auto& [lab, ids] : source.lab_ids()) {
        if (const auto found = lab_ids_.find(lab); found != lab_ids_.end()) {
            for (const auto& id : ids)
                found->second.insert(id);
        }
        else
            lab_ids_.emplace(lab, ids);
    }

} // acmacs::seqdb::v3::scan::sequence_t::merge_from

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::scan::sequence_t::update_subtype(const acmacs::virus::type_subtype_t& subtype)
{
    if (subtype.size() >= 7 && name_->find('/') == 5) { // A(H3) -> A(H3N2)
        name_.get().replace(0, 5, subtype);
        type_subtype_ = subtype;
    }

} // acmacs::seqdb::v3::scan::sequence_t::update_subtype

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
