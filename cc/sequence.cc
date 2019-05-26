#include <cctype>
#include <map>
#include <tuple>
#include <array>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/range.hh"
#include "acmacs-base/algorithm.hh"
#include "seqdb-3/sequence.hh"

// ----------------------------------------------------------------------

static std::vector<std::pair<char, size_t>> symbol_frequences(std::string_view seq);
static std::string translate_nucleotides_to_amino_acids(std::string_view nucleotides, size_t offset);

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::sequence_t::align(std::string_view type_subtype_hint, std::string_view debug_name)
{
    if (aa_.empty())
        return false;
    if (type_subtype_hint == "A(H3N2)")
        return align_h3n2(debug_name) || align_any(debug_name, type_subtype_hint);
    else
        return align_any(debug_name, type_subtype_hint);

} // acmacs::seqdb::v3::sequence_t::align

// ----------------------------------------------------------------------

// Some sequences from CNIC (and perhaps from other labs) have initial
// part of nucleotides with stop codons inside. To figure out correct
// translation we have first to translate with all possible offsets
// (0, 1, 2) and not stoppoing at stop codons, then try to align all
// of them. Most probably just one offset leads to finding correct
// align shift.
void acmacs::seqdb::v3::sequence_t::translate()
{
    constexpr size_t MINIMUM_SEQUENCE_AA_LENGTH = 400; // throw away everything shorter

    if (!nuc_.empty()) {
        using translated_t = std::tuple<std::string, size_t>; // translated longest part, nuc prefix size at which translation started

        const auto transformation = [this](size_t offset) -> translated_t {
            const auto find_longest_part = [offset](std::string aa) -> translated_t {
                auto split_data = acmacs::string::split(aa, "*");
                const auto longest_part = std::max_element(std::begin(split_data), std::end(split_data), [](const auto& e1, const auto& e2) { return e1.size() < e2.size(); });
                return {std::string(*longest_part), offset + static_cast<size_t>(longest_part->data() - aa.data()) * 3};
            };
            return find_longest_part(translate_nucleotides_to_amino_acids(this->nuc_, offset));
        };

        std::array<translated_t, 3> translated;
        std::transform(acmacs::index_iterator(0UL), acmacs::index_iterator(translated.size()), std::begin(translated), transformation);
        const auto longest_translated = std::max_element(std::begin(translated), std::end(translated), [](const auto& e1, const auto& e2) { return std::get<std::string>(e1).size() < std::get<std::string>(e2).size(); });
        if (std::get<std::string>(*longest_translated).size() >= MINIMUM_SEQUENCE_AA_LENGTH) {
            std::tie(aa_, nuc_translation_offset_) = *longest_translated;
        }
    }

} // acmacs::seqdb::v3::sequence_t::translate

// ----------------------------------------------------------------------

#include "acmacs-base/global-constructors-push.hh"
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
#include "acmacs-base/diagnostics-pop.hh"

std::string translate_nucleotides_to_amino_acids(std::string_view nucleotides, size_t offset)
{
    using diff_t = decltype(CODON_TO_PROTEIN)::difference_type;

    std::string result((nucleotides.size() - offset) / 3 + 1, '-');
    auto result_p = result.begin();
    for (auto off = offset; off < nucleotides.size(); off += 3, ++result_p) {
        if (const auto it = CODON_TO_PROTEIN.find(std::string_view(nucleotides.data() + static_cast<diff_t>(off), 3)); it != CODON_TO_PROTEIN.end())
            *result_p = it->second;
        else
            *result_p = 'X';
    }
    result.resize(static_cast<size_t>(result_p - result.begin()));
    return result;

} // translate_nucleotides_to_amino_acids

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::sequence_t::import(std::string_view source)
{
    nuc_.resize(source.size(), '-');
    std::transform(std::begin(source), std::end(source), std::begin(nuc_), [](char c) { return std::toupper(c); });

    auto freq = symbol_frequences(nuc_);
    ranges::sort(freq, [](const auto& e1, const auto& e2) { return e1.second > e2.second; }); // most frequent first
    // std::sort(freq.begin(), freq.end(), [](const auto& e1, const auto& e2) { return e1.second > e2.second; }); // most frequent first

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
    }
    else {
        // if (freq.size() < 12)
        //    fmt::print(stderr, "nuc freq: {} {}\n", nuc_.size(), acmacs::to_string(freq));
        aa_ = nuc_;
        nuc_.clear();
    }

} // acmacs::seqdb::v3::sequence_t::import

// ----------------------------------------------------------------------

std::vector<std::pair<char, size_t>> symbol_frequences(std::string_view seq)
{
    std::vector<std::pair<char, size_t>> result;
    for (auto cc : seq) {
        if (auto found = std::find_if(std::begin(result), std::end(result), [cc](const auto& en) { return en.first == cc; }); found != std::end(result))
            ++found->second;
        else
            result.emplace_back(cc, 1UL);
    }
    return result;

} // symbol_frequences

// ----------------------------------------------------------------------

const char* const sH3N2_align_1 = "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENC";

bool acmacs::seqdb::v3::sequence_t::align_h3n2(std::string_view debug_name)
{
    const std::string_view pattern_1{sH3N2_align_1};
    const auto aa_start = aa_.find(pattern_1.substr(0, 2));
    if (aa_start == std::string::npos) {
        // fmt::print(stderr, "NOT H3? {}\n{}\n", debug_name, std::string_view(aa_.data(), 100));
        return false;
    }
    const auto hamd = ::string::hamming_distance(pattern_1, std::string_view(aa_.data() + aa_start, pattern_1.size()));
    if (hamd > 10)
        fmt::print(stderr, "hamm {} {}\n{}\n", hamd, debug_name, std::string_view(aa_.data(), 120));

    return false;

} // acmacs::seqdb::v3::sequence_t::align_h3n2

// ----------------------------------------------------------------------

bool acmacs::seqdb::v3::sequence_t::align_any(std::string_view debug_name, std::string_view except)
{
    // if (except != "A(H3N2)" && align_h3n2(debug_name))
    //     return true;
    return false;

} // acmacs::seqdb::v3::sequence_t::align_any

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
