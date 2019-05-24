#include <cctype>

#include "acmacs-base/fmt.hh"
#include "acmacs-base/to-string.hh"
#include "acmacs-base/algorithm.hh"
#include "seqdb-3/sequence.hh"

static std::vector<std::pair<char, size_t>> symbol_frequences(std::string_view seq);

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::sequence_t::translate()
{
    if (!nuc_.empty()) {
    }

} // acmacs::seqdb::v3::sequence_t::translate

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::sequence_t::import(std::string_view source)
{
    nuc_.resize(source.size(), '-');
    std::transform(std::begin(source), std::end(source), std::begin(nuc_), [](char c) { return std::toupper(c); });

    auto freq = symbol_frequences(nuc_);
    std::sort(freq.begin(), freq.end(), [](const auto& e1, const auto& e2) { return e1.second > e2.second; }); // most frequent first

    const auto most_freq_are_acgnt = [](const auto& frq) {
        const auto most_frequent_symbols_size = std::min(5UL, frq.size());
        std::string most_frequent_symbols(most_frequent_symbols_size, ' ');
        const auto mfs_end = acmacs::transform_if(frq.begin(), frq.begin() + static_cast<ssize_t>(most_frequent_symbols_size), most_frequent_symbols.begin(), [](const auto& e1) { return e1.second > 5; }, [](const auto& e2) { return e2.first; });
        most_frequent_symbols.resize(static_cast<size_t>(mfs_end - most_frequent_symbols.begin()));
        std::sort(std::begin(most_frequent_symbols), std::end(most_frequent_symbols));
        return std::string_view(most_frequent_symbols.data(), 4) == "ACGT" || most_frequent_symbols == "ACGNT" || most_frequent_symbols == "-ACGT";
    };

    if (freq.size() > 1 && freq.size() < 10 && ((freq[0].second > (nuc_.size() / 4) && freq[1].second > (nuc_.size() / 5)) || most_freq_are_acgnt(freq))) {
        // looks like nuc
    }
    else {
        if (freq.size() < 10)
            fmt::print(stderr, "nuc freq: {} {}\n", nuc_.size(), acmacs::to_string(freq));
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


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
