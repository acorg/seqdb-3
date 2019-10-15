#include <cctype>

#include "acmacs-base/rjson.hh"
#include "acmacs-base/string.hh"
#include "seqdb-3/aa-at-pos.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::amino_acid_at_pos1_list_t acmacs::seqdb::v3::extract_aa_at_pos1(const rjson::value& source)
{
    amino_acid_at_pos1_list_t pos1_aa;
    rjson::for_each(source, [&pos1_aa](const rjson::value& entry) {
        const auto text{entry.to<std::string_view>()};
        if (text.size() >= 2 && text.size() <= 4 && std::isdigit(text.front()) && std::isalpha(text.back()))
            pos1_aa.emplace_back(::string::from_chars<size_t>(text.substr(0, text.size() - 1)), text.back(), true);
        else if (text.size() >= 3 && text.size() <= 5 && text.front() == '!' && std::isdigit(text[1]) && std::isalpha(text.back()))
            pos1_aa.emplace_back(::string::from_chars<size_t>(text.substr(1, text.size() - 2)), text.back(), false);
        else
            throw std::exception{};
    });
    return pos1_aa;

} // acmacs::seqdb::v3::extract_aa_at_pos

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
