#include <cctype>

#include "acmacs-base/rjson.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-split.hh"
#include "seqdb-3/aa-at-pos.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::amino_acid_at_pos1_eq_t acmacs::seqdb::v3::extract_aa_at_pos1_eq(std::string_view source)
{
    if (source.size() >= 2 && source.size() <= 4 && std::isdigit(source.front()) && std::isalpha(source.back()))
        return {acmacs::seqdb::pos1_t{::string::from_chars<size_t>(source.substr(0, source.size() - 1))}, source.back(), true};
    else if (source.size() >= 3 && source.size() <= 5 && source.front() == '!' && std::isdigit(source[1]) && std::isalpha(source.back()))
        return {acmacs::seqdb::pos1_t{::string::from_chars<size_t>(source.substr(1, source.size() - 2))}, source.back(), false};
    else
        throw std::exception{};

} // acmacs::seqdb::v3::extract_aa_at_pos1_eq

// ----------------------------------------------------------------------

acmacs::seqdb::amino_acid_at_pos1_eq_list_t acmacs::seqdb::v3::extract_aa_at_pos1_eq_list(const rjson::value& source)
{
    amino_acid_at_pos1_eq_list_t pos1_aa_eq;
    rjson::for_each(source, [&pos1_aa_eq](const rjson::value& entry) { pos1_aa_eq.push_back(extract_aa_at_pos1_eq(entry.to<std::string_view>())); });
    return pos1_aa_eq;

} // acmacs::seqdb::v3::extract_aa_at_pos_eq_list

// ----------------------------------------------------------------------

acmacs::seqdb::amino_acid_at_pos1_eq_list_t acmacs::seqdb::v3::extract_aa_at_pos1_eq_list(std::string_view source)
{
    // space or comma separated list, e.g. "183P 141E !151K"
    const auto split = [](std::string_view src) {
        if (src.find(",") != std::string_view::npos)
            return acmacs::string::split(src, ",", acmacs::string::Split::RemoveEmpty);
        else
            return acmacs::string::split(src, " ", acmacs::string::Split::RemoveEmpty);
    };

    const auto fields = split(source);
    amino_acid_at_pos1_eq_list_t pos1_aa_eq{fields.size()};
    std::transform(std::begin(fields), std::end(fields), std::begin(pos1_aa_eq), [](std::string_view field) { return extract_aa_at_pos1_eq(field); });
    return pos1_aa_eq;

} // acmacs::seqdb::v3::extract_aa_at_pos_eq_list

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
