#include <cctype>

#include "acmacs-base/rjson-v2.hh"
#include "acmacs-base/string.hh"
#include "acmacs-base/string-split.hh"
#include "acmacs-base/string-from-chars.hh"
#include "seqdb-3/aa-at-pos.hh"

// ----------------------------------------------------------------------

template <typename AA_NUC_LIST> struct list_pos_conv
{
};

template <> struct list_pos_conv<acmacs::seqdb::amino_acid_at_pos1_eq_list_t>
{
    using at_pos1_eq_t = acmacs::seqdb::amino_acid_at_pos1_eq_t;
};

template <> struct list_pos_conv<acmacs::seqdb::nucleotide_at_pos1_eq_list_t>
{
    using at_pos1_eq_t = acmacs::seqdb::nucleotide_at_pos1_eq_t;
};

// ----------------------------------------------------------------------

template <typename R, size_t MIN_SIZE, size_t MAX_SIZE> inline R extract_aa_nuc_at_pos1_eq(std::string_view source)
{
    if (source.size() >= MIN_SIZE && source.size() <= MAX_SIZE && std::isdigit(source.front()) && (std::isalpha(source.back()) || source.back() == '-'))
        return {acmacs::seqdb::pos1_t{acmacs::string::from_chars<size_t>(source.substr(0, source.size() - 1))}, source.back(), true};
    else if (source.size() >= (MIN_SIZE + 1) && source.size() <= (MAX_SIZE + 1) && source.front() == '!' && std::isdigit(source[1]) && (std::isalpha(source.back()) || source.back() == '-'))
        return {acmacs::seqdb::pos1_t{acmacs::string::from_chars<size_t>(source.substr(1, source.size() - 2))}, source.back(), false};
    else
        throw acmacs::seqdb::extract_at_pos_error{fmt::format("invalid aa/nuc-pos: \"{}\" (expected 183P or !183P)", source)};

} // acmacs::seqdb::v3::extract_aa_at_pos1_eq

// ----------------------------------------------------------------------

template <typename R, size_t MIN_SIZE, size_t MAX_SIZE> inline R extract_aa_nuc_at_pos1_eq_list(const rjson::value& source)
{
    return std::visit(
        []<typename T>(T && arg)->R {
            if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
                return extract_aa_nuc_at_pos1_eq_list<R, MIN_SIZE, MAX_SIZE>(std::string_view{arg});
            }
            else if constexpr (std::is_same_v<std::decay_t<T>, rjson::array>) {
                R pos1_aa_eq;
                arg.for_each([&pos1_aa_eq](const rjson::value& entry) { pos1_aa_eq.push_back(extract_aa_nuc_at_pos1_eq<typename list_pos_conv<R>::at_pos1_eq_t, MIN_SIZE, MAX_SIZE>(entry.to<std::string_view>())); });
                return pos1_aa_eq;
            }
            else
                throw acmacs::seqdb::extract_at_pos_error{fmt::format("invalid aa/nuc-at-pos1 list: {}", arg)};
        },
        source.val_());

} // acmacs::seqdb::v3::extract_aa_at_pos_eq_list

// ----------------------------------------------------------------------

template <typename R, size_t MIN_SIZE, size_t MAX_SIZE> inline R extract_aa_nuc_at_pos1_eq_list(std::string_view source)
{
    const auto fields = acmacs::string::split(source, acmacs::string::Split::RemoveEmpty);
    R pos1_aa_eq{fields.size()};
    std::transform(std::begin(fields), std::end(fields), std::begin(pos1_aa_eq),
                   [](std::string_view field) { return extract_aa_nuc_at_pos1_eq<typename list_pos_conv<R>::at_pos1_eq_t, MIN_SIZE, MAX_SIZE>(field); });
    return pos1_aa_eq;

} // acmacs::seqdb::v3::extract_aa_at_pos_eq_list

// ======================================================================

acmacs::seqdb::amino_acid_at_pos1_eq_t acmacs::seqdb::v3::extract_aa_at_pos1_eq(std::string_view source)
{
    return extract_aa_nuc_at_pos1_eq<acmacs::seqdb::amino_acid_at_pos1_eq_t, 2, 4>(source);

} // acmacs::seqdb::v3::extract_aa_at_pos1_eq

// ----------------------------------------------------------------------

acmacs::seqdb::amino_acid_at_pos1_eq_list_t acmacs::seqdb::v3::extract_aa_at_pos1_eq_list(const rjson::value& source)
{
    return extract_aa_nuc_at_pos1_eq_list<acmacs::seqdb::amino_acid_at_pos1_eq_list_t, 2, 4>(source);

} // acmacs::seqdb::v3::extract_aa_at_pos_eq_list

// ----------------------------------------------------------------------

acmacs::seqdb::amino_acid_at_pos1_eq_list_t acmacs::seqdb::v3::extract_aa_at_pos1_eq_list(std::string_view source)
{
    return extract_aa_nuc_at_pos1_eq_list<acmacs::seqdb::amino_acid_at_pos1_eq_list_t, 2, 4>(source);

} // acmacs::seqdb::v3::extract_aa_at_pos_eq_list

// ----------------------------------------------------------------------

acmacs::seqdb::nucleotide_at_pos1_eq_t acmacs::seqdb::v3::extract_nuc_at_pos1_eq(std::string_view source)
{
    return extract_aa_nuc_at_pos1_eq<acmacs::seqdb::nucleotide_at_pos1_eq_t, 2, 5>(source);

} // acmacs::seqdb::v3::extract_nuc_at_pos1_eq

// ----------------------------------------------------------------------

acmacs::seqdb::nucleotide_at_pos1_eq_list_t acmacs::seqdb::v3::extract_nuc_at_pos1_eq_list(const rjson::value& source)
{
    return extract_aa_nuc_at_pos1_eq_list<acmacs::seqdb::nucleotide_at_pos1_eq_list_t, 2, 5>(source);

} // acmacs::seqdb::v3::extract_nuc_at_pos1_eq_list

// ----------------------------------------------------------------------

acmacs::seqdb::nucleotide_at_pos1_eq_list_t acmacs::seqdb::v3::extract_nuc_at_pos1_eq_list(std::string_view source) // space or comma separated list, e.g. 1703A 384C 618C !1010G"
{
    return extract_aa_nuc_at_pos1_eq_list<acmacs::seqdb::nucleotide_at_pos1_eq_list_t, 2, 5>(source);

} // acmacs::seqdb::v3::extract_nuc_at_pos1_eq_list

// ----------------------------------------------------------------------

acmacs::seqdb::pos1_list_t acmacs::seqdb::v3::extract_pos1_list(std::string_view source)
{
    const auto fields = acmacs::string::split(source, acmacs::string::Split::RemoveEmpty);
    pos1_list_t pos1_list{fields.size(), pos1_t{99999}};
    std::transform(std::begin(fields), std::end(fields), std::begin(pos1_list), [](std::string_view field) { return pos1_t{acmacs::string::from_chars<size_t>(field)}; });
    return pos1_list;

} // acmacs::seqdb::v3::extract_aa_at_pos1_list

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
