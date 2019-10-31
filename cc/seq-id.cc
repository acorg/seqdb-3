#include "acmacs-base/range-v3.hh"
#include "seqdb-3/seq-id.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::seq_id_t acmacs::seqdb::v3::make_seq_id(std::string_view designation)
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

    return seq_id_t{ranges::to<std::string>(
        designation
        | ranges::views::remove_if(to_remove)
        | ranges::views::replace('?', 'x')
        | ranges::views::replace_if(to_replace_with_slash, '/') // usually in passages
        | ranges::views::replace_if(to_replace_with_underscore, '_')
                                            )};

} // acmacs::seqdb::v3::make_seq_id

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
