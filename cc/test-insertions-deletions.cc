// #include "acmacs-base/fmt.hh"
#include "seqdb-3/scan-deletions.hh"

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    if (argc != 3) {
        fmt::print(stderr, "Usage {} <master-seq> <to-align-seq>\n", argv[0]);
        return 1;
    }

    const auto res = acmacs::seqdb::scan::deletions_insertions(argv[1], argv[2], acmacs::debug::yes);
    fmt::print("{}\n{}\n", acmacs::seqdb::scan::format_aa(res.insertions, argv[1], '.'), acmacs::seqdb::scan::format_aa(res.deletions, argv[2], '.'));
}

// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
