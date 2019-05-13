#include <iostream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/stream.hh"
#include "seqdb-3/fasta.hh"

// ----------------------------------------------------------------------

using namespace acmacs::argv;
struct Options : public argv
{
    Options(int a_argc, const char* const a_argv[], on_error on_err = on_error::exit) : argv() { parse(a_argc, a_argv, on_err); }

    argument<str_array> filenames{*this, arg_name{"filename"}, mandatory};
};

int main(int argc, char* const argv[])
{
    try {
        Options opt(argc, argv);

        for (const auto& filename : *opt.filenames) {
            const auto data = acmacs::seqdb::fasta_scan(filename);
            std::cout << data << '\n';
        }
        return 0;
    }
    catch (std::exception& err) {
        std::cerr << "ERROR: " << err.what() << '\n';
        return 1;
    }
}

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
