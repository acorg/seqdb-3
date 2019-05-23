#include <iostream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/enumerate.hh"
#include "locationdb/locdb.hh"
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

        for (size_t f_no = 0; f_no < opt.filenames->size(); ++f_no) {
            const auto& filename = (*opt.filenames)[f_no];
            try {
                const std::string file_data_s = acmacs::file::read(filename);
                const std::string_view file_data = file_data_s;
                acmacs::seqdb::v3::fasta::scan_input_t file_input{std::begin(file_data), std::end(file_data)};
                while (!file_input.done()) {
                    acmacs::seqdb::v3::fasta::scan_output_t sequence_ref;
                    try {
                        std::tie(file_input, sequence_ref) = acmacs::seqdb::v3::fasta::scan(file_input);

                        std::optional<acmacs::seqdb::v3::fasta::sequence_t> seq;
                        for (auto parser : {&acmacs::seqdb::v3::fasta::name_gisaid_spaces, &acmacs::seqdb::v3::fasta::name_gisaid_underscores, &acmacs::seqdb::v3::fasta::name_plain}) {
                            seq = (*parser)(sequence_ref.name, {}, filename, file_input.name_line_no);
                            if (seq.has_value())
                                break;
                        }
                        if (seq.has_value()) {
                            try {
                                const auto messages = acmacs::seqdb::v3::fasta::normalize_name(*seq);
                                for (const auto& msg : messages)
                                    std::cerr << "WARNING: " << filename << ':' << file_input.name_line_no << ": " << msg << '\n';

                                std::cout << seq->name;
                                if (!seq->reassortant.empty())
                                    std::cout << " R:" << seq->reassortant;
                                std::cout << '\n';
                            }
                            catch (std::exception& err) {
                                std::cerr << "ERROR: " << filename << ':' << file_input.name_line_no << ": unable to parse name: " << sequence_ref.name << ": " << err.what() << '\n';
                            }
                        }
                        else
                            std::cerr << "WARNING: " << filename << ':' << file_input.name_line_no << ": unable to parse fasta name: " << sequence_ref.name << '\n';
                    }
                    catch (std::exception& err) {
                        throw std::runtime_error(string::concat(file_input.name_line_no, ": ", err.what()));
                    }
                }
            }
            catch (std::exception& err) {
                throw std::runtime_error(string::concat(filename, ':', err.what()));
            }
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
