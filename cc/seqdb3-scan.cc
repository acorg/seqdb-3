#include <iostream>

#include "acmacs-base/argv.hh"
#include "acmacs-base/read-file.hh"
#include "acmacs-base/enumerate.hh"
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

        std::vector<std::vector<acmacs::seqdb::v3::fasta::sequence_t>> sequences_per_file(opt.filenames->size());
#pragma omp parallel for default(shared) schedule(static, 4)
        for (size_t f_no = 0; f_no < opt.filenames->size(); ++f_no) {
            const auto& filename = (*opt.filenames)[f_no];
            try {
                const std::string file_data_s = acmacs::file::read(filename);
                const std::string_view file_data = file_data_s;
                acmacs::seqdb::v3::fasta::scan_input_t file_input{std::begin(file_data), std::end(file_data)};
                while (!file_input.done()) {
                    acmacs::seqdb::v3::fasta::scan_output_t sequence_ref;
                    std::tie(file_input, sequence_ref) = acmacs::seqdb::v3::fasta::scan(file_input);

                    std::optional<acmacs::seqdb::v3::fasta::sequence_t> seq;
                    for (auto parser : {&acmacs::seqdb::v3::fasta::name_gisaid_spaces, &acmacs::seqdb::v3::fasta::name_gisaid_underscores, &acmacs::seqdb::v3::fasta::name_plain}) {
                        seq = (*parser)(sequence_ref.name, filename, file_input.line_no);
                        if (seq.has_value())
                            break;
                    }
                    if (seq.has_value()) {
                        acmacs::seqdb::v3::fasta::normalize_name(*seq, filename, file_input.line_no);
                        seq->sequence = acmacs::seqdb::v3::fasta::normalize_sequence(sequence_ref.sequence, filename, file_input.line_no);
                        sequences_per_file[f_no].push_back(*seq);
                    }
                    else
                        std::cerr << "WARNING: " << filename << ':' << file_input.line_no << ": unable to parse name: " << sequence_ref.name << '\n';
                }
            }
            catch (std::exception& err) {
                throw std::runtime_error(string::concat(filename, err.what()));
            }
        }

        std::vector<acmacs::seqdb::v3::fasta::sequence_t> all_sequences;
        for (auto& per_file : sequences_per_file)
            std::move(per_file.begin(), per_file.end(), std::back_inserter(all_sequences)); 
        std::cout << all_sequences.size() << " sequences read\n";

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
