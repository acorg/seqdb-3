#include "acmacs-base/in-json-parser.hh"
#include "seqdb-3/seqdb-parse.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace local
{
    class labs : public in_json::stack_entry
    {
      public:
        labs(acmacs::seqdb::SeqdbSeq::labs_t& target) : target_{target} {}
        const char* injson_name() override { return "labs"; }
        void injson_put_array() override {}

        void injson_put_key(std::string_view data) override
        {
            in_json::stack_entry::injson_put_key(data);
            target_.emplace_back(data, acmacs::seqdb::SeqdbSeq::lab_ids_t{});
        }

        void injson_put_string(std::string_view data) override { target_.back().second.emplace_back(data); }

      private:
        acmacs::seqdb::SeqdbSeq::labs_t& target_;
    };

    class gisaid_data : public in_json::stack_entry
    {
      public:
        gisaid_data(acmacs::seqdb::SeqdbSeq::gisaid_data_t& target) : target_{target} {}
        const char* injson_name() override { return "gisaid"; }

        void injson_put_array() override
        {
        }

        void injson_pop_array() override
        {
            reset_key();
        }

        void injson_put_string(std::string_view data) override
        {
            switch (key_[0]) {
                case 'i':
                    target_.isolate_ids.emplace_back(data);
                    break;
              case 's':
                    target_.sample_ids_by_sample_provider.emplace_back(data);
                    break;
              case 'S':         // todo
              case 'm':
              case 'o':
              case 'n':
              case 't':
              case 'D':
              case 'd':
                  break;
                default:
                    throw in_json::parse_error(AD_FORMAT("gisaid: unexpected key: \"{}\"", key_));
            }
        }

      private:
        acmacs::seqdb::SeqdbSeq::gisaid_data_t& target_;
    };

    class reference : public in_json::stack_entry
    {
      public:
        reference(acmacs::seqdb::SeqdbSeq::master_ref_t& target) : target_{target} {}
        const char* injson_name() override { return "reference"; }

        void injson_put_string(std::string_view data) override
        {
            switch (key_[0]) {
                case 'N':
                    target_.name = data;
                    break;
                case 'H':
                    target_.hash = data;
                    break;
                // case 'p':
                //     target_.passage = data;
                //     break;
                // case 'r':
                //     target_.reassortant = data;
                //     break;
                // case 'A':
                //     target_.annotations = data;
                //     break;
                default:
                    throw in_json::parse_error(AD_FORMAT("seq reference (\"R\"): unexpected key: \"{}\"", key_));
            }
            reset_key();
        }

      private:
        acmacs::seqdb::SeqdbSeq::master_ref_t& target_;
    };

    class seq : public in_json::stack_entry
    {
      public:
        seq(acmacs::seqdb::SeqdbSeq& target) : target_{target} {}

        const char* injson_name() override { return "seq"; }

        std::unique_ptr<in_json::stack_entry> injson_put_object() override
        {
            switch (key_[0]) {
              case 'l':
                  reset_key();
                  return std::make_unique<labs>(target_.lab_ids);
              case 'G':         // gisaid data
                  reset_key();
                  return std::make_unique<gisaid_data>(target_.gisaid);
              case 'R':         // master with identical sequence reference
                  reset_key();
                  return std::make_unique<reference>(target_.master);
              default:
                  throw in_json::parse_error(AD_FORMAT("seq: unexpected sub-object, key: \"{}\"", key_));
            }
        }

        void injson_put_array() override
        {
            // if (key_.size() != 1 || (key_[0] != 'p' && key_[0] != 'c' && key_[0] != 'h' && key_[0] != 'r'))
            //     throw in_json::parse_error(AD_FORMAT("seq: unexpected array, key: \"{}\"", key_));
        }

        void injson_pop_array() override
        {
            // if (key_.size() != 1 || (key_[0] != 'p' && key_[0] != 'c' && key_[0] != 'h' && key_[0] != 'r'))
            //     throw in_json::parse_error(AD_FORMAT("seq: unexpected array, key: \"{}\"", key_));
            reset_key();
        }

        void injson_put_string(std::string_view data) override
        {
            // else if (key_.size() == 1) {
            switch (key_[0]) {
                case 'p':
                    target_.passages.emplace_back(data);
                    break;
                case 'r':
                    target_.reassortants.emplace_back(data);
                    break;
                case 'c':
                    target_.clades.emplace_back(data);
                    break;
                case 'h':
                    target_.hi_names.emplace_back(data);
                    break;
                // case 'g':
                //     target_.gene = data;
                //     reset_key();
                //     break;
                case 'a':
                    std::get<std::string_view>(target_.amino_acids) = data;
                    reset_key();
                    break;
                case 'n':
                    std::get<std::string_view>(target_.nucs) = data;
                    reset_key();
                    break;
                case 'A':
                    target_.annotations = data;
                    reset_key();
                    break;
                case 'H':
                    target_.hash = data;
                    reset_key();
                    break;
                case 'i':       // issues
                    for (const auto iss : data) {
                        using namespace acmacs::seqdb::sequence;
                        switch (iss) {
                            case 'A':
                                set(target_.issues, issue::not_aligned);
                                break;
                            case 'i':
                                set(target_.issues, issue::has_insertions);
                                break;
                            case 's':
                                set(target_.issues, issue::too_short);
                                break;
                            case 'b':
                                set(target_.issues, issue::garbage_at_the_beginning);
                                break;
                            case 'e':
                                set(target_.issues, issue::garbage_at_the_end);
                                break;
                            default:
                                throw in_json::parse_error(AD_FORMAT("seq: unrecognized issue: {}", iss));
                        }
                    }
                    reset_key();
                    break;
                default:
                    throw in_json::parse_error(AD_FORMAT("seq: unexpected key: \"{}\"", key_));
            }
            // }
            // else
            //     throw in_json::parse_error(AD_FORMAT("seq: unexpected key: \"{}\"", key_));
        }

        void injson_put_integer(std::string_view data) override
        {
            // if (key_.size() == 1) {
            switch (key_[0]) {
                case 's':
                    std::get<acmacs::seqdb::alignment_t>(target_.amino_acids) = acmacs::seqdb::alignment_t{data};
                    reset_key();
                    break;
                case 't':
                    std::get<acmacs::seqdb::alignment_t>(target_.nucs) = acmacs::seqdb::alignment_t{data};
                    reset_key();
                    break;
                default:
                    in_json::stack_entry::injson_put_integer(data);
                    break;
            }
            // }
            // else
            //     in_json::stack_entry::injson_put_integer(data);
        }

      private:
        acmacs::seqdb::SeqdbSeq& target_;
    };

    class entry : public in_json::stack_entry
    {
      public:
        entry(acmacs::seqdb::SeqdbEntry& target) : target_{target} {}

        const char* injson_name() override { return "entry"; }

        std::unique_ptr<in_json::stack_entry> injson_put_object() override
        {
            return std::make_unique<seq>(target_.seqs.emplace_back()); // objects are only under "s"
        }

        void injson_put_array() override
        {
            // fmt::print(stderr, "WARNING: entry::injson_put_array\n");
            // if (key_.size() != 1 || (key_[0] != 'd' && key_[0] != 's'))
            //     throw in_json::parse_error(AD_FORMAT("entry: unexpected array, key: \"{}\"", key_));
        }

        void injson_pop_array() override
        {
            // fmt::print(stderr, "WARNING: entry::injson_pop_array\n");
            // if (key_.size() != 1 || (key_[0] != 'd' && key_[0] != 's'))
            //     throw in_json::parse_error(AD_FORMAT("entry: unexpected array, key: \"{}\"", key_));
            reset_key();
        }

        void injson_put_string(std::string_view data) override
        {
            switch (key_[0]) {
                case 'N':
                    target_.name = data;
                    reset_key();
                    // fmt::print(stderr, "NAME: {}\n", name_);
                    break;
                case 'C':
                    target_.continent = data;
                    reset_key();
                    break;
                case 'c':
                    target_.country = data;
                    reset_key();
                    break;
                case 'd':
                    target_.dates.emplace_back(data);
                    break;
                case 'l':
                    target_.lineage = data;
                    reset_key();
                    break;
                case 'v':
                    target_.virus_type = data;
                    reset_key();
                    break;
                default:
                    throw in_json::parse_error(AD_FORMAT("entry: unexpected key: \"{}\"", key_));
            }
            // }
            // else
            //     throw in_json::parse_error(AD_FORMAT("entry: unexpected key: \"{}\"", key_));
        }

      private:
        acmacs::seqdb::SeqdbEntry& target_;
    };

    class db : public in_json::stack_entry
    {
      public:
        db(std::vector<acmacs::seqdb::SeqdbEntry>& entries) : entries_{entries} {}

        const char* injson_name() override { return "db"; }

        std::unique_ptr<in_json::stack_entry> injson_put_object() override { return std::make_unique<entry>(entries_.emplace_back()); }

        void injson_put_string(std::string_view data) override
        {
            if (key_ == "  version") {
                if (data != "sequence-database-v2" && data != "sequence-database-v3")
                    throw in_json::parse_error(AD_FORMAT("unsupported version: {}", data));
                reset_key();
            }
            else if (key_ == "  date" || key_ == "_")
                reset_key();
            else
                throw in_json::parse_error(AD_FORMAT("unsupported field: \"{}\": {}", key_, data));
        }

        void injson_put_array() override {}
        void injson_pop_array() override { reset_key(); }

      private:
        std::vector<acmacs::seqdb::v3::SeqdbEntry>& entries_;
    };

    using sink = in_json::object_sink<std::vector<acmacs::seqdb::SeqdbEntry>, db>;

} // namespace local

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::parse(std::string_view source, std::vector<SeqdbEntry>& entries)
{
    local::sink sink{entries};
    in_json::parse(sink, std::begin(source), std::end(source));
    // fmt::print("INFO: seqdb entries read: {}\n", entries.size());

} // acmacs::seqdb::v3::parse

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
