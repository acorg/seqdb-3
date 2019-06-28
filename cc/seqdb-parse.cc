#include <stack>
#include <memory>

#include "acmacs-base/in-json.hh"
#include "seqdb-3/seqdb-parse.hh"
#include "seqdb-3/seqdb.hh"

// ----------------------------------------------------------------------

namespace local
{
    class error : public std::runtime_error
    {
      public:
        // template <typename S> error(S&& message) : std::runtime_error(fmt::format("error: {}", message)) {}
        error(std::string_view m1) : std::runtime_error{fmt::format("seqdb parsing error: {}", m1)} {}
        error(std::string_view m1, std::string_view m2) : error(fmt::format("{}{}", m1, m2)) {}
    };

    class stack_entry
    {
      public:
        stack_entry() = default;
        stack_entry(const stack_entry&) = default;
        virtual ~stack_entry() = default;
        virtual const char* injson_name() = 0;
        virtual std::unique_ptr<stack_entry> injson_put_object() = 0;
        virtual void injson_put_string(std::string_view data) { key_ = data; }
        virtual void injson_put_integer(std::string_view /*data*/) { throw error("stack_entry::injson_put_integer"); }
        // virtual void injson_put_real(double /*data*/) { throw error("stack_entry::injson_put_real"); }
        virtual void injson_put_array() { throw error("stack_entry::injson_put_array"); }
        virtual void injson_pop_array() { throw error("stack_entry::injson_pop_array"); }

      protected:
        std::string_view key_{};

        void reset_key() { key_ = std::string_view{}; }

    };

    class labs : public stack_entry
    {
      public:
        labs(seqdb::SeqdbSeq::labs_t& target) : target_{target} {}
        const char* injson_name() override { return "labs"; }
        std::unique_ptr<stack_entry> injson_put_object() override { throw error("labs: unexpected subobject"); }
        void injson_put_array() override {}
        void injson_pop_array() override {}

        void injson_put_string(std::string_view data) override
        {
            if (key_.empty()) {
                stack_entry::injson_put_string(data);
                target_.emplace_back(data, seqdb::SeqdbSeq::lab_ids_t{});
            }
            else
                target_.back().second.emplace_back(data);
        }

      private:
        seqdb::SeqdbSeq::labs_t& target_;
    };

    class seq : public stack_entry
    {
      public:
        seq(seqdb::SeqdbSeq& target) : target_{target} {}

        const char* injson_name() override { return "seq"; }

        std::unique_ptr<stack_entry> injson_put_object() override
        {
            // if (key_[0] == 'l') {
            reset_key();
            return std::make_unique<labs>(target_.lab_ids);
            // }
            // else
            //     throw error("seq: unexpected sub-object, key: ", key_);
        }

        void injson_put_array() override
        {
            // if (key_.size() != 1 || (key_[0] != 'p' && key_[0] != 'c' && key_[0] != 'h' && key_[0] != 'r'))
            //     throw error("seq: unexpected array, key: ", key_);
        }

        void injson_pop_array() override
        {
            // if (key_.size() != 1 || (key_[0] != 'p' && key_[0] != 'c' && key_[0] != 'h' && key_[0] != 'r'))
            //     throw error("seq: unexpected array, key: ", key_);
            reset_key();
        }

        void injson_put_string(std::string_view data) override
        {
            if (key_.empty())
                stack_entry::injson_put_string(data);
            else
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
                        target_.amino_acids = data;
                        reset_key();
                        break;
                    case 'n':
                        target_.nucs = data;
                        reset_key();
                        break;
                    default:
                        throw error("seq: unexpected key: ", data);
                }
            // }
            // else
            //     throw error("seq: unexpected key: ", data);
        }

        void injson_put_integer(std::string_view data) override
        {
            // if (key_.size() == 1) {
            switch (key_[0]) {
                case 's':
                    target_.a_shift = data;
                    reset_key();
                    break;
                case 't':
                    target_.n_shift = data;
                    reset_key();
                    break;
                default:
                    throw error("seq: unexpected integer, key: ", key_);
            }
            // }
            // else
            //     throw error("seq: unexpected integer, key: ", key_);
        }

      private:
        seqdb::SeqdbSeq& target_;
    };

    class entry : public stack_entry
    {
      public:
        entry(seqdb::SeqdbEntry& target) : target_{target} {}

        const char* injson_name() override { return "entry"; }

        std::unique_ptr<stack_entry> injson_put_object() override
        {
            return std::make_unique<seq>(target_.seqs.emplace_back()); // objects are unly under "s"
        }

        void injson_put_array() override
        {
            // fmt::print(stderr, "WARNING: entry::injson_put_array\n");
            // if (key_.size() != 1 || (key_[0] != 'd' && key_[0] != 's'))
            //     throw error("entry: unexpected array, key: " + std::string(key_));
        }

        void injson_pop_array() override
        {
            // fmt::print(stderr, "WARNING: entry::injson_pop_array\n");
            // if (key_.size() != 1 || (key_[0] != 'd' && key_[0] != 's'))
            //     throw error("entry: unexpected array, key: " + std::string(key_));
            reset_key();
        }

        void injson_put_string(std::string_view data) override
        {
            if (key_.empty())
                stack_entry::injson_put_string(data);
            else
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
                        throw error("entry: unexpected key: ", data);
                }
            // }
            // else
            //     throw error("entry: unexpected key: ", data);
        }

      private:
        seqdb::SeqdbEntry& target_;
    };

    class db : public stack_entry
    {
      public:
        db(std::vector<seqdb::SeqdbEntry>& entries) : entries_{entries} {}

        const char* injson_name() override { return "db"; }

        std::unique_ptr<stack_entry> injson_put_object() override { return std::make_unique<entry>(entries_.emplace_back()); }

        void injson_put_string(std::string_view data) override
        {
            if (key_ == "  version") {
                if (data != "sequence-database-v2")
                    throw error("unsupported version: ", data);
                reset_key();
            }
            else if (key_.empty())
                stack_entry::injson_put_string(data);
            else if (key_ == "  date" || key_ == "_")
                reset_key();
            else
                throw error("unsupported field: ", data);
        }

        void injson_put_array() override {}
        void injson_pop_array() override { reset_key(); }

      private:
        std::vector<seqdb::v3::SeqdbEntry>& entries_;
    };

    class sink
    {
      public:
        sink(std::vector<seqdb::SeqdbEntry>& entries) : entries_{entries} {}

        void injson_object_start()
        {
            if (target_.empty())
                target_.push(std::make_unique<db>(entries_));
            else
                target_.push(target_.top()->injson_put_object());
        }

        void injson_object_end() { target_.pop(); }
        void injson_array_start() { target_.top()->injson_put_array(); }
        void injson_array_end() { target_.top()->injson_pop_array(); }
        template <typename Iter> void injson_string(Iter first, Iter last) { target_.top()->injson_put_string({&*first, static_cast<size_t>(last - first)}); }
        template <typename Iter> void injson_integer(Iter first, Iter last) { target_.top()->injson_put_integer({&*first, static_cast<size_t>(last - first)}); }
        template <typename Iter> void injson_real(Iter /*first*/, Iter /*last*/) {}

      private:
        std::vector<seqdb::v3::SeqdbEntry>& entries_;
        std::stack<std::unique_ptr<stack_entry>> target_;
    };

} // namespace local

// ----------------------------------------------------------------------

void seqdb::v3::parse(std::string_view source, std::vector<SeqdbEntry>& entries)
{
    local::sink sink{entries};
    in_json::parse(sink, std::begin(source), std::end(source));
    fmt::print("INFO: seqdb entries read: {}\n", entries.size());

} // seqdb::v3::parse

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
