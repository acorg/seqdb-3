#include "acmacs-base/read-file.hh"
#include "acmacs-base/counter.hh"
#include "acmacs-base/range-v3.hh"
#include "acmacs-base/to-json.hh"
#include "acmacs-chart-2/point-index-list.hh"
#include "seqdb-3/seqdb.hh"
#include "seqdb-3/log.hh"

// ----------------------------------------------------------------------

acmacs::seqdb::seq_id_t acmacs::seqdb::v3::ref::seq_id() const
{
    auto source = acmacs::string::join(acmacs::string::join_space, entry->name, seq().designation());
    if (entry->seqs.size() > 1 && seq_index > 0) {
        // there could be multiple seqs with the same designation, but seq_id must be unique, also garli does not like name duplicates
        std::vector<std::string> designations(entry->seqs.size());
        std::transform(std::begin(entry->seqs), std::end(entry->seqs), std::begin(designations), [](const auto& en) { return en.designation(); });
        if (std::count(std::begin(designations), std::end(designations), designations[seq_index]) > 1)
            source.append(fmt::format("_d{}", seq_index));
    }
    return make_seq_id(source);

} // acmacs::seqdb::v3::ref::seq_id

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::multiple_dates(bool do_filter)
{
    if (do_filter)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.entry->dates.size() < 2; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::multiple_dates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::subtype(const acmacs::uppercase& virus_type)
{
    if (!virus_type.empty()) {
        std::string_view vt = virus_type;
        if (vt == "H1")
            vt = "A(H1N1)";
        else if (vt == "H3")
            vt = "A(H3N2)";
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [vt](const auto& en) { return en.entry->virus_type != vt; }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::subtype

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::lineage(const acmacs::uppercase& lineage)
{
    if (!lineage.empty()) {
        std::string_view lin = lineage;
        switch (lin[0]) {
          case 'V': lin = "VICTORIA"; break;
          case 'Y': lin = "YAMAGATA"; break;
        }
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lin](const auto& en) { return en.entry->lineage != lin; }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::lineage

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::lab(const acmacs::uppercase& lab)
{
    if (!lab.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [lab=static_cast<std::string_view>(lab)](const auto& en) { return !en.has_lab(lab); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::lab

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::whocc_lab(bool do_filter)
{
    if (do_filter)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return ! (en.has_lab("CDC") || en.has_lab("CRICK") || en.has_lab("NIID") || en.has_lab("VIDRL")); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::whocc_lab

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::host(const acmacs::uppercase& host)
{
    if (!host.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [host=static_cast<std::string_view>(host)](const auto& en) { return en.entry->host() != host; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::host

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::continent(const acmacs::uppercase& continent)
{
    if (!continent.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [continent=static_cast<std::string_view>(continent)](const auto& en) { return en.entry->continent != continent; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::continent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::country(const acmacs::uppercase& country)
{
    if (!country.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [country=static_cast<std::string_view>(country)](const auto& en) { return en.entry->country != country; }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::country

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::with_issues(bool keep_with_issues)
{
    if (!keep_with_issues)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return en.has_issues(); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::with_issues

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::clade(const Seqdb& seqdb, const acmacs::uppercase& clade)
{
    if (!clade.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [&seqdb,clade=static_cast<std::string_view>(clade)](const auto& en) { return !en.has_clade(seqdb, clade); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::clade

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent(size_t recent, master_only master)
{
    if (recent > 0) {
        if (master == master_only::yes)
            keep_master_only();
        if (refs_.size() > recent) {
            sort_by_date_recent_first();
            refs_.erase(std::next(std::begin(refs_), static_cast<ssize_t>(recent)), std::end(refs_));
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::recent_matched(const std::vector<size_t>& recent_matched, master_only master)
{
    if (recent_matched.size() > 1 && refs_.size() > recent_matched[0]) {
        if (recent_matched.size() != 2)
            throw std::runtime_error{fmt::format("invalid recent-matched specification: {} {}", recent_matched, recent_matched.size())};
        if (master == master_only::yes)
            keep_master_only();
        if ((recent_matched[0] + recent_matched[1]) < refs_.size()) {
            sort_by_date_recent_first();
            if (master == master_only::yes) {
                // if ref (master) has no hi names and one of its slaves has hi name, replace ref with slave that has hi names and return false
                // if ref (master) has no hi names and none of its slaves has hi name, return true (to remove from refs_)
                size_t number_to_keep = recent_matched[1];
                const auto without_hi_names = [&number_to_keep, remove = true, keep = false](auto& ref) {
                    if (number_to_keep == 0)
                        return remove;
                    if (ref.has_hi_names()) {
                        --number_to_keep;
                        return keep;
                    }
                    const auto& slaves = ref.seq().slaves();
                    if (const auto slave_to_use = std::find_if(std::begin(slaves), std::end(slaves), [](const auto& slave) { return slave.has_hi_names(); }); slave_to_use != std::end(slaves)) {
                        //     ref = *slave_to_use;
                        --number_to_keep;
                        return keep;
                    }
                    else
                        return remove;
                };

                const auto end = std::remove_if(std::next(std::begin(refs_), static_cast<ssize_t>(recent_matched[0])), std::end(refs_), without_hi_names);
                refs_.erase(end, std::end(refs_));
            }
            else {
                const auto usable_size =
                    std::remove_if(std::next(std::begin(refs_), static_cast<ssize_t>(recent_matched[0])), std::end(refs_), [](const auto& en) { return !en.has_hi_names(); }) - std::begin(refs_);
                refs_.erase(std::next(std::begin(refs_), std::min(usable_size, static_cast<ssize_t>(recent_matched[0] + recent_matched[1]))), std::end(refs_));
            }
        }
    }
    return *this;

} // acmacs::seqdb::v3::subset::recent_matched

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::keep_master_only()
{
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.is_master(); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::keep_master_only

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::subset::remove(ref_indexes& to_remove)
{
    std::sort(std::begin(to_remove), std::end(to_remove));
    const auto rm_end = std::unique(std::begin(to_remove), std::end(to_remove));
    auto rm_iter = std::begin(to_remove);
    size_t current_index{0};
    const auto remove_predicate = [&current_index,&rm_iter,rm_end](const auto&) {
        if (rm_iter != rm_end && *rm_iter == current_index++) {
            ++rm_iter;
            return true;
        }
        else
            return false;
    };
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), remove_predicate), std::end(refs_));

} // acmacs::seqdb::v3::subset::remove

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::subset::keep(ref_indexes& to_keep)
{
    std::sort(std::begin(to_keep), std::end(to_keep));
    const auto keep_end = std::unique(std::begin(to_keep), std::end(to_keep));
    auto keep_iter = std::begin(to_keep);
    size_t current_index{0};
    const auto remove_predicate = [&current_index,&keep_iter,keep_end](const auto&) {
        if (keep_iter != keep_end && *keep_iter == current_index++) {
            ++keep_iter;
            return false;
        }
        else
            return true;
    };
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), remove_predicate), std::end(refs_));

} // acmacs::seqdb::v3::subset::keep

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::random(size_t random)
{
    if (random > 0 && refs_.size() > random) {
        std::mt19937 generator{std::random_device()()};
        std::uniform_int_distribution<size_t> distribution(0, refs_.size() - 1);
        ref_indexes to_keep(random);
        std::generate_n(to_keep.begin(), random, [&distribution,&generator]() { return distribution(generator); });
        keep(to_keep);
    }
    return *this;

} // acmacs::seqdb::v3::subset::random

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::subset_every_month(double /*fraction*/)
{
    AD_ERROR("subset::subset_every_month not implemented");
    return *this;

    // if (fraction > 1.0 || (fraction < 0.0 && !float_equal(fraction, -1.0)))
    //     AD_WARNING("unrecognized subset_every_month fraction value: {} (not subsetting)", fraction);
    // if (fraction >= 0.0 && fraction <= 1.0) {
    //     sort_by_date_oldest_first();
    //     ref_indexes to_keep = range_from_0_to(size()) | ranges::to_vector;
    //     std::mt19937 generator{std::random_device()()};
    //     size_t ind_start{0};
    //     std::string_view year_month = refs_[to_keep[ind_start]].entry->date().substr(0, 7);
    //     for (const size_t ind : range_from_to(1ul, size())) {
    //     // std::uniform_int_distribution<size_t> distribution(0, refs_.size() - 1);
    //     }
    // }
    // return *this;

} // acmacs::seqdb::v3::subset::subset_every_month

// ----------------------------------------------------------------------

void acmacs::seqdb::v3::remove_nuc_duplicates(subset::refs_t& refs, bool keep_hi_matched)
{
    if (keep_hi_matched) {
        // master sequences and hi matched (if requested) in the [std::begin(refs), to_remove_candidates_start] range
        const auto to_remove_canditates_start =
            std::partition(std::begin(refs), std::end(refs), [](const auto& ref) { return ref.is_master() || ref.is_hi_matched(); });
        // move slave seq from [to_remove_canditates_start, std::end(refs)] that reference to
        // a sequence in [std::begin(refs), to_remove_candidates_start]
        // to the [to_remove_start, std::end(refs)] range
        const auto to_remove_start = std::partition(to_remove_canditates_start, std::end(refs), [beg = std::begin(refs), end = to_remove_canditates_start](const auto& ref1) {
            return std::find_if(beg, end, [&ref1](const auto& ref2) { return ref2.matches(ref1.seq().master); }) == end;
        });

        refs.erase(to_remove_start, std::end(refs));
    }
    else {
        refs.erase(std::remove_if(std::begin(refs), std::end(refs), [](const auto& ref) { return !ref.is_master(); }), std::end(refs));
    }

} // acmacs::seqdb::v3::remove_nuc_duplicates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_nuc_duplicates(bool do_remove, bool keep_hi_matched)
{
    if (do_remove)
        acmacs::seqdb::v3::remove_nuc_duplicates(refs_, keep_hi_matched);
    return *this;

} // acmacs::seqdb::v3::subset::remove_nuc_duplicates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_empty(const Seqdb& seqdb, bool nuc)
{
    const auto is_empty = [&seqdb, nuc](const auto& ref) {
        const auto& seq = ref.seq_with_sequence(seqdb);
        // AD_LOG(acmacs::log::sequences, "      master aa:{} nuc:{} orig:{}", seq.aa_aligned_length_master(), seq.nuc_aligned_length_master(), ref.seq_id());
        return nuc ? seq.nuc_aligned_length_master() == 0 : seq.aa_aligned_length_master() == 0;
    };

    AD_LOG(acmacs::log::sequences, "removing empty ({}) from {} sequences", nuc ? "nuc" : "aa", refs_.size());
    refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), is_empty), std::end(refs_));
    AD_LOG(acmacs::log::sequences, "    {} sequences left", refs_.size());
    return *this;

} // acmacs::seqdb::v3::subset::remove_empty

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_nuc_duplicates(const Seqdb& seqdb, bool do_remove, bool keep_hi_matched)
// {
//     if (do_remove) {
//         const acmacs::Counter counter_nuc_length(refs_, [&seqdb](const auto& en) { return en.nuc_aligned_length(seqdb); });
//         const auto nuc = [&seqdb, nuc_common_length = counter_nuc_length.max().first](const auto& en) { return en.nuc_aligned(seqdb, nuc_common_length); };
//         const auto hi_names = [](const auto& en) { return en.seq().hi_names.size(); };
//         std::sort(std::begin(refs_), std::end(refs_), [=](const auto& e1, const auto& e2) {
//             if (const auto n1 = nuc(e1), n2 = nuc(e2); n1 == n2)
//                 return hi_names(e1) > hi_names(e2);
//             else
//                 return n1 < n2;
//         });
//         if (keep_hi_matched) {
//             refs_.erase(std::unique(std::begin(refs_), std::end(refs_), [=](const auto& e1, const auto& e2) { return nuc(e1) == nuc(e2) && (hi_names(e1) == 0 || hi_names(e2) == 0); }),
//                         std::end(refs_));
//         }
//         else {
//             refs_.erase(std::unique(std::begin(refs_), std::end(refs_), [=](const auto& e1, const auto& e2) { return nuc(e1) == nuc(e2); }), std::end(refs_));
//         }
//     }
//     return *this;

// } // acmacs::seqdb::v3::subset::remove_nuc_duplicates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::refs_t::iterator acmacs::seqdb::v3::subset::most_recent_with_hi_name()
{
    auto result = std::end(refs_);
    std::string_view date;
    for (auto refp = std::begin(refs_); refp != std::end(refs_); ++refp) {
        if (refp->has_hi_names() && refp->entry->date() > date) { // refs_[no].seq().reassortants.empty() &&
            result = refp;
            date = refp->entry->date();
            // fmt::print(stderr, "DEBUG: [{}] {}\n", date, result->full_name());
        }
    }
    return result;

} // acmacs::seqdb::v3::subset::most_recent_with_hi_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::with_hi_name(bool with_hi_name)
{
    if (with_hi_name)
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [](const auto& en) { return !en.has_hi_names(); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::with_hi_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::aa_at_pos(const Seqdb& seqdb, const amino_acid_at_pos1_eq_list_t& aa_at_pos)
{
    if (!aa_at_pos.empty()) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&aa_at_pos,&seqdb](const auto& en) {
                                       try {
                                           const auto& seq = en.seq().with_sequence(seqdb);
                                           return seq.amino_acids.empty() || !seq.matches(aa_at_pos); // true to remove
                                       }
                                       catch (std::exception& err) {
                                           throw std::runtime_error{fmt::format("{}, full_name: {}", err, en.full_name())};
                                       }
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::aa_at_pos

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::nuc_at_pos(const Seqdb& seqdb, const nucleotide_at_pos1_eq_list_t& nuc_at_pos)
{
    if (!nuc_at_pos.empty()) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&nuc_at_pos,&seqdb](const auto& en) {
                                       try {
                                           const auto& seq = en.seq().with_sequence(seqdb);
                                           return seq.nucs.empty() || !seq.matches(nuc_at_pos); // true to remove
                                       }
                                       catch (std::exception& err) {
                                           throw std::runtime_error{fmt::format("{}, full_name: {}", err, en.full_name())};
                                       }
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::nuc_at_pos

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::min_aa_length(const Seqdb& seqdb, size_t length)
{
    if (length) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [length, &seqdb](const auto& en) { return en.aa_aligned_length(seqdb) < length; }), std::end(refs_));
    }
    return *this;


} // acmacs::seqdb::v3::subset::min_aa_length

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::min_nuc_length(const Seqdb& seqdb, size_t length)
{
    if (length) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [length, &seqdb](const auto& en) { return en.nuc_aligned_length(seqdb) < length; }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::min_nuc_length

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::remove_with_front_back_deletions(const Seqdb& seqdb, bool remove, size_t length)
{
    if (remove) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [length, &seqdb](const auto& en) {
            const auto nucs = en.nuc_aligned(seqdb);
            if (nucs.at(pos1_t{1}) == '-')
                return true;
            if (length > 0 && (nucs.size() < pos0_t{length} || nucs.at(pos1_t{length}) == '-'))
                return true;    // too short or has deletion in the last nuc
            return false;
        }), std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::remove_with_front_back_deletions

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::names_matching_regex(const std::vector<std::string_view>& regex_list)
{
    if (!regex_list.empty()) {
        std::vector<std::regex> re_list(regex_list.size());
        std::transform(std::begin(regex_list), std::end(regex_list), std::begin(re_list),
                       [](const auto& regex_s) { return std::regex(std::begin(regex_s), std::end(regex_s), std::regex_constants::icase); });
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&re_list](const auto& en) {
                                       return std::none_of(std::begin(re_list), std::end(re_list), [full_name = en.full_name()](const auto& re) { return std::regex_search(full_name, re); });
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::names_matching_regex

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::exclude(const std::vector<std::string_view>& seq_ids)
{
    if (!seq_ids.empty()) {
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_),
                                   [&seq_ids](const auto& en) {
                                       return std::any_of(std::begin(seq_ids), std::end(seq_ids), [seq_id = en.seq_id()](const auto& si) { return si == seq_id; });
                                   }),
                    std::end(refs_));
    }
    return *this;

} // acmacs::seqdb::v3::subset::exclude

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::dates(std::string_view start, std::string_view end)
{
    if (!start.empty() || !end.empty())
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), [start,end](const auto& en) { return !en.entry->date_within(start, end); }), std::end(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::dates

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend(std::string_view seq_id, const Seqdb& seqdb)
{
    if (!seq_id.empty()) {
        auto candidates = seqdb.select_by_seq_id(seq_id);
        if (candidates.empty())
            throw std::runtime_error{fmt::format("no sequences with seq-id \"{}\" found (seqdb::v3::subset::prepend)", seq_id)};
        refs_.erase(std::remove(std::begin(refs_), std::end(refs_), candidates.front()), std::end(refs_)); // remove it, if selected earlier
        refs_.insert(std::begin(refs_), candidates.front());
    }
    return *this;

} // prepend

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend(const std::vector<std::string_view>& seq_ids, const Seqdb& seqdb)
{
    if (!seq_ids.empty()) {
        auto candidates = seqdb.select_by_seq_id(seq_ids);
        if (candidates.empty())
            throw std::runtime_error{fmt::format("no sequences by seq-ids found to prepend")};
        const auto select_to_remove = [&candidates](const auto& ref) { return std::find(std::begin(candidates), std::end(candidates), ref) != std::end(candidates); };
        refs_.erase(std::remove_if(std::begin(refs_), std::end(refs_), select_to_remove), std::end(refs_)); // remove it, if selected earlier
        refs_.insert(std::begin(refs_), std::begin(candidates), std::end(candidates));
    }
    return *this;

} // prepend

// ----------------------------------------------------------------------

// acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::prepend_single_matching(std::string_view re, const Seqdb& seqdb)
// {
//     if (!re.empty()) {
//         auto candidates = seqdb.select_by_regex(re);
//         if (candidates.size() != 1) {
//             fmt::print(stderr, "WARNING: Selected sequences: {}\n", candidates.size());
//             for (const auto& candidate : candidates)
//                 fmt::print(stderr, "    {}\n", candidate.seq_id());
//             throw std::runtime_error{fmt::format("regular expression must select single sequence: \"{}\", selected: {}", re, candidates.size())};
//         }
//         refs_.erase(std::remove(std::begin(refs_), std::end(refs_), candidates.front()), std::end(refs_)); // remove it, if selected earlier
//         refs_.insert(std::begin(refs_), candidates.front());
//     }
//     return *this;

// } // acmacs::seqdb::v3::subset::prepend_single_matching

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::sort(sorting srt)
{
    switch (srt) {
        case sorting::none:
            break;
        case sorting::name_asc:
            sort_by_name_asc();
            break;
        case sorting::name_desc:
            sort_by_name_desc();
            break;
        case sorting::date_asc:
            sort_by_date_oldest_first();
            break;
        case sorting::date_desc:
            sort_by_date_recent_first();
            break;
    }
    return *this;

} // acmacs::seqdb::v3::subset::sort

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_stat(const Seqdb& seqdb, bool do_report)
{
    if (do_report) {
        if (!refs_.empty()) {
            size_t with_hi_names = 0;
            std::string_view min_date = refs_.front().entry->date(), max_date = min_date;
            Counter<std::string> by_year;
            Counter<size_t> aa_length, nuc_length;
            for (const auto& ref : refs_) {
                const auto date = ref.entry->date();
                if (date < min_date)
                    min_date = date;
                else if (date > max_date)
                    max_date = date;
                if (date.size() >= 4)
                    by_year.count(date.substr(0, 4));
                if (!ref.seq().hi_names.empty())
                    ++with_hi_names;
                aa_length.count(ref.seq_with_sequence(seqdb).aa_aligned_length_master());
                nuc_length.count(ref.seq_with_sequence(seqdb).nuc_aligned_length_master());
            }
            fmt::print(stderr, "Selected sequences: {:6d}\n      HiDb matches: {:6d}\n        Date range: {} - {}\n", refs_.size(), with_hi_names, min_date, max_date);
            constexpr const size_t limit{10};
            fmt::print(stderr, "AA length:\n{}    {:4d} more lengths\nNucleotide lengths:\n{}    {:4d} more lengths\nBy year:\n{}",                                               //
                       aa_length.report_sorted_max_first("    {value:4d}  {counter:6d}  {counter_percent:3.0f}%\n", limit), aa_length.size() > limit ? aa_length.size() - limit : 0, //
                       nuc_length.report_sorted_max_first("    {value:4d}  {counter:6d}  {counter_percent:3.0f}%\n", limit),
                       nuc_length.size() > limit ? nuc_length.size() - limit : 0, //
                       by_year.report("    {value}  {counter:6d}  {counter_percent:3.0f}%\n"));
        }
        else {
            fmt::print(stderr, "No sequences selected\n");
        }
    }

    return *this;

} // acmacs::seqdb::v3::subset::report_stat

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_stat_month_region(bool do_report)
{
    if (do_report) {
        if (!refs_.empty()) {

            using namespace std::string_view_literals;
            constexpr static const std::array continents{"AFRICA"sv, "NORTH-AMERICA"sv, "CENTRAL-AMERICA"sv, "SOUTH-AMERICA"sv, "ASIA"sv, "AUSTRALIA-OCEANIA"sv, "MIDDLE-EAST"sv, "EUROPE"sv, "RUSSIA"sv, "UNKNOWN"sv};
            struct MonthEntry
            {
                size_t total{0};
                std::array<size_t, continents.size()> per_region{0};
            };

            std::map<std::string, MonthEntry> stat;
            for (const auto& ref : refs_) {
                std::string date{ref.entry->date()};
                if (date.size() > 7)
                    date.resize(7);
                else if (date.size() == 4)
                    date += "-??";
                auto& en = stat[date];
                ++en.total;
                if (const auto continent = std::find(std::begin(continents), std::end(continents), ref.entry->continent); continent != std::end(continents)) {
                    ++en.per_region[static_cast<size_t>(continent - std::begin(continents))];
                }
                else {
                    if (!ref.entry->continent.empty())
                        AD_WARNING("Continent name not found: \"{}\"", ref.entry->continent);
                    ++en.per_region.back();
                }
            }

            fmt::print("             Africa   N.America C.America S.America   Asia     Oceania  Mid.East   Europe    Russia   Unknown    TOTAL\n");
            for (const auto& [date, data] : stat) {
                fmt::print("{}  ", date);
                for (size_t ind{0}; ind < data.per_region.size(); ++ind) {
                    if (data.per_region[ind])
                        fmt::print("  {:6d}  ", data.per_region[ind]);
                    else
                        fmt::print("          ");
                }
                fmt::print("  {:6d}\n", data.total);
            }
        }
        else {
            fmt::print(stderr, "No sequences selected\n");
        }
    }

    return *this;

} // acmacs::seqdb::v3::subset::report_stat_month_region

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::report_aa_at(const Seqdb& seqdb, const pos1_list_t& pos1_list)
{
    if (!pos1_list.empty() && !refs_.empty()) {
        std::vector<CounterChar> counters(pos1_list.size());
        for (const auto& ref : refs_) {
            for (auto index : acmacs::range(pos1_list.size()))
                counters[index].count(ref.aa_at_pos(seqdb, pos1_list[index]));
        }
        fmt::print(stderr, "AA at pos stat:\n");
        for (auto index : acmacs::range(pos1_list.size()))
            fmt::print(stderr, "  {}\n{}", pos1_list[index], counters[index].report_sorted_max_first(fmt::format("    {:3d}{{first}}  {{second:5d}}\n", pos1_list[index])));
    }
    return *this;

} // acmacs::seqdb::v3::subset::report_aa_at

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::export_sequences(std::string_view filename, const Seqdb& seqdb, const export_options& options)
{
    if (!filename.empty()) {
        auto to_export = export_collect(seqdb, options);

        if (options.e_most_common_length == export_options::most_common_length::yes) {
            const acmacs::Counter counter(to_export, [](const auto& en) { return en.sequence.size(); });
            const auto most_common_length = counter.max().first;
            AD_LOG(acmacs::log::fasta, "most common length: {}", most_common_length);
            ranges::for_each(to_export, [most_common_length](auto& en) { en.sequence.resize(most_common_length, '-'); });
        }
        else if (options.e_length > 0) {
            AD_LOG(acmacs::log::fasta, "sequence length for exporting: {}", options.e_length);
            ranges::for_each(to_export, [length = options.e_length](auto& en) { en.sequence.resize(length, '-'); });
        }

        ranges::for_each(to_export, [deletion_report_threshold = options.e_deletion_report_threshold](auto& en) {
            const auto dels = static_cast<size_t>(ranges::count_if(en.sequence, [](char nuc_aa) { return nuc_aa == '-' || nuc_aa == 'X'; }));
            const auto dels_at_the_end = en.sequence.back() == '-' || en.sequence.back() == 'X';
            if (dels_at_the_end || dels > deletion_report_threshold)
                AD_WARNING("{}: {} deletions or unknown AAs or deletions at the end", en.seq_id, dels);
        });

        AD_LOG(acmacs::log::fasta, "writing {} sequences to {}", to_export.size(), filename);
        acmacs::file::write(filename, export_fasta(to_export, options));
    }
    return *this;

} // acmacs::seqdb::v3::subset::export_sequences

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::export_json_sequences(std::string_view filename, const Seqdb& seqdb, const export_options& options)
{
    if (!filename.empty()) {
        auto to_export = export_collect(seqdb, options);

        if (options.e_most_common_length == export_options::most_common_length::yes) {
            const acmacs::Counter counter(to_export, [](const auto& en) { return en.sequence.size(); });
            const auto most_common_length = counter.max().first;
            AD_LOG(acmacs::log::fasta, "most common length: {}", most_common_length);
            ranges::for_each(to_export, [most_common_length](auto& en) { en.sequence.resize(most_common_length, '-'); });
        }
        else if (options.e_length > 0) {
            AD_LOG(acmacs::log::fasta, "sequence length for exporting: {}", options.e_length);
            ranges::for_each(to_export, [length = options.e_length](auto& en) { en.sequence.resize(length, '-'); });
        }

        AD_LOG(acmacs::log::fasta, "writing {} sequences to {}", to_export.size(), filename);
        acmacs::file::write(filename, export_json(to_export, options));
    }
    return *this;

} // acmacs::seqdb::v3::subset::export_json_sequences

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::make_name(const Seqdb& seqdb, std::string_view name_format, const ref& entry) const
{
    const auto nf = ::string::replace(::string::replace(name_format, "\\t", "\t"), "\\n", "\n");
    return fmt::format(nf,
                       fmt::arg("seq_id", entry.seq_id()),
                       fmt::arg("full_name", entry.full_name()),
                       fmt::arg("hi_name_or_full_name", entry.hi_name_or_full_name()),
                       fmt::arg("hi_names", entry.seq().hi_names),
                       fmt::arg("hi_name", !entry.seq().hi_names.empty() ? entry.seq().hi_names.front() : std::string_view{}),
                       fmt::arg("lineage", entry.entry->lineage),
                       fmt::arg("name", entry.entry->name),
                       fmt::arg("date", entry.entry->date()),
                       fmt::arg("dates", entry.entry->dates),
                       fmt::arg("lab_id", entry.seq().lab_id()),
                       fmt::arg("passage", entry.seq().passage()),
                       fmt::arg("clades", entry.seq().with_sequence(seqdb).clades),
                       fmt::arg("lab", entry.seq().lab()),
                       fmt::arg("country", entry.entry->country),
                       fmt::arg("continent", entry.entry->continent),
                       fmt::arg("group_no", entry.group_no ? fmt::format("group:{}", entry.group_no) : std::string{}),
                       fmt::arg("hamming_distance", entry.hamming_distance),
                       fmt::arg("nuc_length", entry.seq().nuc_aligned_length_master()),
                       fmt::arg("aa_length", entry.seq().aa_aligned_length_master()),
                       fmt::arg("gisaid_accession_numbers", acmacs::string::join(acmacs::string::join_sep_t{"|"}, entry.seq().gisaid.isolate_ids)),
                       fmt::arg("ncbi_accession_numbers", acmacs::string::join(acmacs::string::join_sep_t{"|"}, entry.seq().gisaid.sample_ids_by_sample_provider))
                       );

} // acmacs::seqdb::v3::subset::make_name

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset::collected_t acmacs::seqdb::v3::subset::export_collect(const Seqdb& seqdb, const export_options& options) const
{
    const auto get_seq = [&options,&seqdb](const auto& entry) -> std::string_view {
        const auto& seq = entry.seq().with_sequence(seqdb);
        AD_LOG(acmacs::log::fasta, "{} has-seq:{}", entry.seq_id(), entry.is_master());
        if (!entry.is_master())
            AD_LOG(acmacs::log::fasta, "    ref:({} {})", entry.seq().master.name, entry.seq().master.hash);
        AD_LOG(acmacs::log::fasta, "    aa:{} nuc:{}", seq.aa_aligned_length_master(), seq.nuc_aligned_length_master());
        if (options.e_format == export_options::format::fasta_aa) {
            if (options.e_aligned == export_options::aligned::yes)
                return *seq.aa_aligned_master();
            else
                return std::get<std::string_view>(seq.amino_acids);
        }
        else {
            if (options.e_aligned == export_options::aligned::yes)
                return *seq.nuc_aligned_master();
            else
                return std::get<std::string_view>(seq.nucs);
        }
    };

    collected_t result(refs_.size()); // {seq_id, sequence}
    std::transform(std::begin(refs_), std::end(refs_), std::begin(result),
                   [this, &options, &get_seq, &seqdb](const auto& en) -> collected_entry_t { return {make_name(seqdb, options.e_name_format, en), std::string{get_seq(en)}}; });
    // remove entries with empty sequences
    result.erase(std::remove_if(std::begin(result), std::end(result), [](const auto& en) { return en.sequence.empty(); }), std::end(result));
    AD_LOG(acmacs::log::fasta, "collected for exporting: {}", result.size());
    return result;

} // acmacs::seqdb::v3::subset::export_collect

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::export_fasta(const collected_t& entries, const export_options& options)
{
    std::string output;
    const auto output_size =
        std::accumulate(std::begin(entries), std::end(entries), 0UL, [](size_t size, const auto& en) { return size + en.seq_id.size() + en.sequence.size() + 2 + en.sequence.size() / 40; });
    output.reserve(output_size);
    for (const auto& en : entries) {
        output.append(1, '>');
        output.append(en.seq_id);
        output.append(1, '\n');
        if (options.e_wrap_at == 0 || options.e_wrap_at >= en.sequence.size()) {
            output.append(en.sequence);
            output.append(1, '\n');
        }
        else {
            for (const auto chunk : en.sequence | ranges::views::chunk(options.e_wrap_at)) {
                output.append(ranges::to<std::string>(chunk));
                output.append(1, '\n');
            }
        }
    }
    fmt::print("INFO: exported to fasta: {}\n", entries.size());
    return output;

} // acmacs::seqdb::v3::subset::export_fasta

// ----------------------------------------------------------------------

std::string acmacs::seqdb::v3::subset::export_json(const collected_t& entries, const export_options& options)
{
    to_json::array arr;
    for (const auto& en : entries) {
        arr << to_json::object{
            to_json::key_val{"N", en.seq_id},
            to_json::key_val{"S", en.sequence},
        };
    }
    return fmt::format("{}\n", arr);

} // acmacs::seqdb::v3::subset::export_json

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::append(const subset& another)
{
    std::copy(std::begin(another), std::end(another), std::back_inserter(refs_));
    return *this;

} // acmacs::seqdb::v3::subset::append

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset acmacs::seqdb::v3::subset::filter_by_indexes(const acmacs::chart::PointIndexList& indexes, enum matched_only matched_only) const
{
    subset result;
    for (auto index : indexes) {
        if (index < refs_.size() && (matched_only == matched_only::no || refs_[index]))
            result.refs_.push_back(refs_[index]);
    }
    return result;

} // acmacs::seqdb::v3::subset::filter_by_indexes

// ----------------------------------------------------------------------

acmacs::seqdb::v3::subset& acmacs::seqdb::v3::subset::print(const Seqdb& seqdb, std::string_view name_format, bool do_print)
{
    if (do_print) {
        for (const auto& ref : *this)
            fmt::print("{}\n", make_name(seqdb, name_format, ref));
    }
    return *this;

} // acmacs::seqdb::v3::subset::print

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
