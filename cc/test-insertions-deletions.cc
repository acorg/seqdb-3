#include "acmacs-base/fmt.hh"

// ----------------------------------------------------------------------

static void detect_deletions(std::string_view master, std::string_view to_align);

// ----------------------------------------------------------------------

static constexpr inline bool common(char a, char b) { return a == b && a != 'X' && a != '-'; }

static inline size_t number_of_common(std::string_view s1, std::string_view s2)
{
    size_t num = 0;
    for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
        if (common(*f1, *f2))
            ++num;
    }
    return num;
}

static inline size_t first_common(std::string_view s1, std::string_view s2)
{
    for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
        if (common(*f1, *f2))
            return static_cast<size_t>(f1 - s1.begin());
    }
    return std::string_view::npos;
}

static inline size_t find_head(std::string_view s1, std::string_view s2, size_t threshold)
{
    auto last_common = s1.end();
    for (auto f1 = s1.begin(), f2 = s2.begin(); f1 != s1.end() && f2 != s2.end(); ++f1, ++f2) {
        if (common(*f1, *f2))
            last_common = f1;
        else if (static_cast<size_t>(f1 - last_common) >= threshold)
            return static_cast<size_t>(last_common - s1.begin() + 1);
    }
    return std::min(s1.size(), s2.size());
}

static inline size_t find_tail(std::string_view s1, std::string_view s2, size_t threshold)
{
    auto last_common = s1.rend();
    for (auto f1 = s1.rbegin(), f2 = s2.rbegin(); f1 != s1.rend() && f2 != s2.rend(); ++f1, ++f2) {
        if (common(*f1, *f2))
            last_common = f1;
        else if (static_cast<size_t>(f1 - last_common) >= threshold) {
            fmt::print(stderr, "last_common {}\n", *last_common);
            return static_cast<size_t>(last_common - s1.rbegin() + 1);
        }
    }
    return std::min(s1.size(), s2.size());
}

// ----------------------------------------------------------------------

int main(int argc, char* const argv[])
{
    detect_deletions(argv[1], argv[2]);
}

// ----------------------------------------------------------------------

void detect_deletions(std::string_view master, std::string_view to_align)
{
    fmt::print(stderr, "{}\n{}\n\n", master, to_align);

    auto head = find_head(master, to_align, 3);
    auto tail = find_tail(master.substr(head), to_align.substr(head), 3);
    fmt::print(stderr, "head:{} tail:{} size:{}\n", head, tail, master.size());
    fmt::print(stderr, "{} {} {}\n{} {} {}\n",
               master.substr(0, head), master.substr(head, master.size() - head - tail), master.substr(master.size() - tail),
               to_align.substr(0, head), to_align.substr(head, to_align.size() - head - tail), to_align.substr(to_align.size() - tail));

    // fmt::print(stderr, "common all {} of {}\n", number_of_common(master, to_align), master.size());
    // auto chunk1 = number_of_common(master, to_align);
    // while (master[chunk1 - 1] != to_align[chunk1 - 1])
    //     --chunk1;
    // while (master[chunk1] == to_align[chunk1])
    //     ++chunk1;
    // const auto fc = first_common(master.substr(chunk1), to_align.substr(chunk1));
    // fmt::print(stderr, "first_common {}\n", fc);
    // if (fc < 4) {
    //     chunk1 += fc;
    //     while (master[chunk1] == to_align[chunk1])
    //         ++chunk1;
    // }
    // fmt::print(stderr, "common chunk1 left {} of {}\n", number_of_common(master.substr(0, chunk1), to_align.substr(0, chunk1)), chunk1);
    // fmt::print(stderr, "{}\n{}\n", master.substr(0, chunk1), to_align.substr(0, chunk1));
    // fmt::print(stderr, "common chunk1 right {}\n", number_of_common(master.substr(chunk1), to_align.substr(chunk1)));
    // fmt::print(stderr, "{}\n{}\n", master.substr(chunk1), to_align.substr(chunk1));

    // for (size_t ins = 1; ins < 10; ++ins) {
    //     fmt::print(stderr, "common chunk1+{} right {} of {}\n", ins, number_of_common(master.substr(chunk1 + ins), to_align.substr(chunk1)), master.size() - chunk1 - ins);
    //     fmt::print(stderr, "{}\n{}\n", master.substr(chunk1 + ins), to_align.substr(chunk1));
    // }

} // detect_deletions

// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
