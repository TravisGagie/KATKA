// tokwrite: writes the run-length encoded tag array as 32-bit integers for bigrepair -i.
// usage: tokwrite prefix.tagruns both|genus|len out.int32
//   both: len_1, g_1, len_2, g_2, ... with lengths as dense ranks 0..L-1 and genera as L + g (disjoint values)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
static void die(const char *m) { fprintf(stderr, "tokwrite: %s\n", m); exit(1); }
struct Runs {   // streams the runs of a .tagruns file: run j has start rs[j], length rs[j+1]-rs[j] (last: n-rs[j]), genus tg[j]
    std::string f; uint64_t n = 0, r = 0;
    explicit Runs(const std::string &file) : f(file) {
        FILE *h = fopen(f.c_str(), "rb"); if (!h) die("cannot open the .tagruns file");
        if (fread(&n, 8, 1, h) != 1 || fread(&r, 8, 1, h) != 1) die("bad header"); fclose(h);
    }
    template <class F> void each(F fn) const {   // fn(j, length, genus)
        FILE *a = fopen(f.c_str(), "rb"), *b = fopen(f.c_str(), "rb");
        fseek(a, 16, SEEK_SET); fseek(b, 16 + 4 * r, SEEK_SET);
        const uint64_t B = 1 << 20; std::vector<uint32_t> s(B + 1); std::vector<uint16_t> g(B);
        uint32_t next; if (fread(&next, 4, 1, a) != 1) die("read");
        for (uint64_t j0 = 0; j0 < r; j0 += B) {
            uint64_t m = std::min(B, r - j0);
            s[0] = next;
            if (j0 + m < r) { if (fread(&s[1], 4, m, a) != m) die("read"); next = s[m]; }
            else { if (m > 1 && fread(&s[1], 4, m - 1, a) != m - 1) die("read"); s[m] = (uint32_t)n; }
            if (fread(g.data(), 2, m, b) != m) die("read");
            for (uint64_t i = 0; i < m; ++i) fn(j0 + i, (uint64_t)s[i + 1] - s[i], g[i]);
        }
        fclose(a); fclose(b);
    }
};

int main(int argc, char **argv) {
    if (argc != 4) die("usage: tokwrite prefix.tagruns both|genus|len out.int32");
    std::string mode = argv[2]; Runs R(argv[1]);
    std::vector<uint32_t> lens; uint32_t maxg = 0;
    R.each([&](uint64_t, uint64_t l, uint16_t g) { if (mode != "genus") lens.push_back((uint32_t)l); maxg = std::max<uint32_t>(maxg, g); });
    std::sort(lens.begin(), lens.end()); lens.erase(std::unique(lens.begin(), lens.end()), lens.end()); lens.shrink_to_fit();
    uint32_t L = lens.size();
    FILE *o = fopen(argv[3], "wb"); if (!o) die("cannot write");
    std::vector<uint32_t> buf; buf.reserve(1 << 20); uint64_t N = 0;
    auto put = [&](uint32_t x) { buf.push_back(x); ++N; if (buf.size() == buf.capacity()) { fwrite(buf.data(), 4, buf.size(), o); buf.clear(); } };
    R.each([&](uint64_t, uint64_t l, uint16_t g) {
        uint32_t lr = mode == "genus" ? 0 : (uint32_t)(std::lower_bound(lens.begin(), lens.end(), (uint32_t)l) - lens.begin());
        if (mode == "both") { put(lr); put(L + g); } else if (mode == "genus") put(g); else put(lr);
    });
    fwrite(buf.data(), 4, buf.size(), o); fclose(o);
    fprintf(stderr, "%s: %lu tokens, alphabet %u\n", argv[3], N, mode == "both" ? L + maxg + 1 : mode == "genus" ? maxg + 1 : L);
    return 0;
}
