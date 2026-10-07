// tagz: size of the LZ77 parse of the run-length encoded tag array.
// usage: tagz prefix.tagruns [both|genus|len]
//   both  (default): the sequence len_1, g_1, len_2, g_2, ... of run lengths and genus IDs, one token each
//                    (lengths and genus IDs come from disjoint alphabets, so a phrase cannot match one with the other)
//   genus:           the genus IDs of the runs only
//   len:             the run lengths only
// Prints the number of tokens N, the alphabet size, the number z of LZ77 phrases (greedy, self-referential, each
// phrase the longest previous match or one new token) and N / z.  Suffix array with libsais; LZ77 with the KKP3
// algorithm (previous/next smaller values of the suffix array).  Memory: about 12N bytes at the peak.
// .tagruns format (rz-tagbuild stage 1): u64 n, u64 r, u32 run_start[r], u16 genus[r].
#include "libsais.h"
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static void die(const char *m) { fprintf(stderr, "tagz: %s\n", m); exit(1); }

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
    if (argc < 2) die("usage: tagz prefix.tagruns [both|genus|len]");
    std::string mode = argc > 2 ? argv[2] : "both";
    if (mode != "both" && mode != "genus" && mode != "len") die("mode must be both, genus or len");
    double t0 = now();
    Runs R(argv[1]);
    uint64_t r = R.r, N = mode == "both" ? 2 * r : r;
    if (N >= (1ull << 31) - 16) die("too many tokens for 32-bit libsais");
    // dense ranks of the run lengths, and the number of genera
    std::vector<uint32_t> lens; uint32_t maxg = 0;
    if (mode != "genus") { lens.reserve(r); }
    R.each([&](uint64_t, uint64_t l, uint16_t g) { if (mode != "genus") lens.push_back((uint32_t)l); maxg = std::max<uint32_t>(maxg, g); });
    std::sort(lens.begin(), lens.end()); lens.erase(std::unique(lens.begin(), lens.end()), lens.end()); lens.shrink_to_fit();
    uint64_t L = lens.size(), G = maxg + 1;
    auto lrank = [&](uint64_t l) { return (int32_t)(std::lower_bound(lens.begin(), lens.end(), (uint32_t)l) - lens.begin()); };
    auto build = [&](int32_t *T) {
        R.each([&](uint64_t j, uint64_t l, uint16_t g) {
            if (mode == "both") { T[2 * j] = lrank(l); T[2 * j + 1] = (int32_t)(L + g); }
            else if (mode == "genus") T[j] = g;
            else T[j] = lrank(l);
        });
    };
    int64_t k = mode == "both" ? L + G : mode == "genus" ? G : L;
    fprintf(stderr, "n=%lu runs=%lu tokens=%lu alphabet=%ld (lengths %lu, genera %lu)\n", R.n, r, N, k, L, G);
    // suffix array
    const int32_t fs = 6 * (int32_t)std::min<int64_t>(k, 1 << 24);
    int32_t *T = (int32_t *)malloc(4 * N), *SA = (int32_t *)malloc(4 * (N + fs));
    if (!T || !SA) die("out of memory (T, SA)");
    build(T);
    if (libsais_int(T, SA, (int32_t)N, (int32_t)k, fs) != 0) die("libsais_int failed");
    fprintf(stderr, "suffix array: %.0f s\n", now() - t0);
    free(T);
    // KKP3: previous and next smaller values (text positions) of each suffix in suffix-array order
    int32_t *PSV = (int32_t *)malloc(4 * N), *NSV = (int32_t *)malloc(4 * N);
    if (!PSV || !NSV) die("out of memory (PSV, NSV)");
    int32_t top = -1;
    for (uint64_t i = 0; i <= N; ++i) {
        int32_t x = i < N ? SA[i] : -1;
        while (top > x) { NSV[top] = x; top = PSV[top]; }
        if (i < N) { PSV[x] = top; top = x; }
    }
    free(SA);
    T = (int32_t *)malloc(4 * N); if (!T) die("out of memory (T again)");
    build(T);
    fprintf(stderr, "PSV/NSV: %.0f s\n", now() - t0);
    // greedy LZ77 factorization
    auto lcp = [&](uint64_t i, int32_t j) -> uint64_t {
        if (j < 0) return 0;
        uint64_t l = 0; while (i + l < N && T[j + l] == T[i + l]) ++l; return l;
    };
    uint64_t z = 0, lit = 0, i = 0, longest = 0;
    while (i < N) {
        uint64_t l = std::max(lcp(i, PSV[i]), lcp(i, NSV[i]));
        if (l == 0) { ++lit; l = 1; }
        longest = std::max(longest, l);
        ++z; i += l;
    }
    printf("%s: tokens %lu, alphabet %ld, LZ77 phrases z = %lu (%lu literals), tokens per phrase %.1f, longest phrase %lu\n",
           mode.c_str(), N, k, z, lit, (double)N / z, longest);
    fprintf(stderr, "total %.0f s\n", now() - t0);
    return 0;
}
