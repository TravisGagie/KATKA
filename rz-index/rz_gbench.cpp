// Benchmark: counting genera in BWT intervals with the run-length tag array (tag.hpp) versus with the run starts of the
// tag array plus a RePair grammar over the genus IDs of its runs (bigrepair -i on the genus sequence, see
// tools/tagz/tokwrite.cpp).  The grammar is stored plainly here, with every rule's expansion length (in runs) in an
// array; a compact encoding would change the space, not the access pattern.
// usage: rz-gbench bac.s1.tag tags_genus.int32 intervals.bin [reads]
//   intervals.bin: pairs of u64 [sp, ep), written by rz-classify with RZ_IVLOG=file (RZ_COUNTS=1, -T, one thread)
#include "tag.hpp"
#include <chrono>
#include <cstdio>
using namespace rz;
struct Grammar {
    int32_t alpha = 0; std::vector<int32_t> Lc, Rc, C; std::vector<uint32_t> len; std::vector<uint64_t> pre; int height = 0;
    std::vector<int32_t> st;
    bool load(const std::string &p) {
        FILE *f = fopen((p + ".R").c_str(), "rb"); if (!f) return false;
        fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
        if (fread(&alpha, 4, 1, f) != 1) return false; long R = (sz - 4) / 8;
        std::vector<int32_t> rl(2 * R); if (fread(rl.data(), 8, R, f) != (size_t)R) return false; fclose(f);
        Lc.resize(R); Rc.resize(R); for (long i = 0; i < R; ++i) { Lc[i] = rl[2 * i]; Rc[i] = rl[2 * i + 1]; }
        len.assign(alpha + R, 1); std::vector<uint8_t> h(alpha + R, 0);
        for (long i = 0; i < R; ++i) { len[alpha + i] = len[Lc[i]] + len[Rc[i]]; h[alpha + i] = 1 + std::max(h[Lc[i]], h[Rc[i]]); height = std::max<int>(height, h[alpha + i]); }
        f = fopen((p + ".C").c_str(), "rb"); if (!f) return false; fseek(f, 0, SEEK_END); long n = ftell(f) / 4; fseek(f, 0, SEEK_SET);
        C.resize(n); if (fread(C.data(), 4, n, f) != (size_t)n) return false; fclose(f);
        pre.resize(n + 1); pre[0] = 0; for (long i = 0; i < n; ++i) pre[i + 1] = pre[i] + len[C[i]];
        st.reserve(256); return true;
    }
    // calls f(genus) for the runs a, a+1, ..., a+cnt-1
    template <class F> inline void decode(uint64_t a, uint64_t cnt, F f) {
        uint64_t k = std::upper_bound(pre.begin(), pre.end(), a) - pre.begin() - 1, off = a - pre[k];
        st.clear(); int32_t x = C[k];
        while (x >= alpha) { int32_t l = Lc[x - alpha], r = Rc[x - alpha]; if (off < len[l]) { st.push_back(r); x = l; } else { off -= len[l]; x = r; } }
        f(x);
        while (--cnt) {
            if (st.empty()) x = C[++k]; else { x = st.back(); st.pop_back(); }
            while (x >= alpha) { st.push_back(Rc[x - alpha]); x = Lc[x - alpha]; }
            f(x);
        }
    }
};
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: rz-gbench bac.s1.tag tags_genus.int32 intervals.bin [reads]\n"); return 1; }
    tag_index X; if (!X.load(argv[1], false)) { fprintf(stderr, "cannot load %s\n", argv[1]); return 1; }
    Grammar G; double t0 = now(); if (!G.load(argv[2])) { fprintf(stderr, "cannot load grammar %s\n", argv[2]); return 1; }
    if (G.pre.back() != X.rho) { fprintf(stderr, "grammar expands to %lu runs, tag array has %lu\n", G.pre.back(), X.rho); return 1; }
    fprintf(stderr, "grammar: %zu rules, |C| %zu, height %d, loaded in %.0f s\n", G.Lc.size(), G.C.size(), G.height, now() - t0);
    std::vector<u64> iv; { FILE *f = fopen(argv[3], "rb"); u64 v[2]; while (fread(v, 8, 2, f) == 2) { iv.push_back(v[0]); iv.push_back(v[1]); } fclose(f); }
    size_t m = iv.size() / 2; double reads = argc > 4 ? atof(argv[4]) : 0;
    // runs per interval
    uint64_t runs = 0, maxr = 0; for (size_t i = 0; i < m; ++i) { u64 r = X.run_of(iv[2 * i + 1] - 1) - X.run_of(iv[2 * i]) + 1; runs += r; maxr = std::max(maxr, r); }
    printf("%zu intervals, %.2f runs each on average (max %lu)\n", m, (double)runs / m, maxr);
    std::vector<std::pair<u64, u64>> out; std::vector<uint32_t> cnt(1 << 16, 0); std::vector<u64> marked; uint64_t chkA = 0, chkB = 0;
    for (int rep = 0; rep < 2; ++rep) {   // second repetition is the one reported (warm caches)
        double ta = now();
        for (size_t i = 0; i < m; ++i) { X.count(iv[2 * i], iv[2 * i + 1], out); for (auto &p : out) chkA += p.first * 31 + p.second; }
        ta = now() - ta;
        double tb = now();
        for (size_t i = 0; i < m; ++i) {
            u64 sp = iv[2 * i], ep = iv[2 * i + 1], a = X.run_of(sp), b = X.run_of(ep - 1), j = a;
            out.clear();
            G.decode(a, b - a + 1, [&](int32_t g) {
                u64 st = std::max(sp, (u64)X.RBs(j + 1)), en = std::min(ep, j + 1 < X.rho ? (u64)X.RBs(j + 2) : X.n);
                if (cnt[g] == 0) marked.push_back(g); cnt[g] += (uint32_t)(en - st); ++j;
            });
            for (u64 t : marked) { out.push_back({t, cnt[t]}); cnt[t] = 0; } marked.clear();
            std::sort(out.begin(), out.end());
            for (auto &p : out) chkB += p.first * 31 + p.second;
        }
        tb = now() - tb;
        if (rep == 1) {
            printf("tag array: %.1f ns per interval%s\n", 1e9 * ta / m, reads ? (", " + std::to_string(1e6 * ta / reads).substr(0, 5) + " us per read").c_str() : "");
            printf("grammar:   %.1f ns per interval%s\n", 1e9 * tb / m, reads ? (", " + std::to_string(1e6 * tb / reads).substr(0, 5) + " us per read").c_str() : "");
            printf("answers %s\n", chkA == chkB ? "identical" : "DIFFER");
        }
        chkA = chkB = 0;
    }
    return 0;
}
