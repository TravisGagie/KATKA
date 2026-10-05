// Build a subsampled r-index (sr.hpp) from the full baseline r-index and the move structure.
// usage: sr-build <index.rix> <index.mv> <s> <out.sri>
#include "sr.hpp"
#include <fstream>
#include <numeric>
using namespace rz; using std::vector;
int main(int argc, char **argv) {
    // memory: the full r-index plus 4 bytes per phi entry (sort order) and the output
    if (argc != 5) { std::cerr << "usage: sr-build <index.rix> <index.mv (unused)> <s[,s2,...]> <out.sri | out-prefix>\n"
                                 "  with several values of s, writes <out-prefix>.s<s>.sri for each\n"; return 1; }
    double t0 = now();
    rz::rindex F; { std::ifstream in(argv[1], std::ios::binary); F.load(in); }
    const rlbwt &B = F.bwt;
    u64 n = F.n, R = B.R;
    vector<u64> svals;
    { std::string a = argv[3]; size_t p = 0; while (p <= a.size()) { size_t q = a.find(',', p); if (q == std::string::npos) q = a.size(); svals.push_back(std::stoull(a.substr(p, q - p))); p = q + 1; } }
    bool multi = svals.size() > 1;
    auto chr = [&](u64 k) -> unsigned char { return B.heads[k]; };
    auto rstart = [&](u64 k) -> u64 { return k == R ? n : B.starts.select(k + 1); };
    u64 M = F.keys.rank(n);
    auto val = [&](u64 j) -> u64 { u64 v = F.vals[j]; return v < R ? (u64)F.esa[v] : (u64)F.extra[v - R]; };
    auto req = [&](u64 j) -> bool { u64 v = F.vals[j]; return v >= R || chr(v) == 1; };
    // entries in text order of their (distinct) values: mark the values in a bitvector over [0, n)
    // and rank them, instead of sorting with random accesses
    vector<uint32_t> ord(M);
    {
        bit_vector vb(n, 0);
        for (u64 j = 0; j < M; ++j) vb[val(j)] = 1;
        rank_support_v5<1> vr(&vb);
        for (u64 j = 0; j < M; ++j) ord[vr(val(j))] = (uint32_t)j;
    }
    printf("ordered %lu entries in %.1fs\n", M, now() - t0); fflush(stdout);
    // kept entries for every s (M bits each); then the order is no longer needed
    vector<bit_vector> keptv;
    for (u64 s : svals) {
    keptv.emplace_back(M, 0);
    bit_vector &kept = keptv.back();
    if (M) {
        kept[ord[0]] = 1;
        u64 lastv = val(ord[0]);
        for (u64 i = 2; i < M; ++i) {
            u64 prev = ord[i - 1], cur = ord[i];
            if (val(cur) - lastv > s || req(prev)) { kept[prev] = 1; lastv = val(prev); }
        }
        kept[ord[M - 1]] = 1;
    }
    }
    vector<uint32_t>().swap(ord);
    for (u64 si = 0; si < svals.size(); ++si) {
    u64 s = svals[si];
    double t1 = now();
    bit_vector kept;
    kept.swap(keptv[si]);
    srindex X;
    X.n = n; X.s = s; X.R = R;
    u64 Mk = 0; for (u64 j = 0; j < M; ++j) Mk += kept[j];
    {
        sd_vector_builder kb(n, Mk);
        X.vals = int_vector<>(Mk, 0, bits::hi(n) + 1);
        X.trusted = bit_vector(Mk, 0);
        vector<u64> ar;
        u64 t = 0, knext = M ? F.keys.select(1) : 0;
        for (u64 j = 0; j < M; ++j) {
            u64 kj = knext;
            if (j + 1 < M) knext = F.keys.select(j + 2);
            if (!kept[j]) continue;
            kb.set(kj); X.vals[t] = val(j);
            bool tr = j + 1 == M || kept[j + 1];
            X.trusted[t] = tr;
            if (!tr) ar.push_back(knext - kj);
            ++t;
        }
        X.keys.build_from(kb);
        u64 mxa = 1; for (u64 a : ar) mxa = std::max(mxa, a);
        X.area = int_vector<>(ar.size(), 0, bits::hi(mxa) + 1);
        for (u64 i = 0; i < ar.size(); ++i) X.area[i] = ar[i];
        util::init_support(X.trusted_r0, &X.trusted);
    }
    X.ekept = bit_vector(R, 0);
    for (u64 j = 0; j < M; ++j) if (kept[j]) { u64 v = F.vals[j]; if (v < R) X.ekept[v] = 1; }
    X.ekept[R - 1] = 1;
    util::init_support(X.ekept_r, &X.ekept);
    u64 E = X.ekept_r(R);
    X.esav = int_vector<>(E, 0, bits::hi(n) + 1);
    for (u64 k = 0, t = 0; k < R; ++k) if (X.ekept[k]) X.esav[t++] = F.esa[k];
    // terminator rows: SA of the last row of each terminator run is a run-end sample; the rows
    // above it follow by (exact) phi
    for (u64 k = 0; k < R; ++k) if (chr(k) == 1) {
        u64 st = rstart(k), l = rstart(k + 1) - st;
        vector<u64> sa(l);
        sa[l - 1] = F.esa[k];
        for (u64 x = l - 1; x-- > 0;) sa[x] = F.phi(sa[x + 1]);
        for (u64 x = 0; x < l; ++x) { X.trow.push_back(st + x); X.tpos.push_back(sa[x]); }
    }
    {
        u64 K = F.dstart.rank(n);
        vector<u64> ds(K);
        for (u64 j = 0; j < K; ++j) ds[j] = F.dstart.select(j + 1);
        X.dstart.build(ds, n);
    }
    std::string of = multi ? std::string(argv[4]) + ".s" + std::to_string(s) + ".sri" : std::string(argv[4]);
    std::ofstream o(of, std::ios::binary);
    u64 b = X.serialize(o);
    printf("s=%lu entries kept %lu / %lu, run-end samples kept %lu / %lu, untrusted %lu, bytes=%lu time=%.1fs\n",
           s, Mk, M, E, R, (u64)X.area.size(), b, now() - t1);
    fflush(stdout);
    }
}
