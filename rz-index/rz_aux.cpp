// Build the auxiliary structures (union bitvectors, scan arrays) for an rz-index.
// usage: rz-aux <index.rz> <out.aux>
#include "rz_index.hpp"
#include <fstream>
using namespace rz; using std::vector;

static void build_union(const sparse_bv &A, const sparse_bv &B, u64 n, union_bv &U) {
    u64 za = A.rank(n), zb = B.rank(n);
    vector<u64> u; u.reserve(za + zb);
    vector<char> fa, fb; fa.reserve(za + zb); fb.reserve(za + zb);
    u64 i = 0, j = 0;
    u64 a = za ? A.select(1) : NONE, b = zb ? B.select(1) : NONE;
    while (i < za || j < zb) {
        u64 x = std::min(a, b);
        u.push_back(x); fa.push_back(a == x); fb.push_back(b == x);
        if (a == x) { ++i; a = i < za ? A.select(i + 1) : NONE; }
        if (b == x) { ++j; b = j < zb ? B.select(j + 1) : NONE; }
    }
    U.U.build(u, n);
    // one bit per union element; rank over [0, j) = number of left/right marks among the first j
    U.fl = bit_vector(u.size() + 1, 0); U.fr = bit_vector(u.size() + 1, 0);
    for (u64 k = 0; k < u.size(); ++k) { U.fl[k] = fa[k]; U.fr[k] = fb[k]; }
    util::init_support(U.rl, &U.fl); util::init_support(U.rr, &U.fr);
}
static void build_scan(const grid_rmq &g, scan_arrays &S) {
    // follow every point level by level (stable partition, as the wavelet matrix was built)
    u64 z = g.z, w = bits::hi(std::max<u64>(1, z - 1)) + 1;
    vector<uint32_t> cx(z), cy(z, 0), nx(z), ny(z);
    for (u64 p = 0; p < z; ++p) cx[p] = (uint32_t)p;
    for (u64 l = 0; l < g.L; ++l) {
        u64 a = 0, b = g.Z[l];
        const bit_vector &bv = g.bv[l];
        for (u64 p = 0; p < z; ++p) {
            bool bit = bv[p];
            uint32_t y = (cy[p] << 1) | bit;
            if (bit) { nx[b] = cx[p]; ny[b] = y; ++b; } else { nx[a] = cx[p]; ny[a] = y; ++a; }
        }
        cx.swap(nx); cy.swap(ny);
    }
    S.yx = int_vector<>(z, 0, w); S.xy = int_vector<>(z, 0, w); S.kx = int_vector<>(z, 0, w);
    for (u64 p = 0; p < z; ++p) { u64 x = cx[p], y = cy[p]; S.yx[x] = y; S.xy[y] = x; S.kx[x] = g.val[g.L][p]; }
}
int main(int argc, char **argv) {
    if (argc != 3) { std::cerr << "usage: rz-aux <index.rz> <out.aux>\n"; return 1; }
    double t0 = now();
    rz::index Z; { std::ifstream in(argv[1], std::ios::binary); Z.load(in); }
    u64 n = Z.bwt.n;
    aux_index A;
    build_union(Z.L.Bf, Z.R.Bf, n, A.F);
    build_union(Z.L.Br, Z.R.Br, n, A.Rr);
    build_scan(Z.L.grid, A.SL);
    build_scan(Z.R.grid, A.SR);
    std::ofstream o(argv[2], std::ios::binary);
    u64 b = A.serialize(o);
    struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
        std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } };
    counter c1, c2, c3, c4; std::ostream o1(&c1), o2(&c2), o3(&c3), o4(&c4);
    A.F.serialize(o1); A.Rr.serialize(o2);
    Z.L.Bf.serialize(o3); Z.R.Bf.serialize(o3); Z.L.Br.serialize(o4); Z.R.Br.serialize(o4);
    printf("aux bytes=%lu: unions %lu (replacing %lu of Bf/Br), scan arrays %lu; time=%.1fs\n",
           b, c1.c + c2.c, c3.c + c4.c, b - c1.c - c2.c, now() - t0);
}
