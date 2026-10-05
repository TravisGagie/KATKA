// Profile rz-index grid queries: where does the left/right time go?
// usage: rz-prof <index.rz> <patterns>
// Replays the rz-index query with timers around each component:
//   bs      : the two backward searches
//   map     : sparse-bitvector ranks mapping BWT intervals to grid ranges (every split)
//   trav    : wavelet-tree traversal finding the covering nodes (ranks on the way down)
//   rmq     : 1D RMQ at each covering node
//   descent : walking from the covering node down to the leaves to read the value (what -s shortens)
//   select  : phrase-rank -> text position
// and, for step s = 1, 2, 4, 8, 16, how many descent levels would remain if values were stored
// every s levels (levels are counted exactly; time is estimated with the measured cost per level).
#include "rz_index.hpp"
#include <fstream>
#include <x86intrin.h>
using namespace rz;
using std::string; using std::vector;

static vector<string> read_patterns(const string &f) {
    std::ifstream in(f, std::ios::binary);
    string header; std::getline(in, header);
    auto get = [&](const string &key) { size_t p = header.find(key + "="); return std::stoull(header.substr(p + key.size() + 1)); };
    u64 N = get("number"), m = get("length");
    vector<string> P(N, string(m, 0));
    for (u64 i = 0; i < N; ++i) in.read(&P[i][0], m);
    return P;
}

struct prof {
    u64 hy[8] = {0}, hx[8] = {0}, hmin[8] = {0}, empty_q = 0;   // widths: 1,2,3-4,5-8,9-16,17-64,65-1024,>1024
    u64 c_map = 0, c_trav = 0, c_rmq = 0, c_desc = 0, c_sel = 0, c_bs = 0;
    u64 nodes = 0, queries = 0, desc_levels = 0, trav_levels = 0;
    u64 steps_if[5] = {0, 0, 0, 0, 0};           // descent levels with s = 1,2,4,8,16
};
static const u64 SS[5] = {1, 2, 4, 8, 16};

struct pgrid {
    const grid_rmq &g; prof &P;
    void rec(u64 l, u64 prefix, u64 s, u64 e, u64 ya, u64 yb, u64 &best) {
        if (s >= e) return;
        u64 L = g.L, lo = prefix << (L - l), hi = (prefix + 1) << (L - l);
        if (hi <= ya || lo >= yb) return;
        if (ya <= lo && hi <= yb) {
            ++P.nodes;
            u64 t0 = __rdtsc();
            u64 p = (l == L || e - s == 1) ? s : g.rmq[l](s, e - 1);
            u64 t1 = __rdtsc();
            u64 ll = l, pp = p;
            while (g.val[ll].size() == 0) { u64 r1 = g.rk[ll](pp); pp = g.bv[ll][pp] ? g.Z[ll] + r1 : pp - r1; ++ll; }
            u64 v = g.val[ll][pp];
            u64 t2 = __rdtsc();
            P.c_rmq += t1 - t0; P.c_desc += t2 - t1; P.desc_levels += ll - l;
            for (int k = 0; k < 5; ++k) { u64 nx = std::min(L, (l + SS[k] - 1) / SS[k] * SS[k]); P.steps_if[k] += nx - l; }
            if (v < best) best = v;
            return;
        }
        ++P.trav_levels;
        u64 rs = g.rk[l](s), re = g.rk[l](e);
        rec(l + 1, prefix << 1, s - rs, e - re, ya, yb, best);
        rec(l + 1, (prefix << 1) | 1, g.Z[l] + rs, g.Z[l] + re, ya, yb, best);
    }
    static int bucket(u64 w) { return w <= 1 ? 0 : w <= 2 ? 1 : w <= 4 ? 2 : w <= 8 ? 3 : w <= 16 ? 4 : w <= 64 ? 5 : w <= 1024 ? 6 : 7; }
    u64 query(u64 x1, u64 x2, u64 ya, u64 yb) {
        ++P.queries;
        P.hy[bucket(yb - ya)]++; P.hx[bucket(x2 - x1)]++; P.hmin[bucket(std::min(yb - ya, x2 - x1))]++;
        u64 best = NONE;
        u64 t0 = __rdtsc(), c0 = P.c_rmq + P.c_desc;
        rec(0, 0, x1, x2, ya, yb, best);
        if (best == NONE) ++P.empty_q;
        P.c_trav += (__rdtsc() - t0) - (P.c_rmq + P.c_desc - c0);
        return best;
    }
};

int main(int argc, char **argv) {
    if (argc != 3) { std::cerr << "usage: rz-prof <index.rz> <patterns>\n"; return 1; }
    vector<string> Ps = read_patterns(argv[2]);
    rz::index Z; { std::ifstream in(argv[1], std::ios::binary); Z.load(in); }
    double w0 = now(); u64 k0 = __rdtsc();
    prof P;
    pgrid GL{Z.L.grid, P}, GR{Z.R.grid, P};
    vector<u64> fs, fe, xs, xe;
    for (auto &Pt : Ps) {
        u64 m = Pt.size(), t = __rdtsc();
        fs.assign(m + 1, 0); fe = fs; xs = fs; xe = fs;
        u64 sp = 0, ep = Z.bwt.n; fs[m] = 0; fe[m] = ep;
        bool ok = true;
        for (u64 i = m; i-- > 0;) { Z.bwt.extend(sp, ep, (unsigned char)Pt[i]); if (sp >= ep) { ok = false; break; } fs[i] = sp; fe[i] = ep; }
        if (!ok) continue;
        sp = 0; ep = Z.bwt.n; xs[0] = 0; xe[0] = ep;
        for (u64 i = 1; i <= m; ++i) { Z.bwt.extend(sp, ep, comp((unsigned char)Pt[i - 1])); xs[i] = sp; xe[i] = ep; }
        P.c_bs += __rdtsc() - t;
        for (u64 i = 1; i <= m; ++i) {
            t = __rdtsc();
            u64 qa = Z.L.Br.rank(xs[i]), qb = Z.L.Br.rank(xe[i]), qc = 0, qd = 0;
            if (qa < qb) { qc = Z.L.Bf.rank(fs[i]); qd = Z.L.Bf.rank(fe[i]); }
            P.c_map += __rdtsc() - t;
            if (qa >= qb || qc >= qd) continue;
            u64 k = GL.query(qa, qb, qc, qd);
            t = __rdtsc();
            if (k != NONE) { volatile u64 s = Z.L.Bp.select(k + 1); (void)s; }
            P.c_sel += __rdtsc() - t;
        }
        for (u64 i = 0; i < m; ++i) {
            t = __rdtsc();
            u64 qa, qb, qc = 0, qd = 0;
            if (i == 0) { qa = 0; qb = Z.R.z; } else { qa = Z.R.Br.rank(xs[i]); qb = Z.R.Br.rank(xe[i]); }
            if (qa < qb) { qc = Z.R.Bf.rank(fs[i]); qd = Z.R.Bf.rank(fe[i]); }
            P.c_map += __rdtsc() - t;
            if (qa >= qb || qc >= qd) continue;
            u64 v = GR.query(qa, qb, qc, qd);
            t = __rdtsc();
            if (v != NONE) { volatile u64 s = Z.R.Bp.select(Z.R.z - 1 - v + 1); (void)s; }
            P.c_sel += __rdtsc() - t;
        }
    }
    double hz = (__rdtsc() - k0) / (now() - w0);        // TSC ticks per second
    double us = 1e6 / hz / Ps.size();
    double lev = P.desc_levels ? P.c_desc / (double)P.desc_levels : 0;   // ticks per descent level
    u64 w = Z.L.grid.val[Z.L.grid.L].width();
    printf("m=%zu | us/pattern: bs=%.1f map=%.1f trav=%.1f rmq=%.1f descent=%.1f select=%.1f | "
           "per pattern: grid queries %.1f, covering nodes %.1f, descent levels %.1f (%.2f us/level)\n",
           Ps[0].size(), P.c_bs * us, P.c_map * us, P.c_trav * us, P.c_rmq * us, P.c_desc * us, P.c_sel * us,
           (double)P.queries / Ps.size(), (double)P.nodes / Ps.size(), (double)P.desc_levels / Ps.size(), lev / hz * 1e6);
    auto pr = [&](const char *nm, u64 *h) { printf("   %s width 1/2/<=4/<=8/<=16/<=64/<=1024/>1024:", nm); for (int k = 0; k < 8; ++k) printf(" %.1f%%", 100.0 * h[k] / std::max<u64>(1, P.queries)); printf("\n"); };
    pr("y", P.hy); pr("x", P.hx); pr("min", P.hmin);
    printf("   empty grid queries: %.1f%%\n", 100.0 * P.empty_q / std::max<u64>(1, P.queries));
    printf("   -s estimate (L=%lu/%lu levels, %lu-bit values):", Z.L.grid.L, Z.R.grid.L, w);
    for (int k = 0; k < 5; ++k) {
        u64 stored = (Z.L.grid.L + SS[k] - 1) / SS[k] + (Z.R.grid.L + SS[k] - 1) / SS[k] - 2;   // extra stored levels, both grids
        double extraGB = stored * (double)Z.L.z * w / 8 / 1e9;
        printf("  s=%lu: descent %.1f us, +%.2f GB", SS[k], P.steps_if[k] * lev * us, extraGB);
    }
    printf("\n");
}
