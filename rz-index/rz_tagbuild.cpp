// Build the tag index (tag.hpp) of a document collection S built by rz-prep (one document per genus).
// usage: rz-tagbuild <prefix> <s> <out.tag> [-n]
//   writes out.tag and, unless -n, out.tag.rmq, the RMQ that document listing needs (counting does not).
//   reads prefix.tbl, prefix.rix (rz-build -A: the r-index's SA samples at BWT run boundaries) and, for s > 1,
//   prefix.rz (its RLBWT, for LF); s = 1: every run's tag is stored.
//   Stage 1 (the runs and their tags) is cached in prefix.tagruns.  It builds no suffix array: it streams the
//   suffix array from the bottom of the BWT up with phi (SA[x-1] = phi(SA[x])), in memory proportional to the
//   number of BWT runs plus the output, and maps each entry to its document by a predecessor search.
#include "tag.hpp"
#include <fstream>
#include <sstream>
#include <deque>
#include <chrono>
using namespace rz;
static double secs(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }
int main(int argc, char **argv) {
    bool norm = argc == 5 && std::string(argv[4]) == "-n";   // -n: no RMQ (counting only; listing needs it)
    if (argc != 4 && !norm) { fprintf(stderr, "usage: rz-tagbuild prefix s out.tag [-n]\n  writes out.tag and, unless -n, out.tag.rmq (for listing)\n"); return 1; }
    std::string p = argv[1]; u64 s = std::stoull(argv[2]); auto T0 = std::chrono::steady_clock::now();
    std::vector<uint32_t> rs; std::vector<uint16_t> tg; u64 n = 0;
    {   // stage 1
        std::ifstream c(p + ".tagruns", std::ios::binary);
        if (c) {
            u64 r; c.read((char *)&n, 8); c.read((char *)&r, 8); rs.resize(r); tg.resize(r);
            c.read((char *)rs.data(), 4 * r); c.read((char *)tg.data(), 2 * r);
        } else {
            std::vector<u64> start; std::ifstream tb(p + ".tbl"); std::string line;
            while (std::getline(tb, line)) { std::istringstream ss(line); std::string x; for (int i = 0; i < 5; ++i) std::getline(ss, x, '\t'); start.push_back(std::stoull(x)); }
            std::ifstream ri_in(p + ".rix", std::ios::binary);
            if (!ri_in) { fprintf(stderr, "rz-tagbuild: %s.rix not found (build it with rz-build -A)\n", p.c_str()); return 1; }
            rz::rindex RI; RI.load(ri_in);
            n = RI.n; u64 N = n - 1;                    // row 0 is the terminator
            if (RI.dstart.rank(n) != 1) { fprintf(stderr, "rz-tagbuild: expected one string (one terminator)\n"); return 1; }
            fprintf(stderr, "r-index samples loaded (%.0f s); streaming the suffix array\n", secs(T0));
            u64 d = RI.esa[RI.bwt.R - 1];               // SA[n-1]
            int prev = -1; u64 g0 = 0, g1 = start.size() > 1 ? start[1] : ~0ULL;   // cached document [g0, g1)
            int cd = 0;
            for (u64 row = n; row-- > 0;) {
                u64 q = d < N ? d : N - 1;              // row 0: give it the last document
                if (q < g0 || q >= g1) {
                    cd = std::upper_bound(start.begin(), start.end(), q) - start.begin() - 1;
                    g0 = start[cd]; g1 = cd + 1 < (int)start.size() ? start[cd + 1] : ~0ULL;
                }
                if (cd != prev) { rs.push_back(row); tg.push_back(cd); prev = cd; } else rs.back() = row;   // runs found bottom-up
                if (row) d = RI.phi(d);
            }
            std::reverse(rs.begin(), rs.end()); std::reverse(tg.begin(), tg.end());
            if (rs.empty() || rs[0] != 0) { fprintf(stderr, "rz-tagbuild: internal error: no run at row 0\n"); return 1; }
            std::ofstream o(p + ".tagruns", std::ios::binary); u64 r = rs.size();
            o.write((char *)&n, 8); o.write((char *)&r, 8); o.write((char *)rs.data(), 4 * r); o.write((char *)tg.data(), 2 * r);
        }
    }
    u64 rho = rs.size();
    fprintf(stderr, "n=%lu runs=%lu (stage 1 %.0f s)\n", n, rho, secs(T0));
    tag_index X; X.n = n; X.rho = rho; X.s = s;
    { sdsl::bit_vector b(n, 0); for (u64 r : rs) b[r] = 1; X.RB = sdsl::sd_vector<>(b); }
    sdsl::util::init_support(X.RBr, &X.RB); sdsl::util::init_support(X.RBs, &X.RB);
    if (!norm) {   // RMQ over C (previous run with the same tag, +1; 0 if none)
        sdsl::int_vector<> C(rho, 0, sdsl::bits::hi(rho + 1) + 1); std::vector<uint32_t> last(1 << 16, 0);
        for (u64 i = 0; i < rho; ++i) { C[i] = last[tg[i]]; last[tg[i]] = i + 1; }
        X.rmq = sdsl::rmq_succinct_sct<>(&C); X.has_rmq = true;
    }
    u64 nroot = 0, ncyc = 0, nsamp = rho;
    if (s == 1) {
        X.L = sdsl::int_vector<>(rho, 0, 14); for (u64 i = 0; i < rho; ++i) X.L[i] = tg[i];
    } else {
        rz::index Z; { std::ifstream in(p + ".rz", std::ios::binary); Z.load(in); }
        const rlbwt &B = Z.bwt;
        std::vector<uint32_t> e(rho); std::vector<char> root(rho, 0);
        for (u64 i = 0; i < rho; ++i) {
            u64 b = rs[i]; unsigned char c = B.at(b);
            u64 lf = B.C[c] + B.rank(b, c);
            e[i] = std::upper_bound(rs.begin(), rs.end(), (uint32_t)lf) - rs.begin() - 1;
            if (i == 0 || tg[e[i]] != tg[i]) { root[i] = 1; ++nroot; }
        }
        fprintf(stderr, "LF of run heads done (%.0f s), roots %lu\n", secs(T0), nroot);
        // greedy sampling on the out-degree-1 graph i -> e[i] (roots have no out-edge)
        std::vector<uint32_t> indeg(rho, 0); std::vector<uint8_t> need(rho, 0); std::vector<char> samp(rho, 0), done(rho, 0);
        for (u64 i = 0; i < rho; ++i) if (!root[i]) ++indeg[e[i]];
        std::deque<u64> Q; for (u64 i = 0; i < rho; ++i) if (!indeg[i]) Q.push_back(i);
        auto process = [&](u64 v, bool force) {
            done[v] = 1;
            if (root[v] || force || need[v] + 1 >= s) { samp[v] = 1; return; }
            u64 w = e[v]; if (need[v] + 1 > need[w]) need[w] = need[v] + 1;
        };
        while (!Q.empty()) {
            u64 v = Q.front(); Q.pop_front();
            process(v, false);
            if (!root[v] && --indeg[e[v]] == 0) Q.push_back(e[v]);
        }
        for (u64 c = 0; c < rho; ++c) {               // what is left lies on cycles
            if (done[c]) continue;
            ++ncyc; process(c, true);
            for (u64 v = e[c]; v != c; v = e[v]) process(v, false);
        }
        nsamp = 0; X.SB = sdsl::bit_vector(rho, 0);
        for (u64 i = 0; i < rho; ++i) if (samp[i]) { X.SB[i] = 1; ++nsamp; }
        sdsl::util::init_support(X.SBr, &X.SB);
        X.L = sdsl::int_vector<>(nsamp, 0, 14);
        for (u64 i = 0, j = 0; i < rho; ++i) if (samp[i]) X.L[j++] = tg[i];
        // check every run's tag
        X.bwt = &B; u64 bad = 0, maxsteps = 0;
        u64 stepc = s > 4 ? 97 : 1, nchk = 0;
        for (u64 i = 0; i < rho; i += stepc) { ++nchk; u64 b0 = X.lf_steps; if (X.tag(i) != tg[i]) ++bad; maxsteps = std::max(maxsteps, X.lf_steps - b0); }
        if (bad) { fprintf(stderr, "ERROR: %lu runs get the wrong tag\n", bad); return 1; }
        fprintf(stderr, "%lu tags checked (every %lu-th run); max LF steps %lu, mean %.2f\n", nchk, stepc, maxsteps, (double)X.lf_steps / nchk);
    }
    X.save(argv[3]);
    fprintf(stderr, "s=%lu: runs %lu, roots %lu, cycles %lu, sampled %lu (%.1f%%) | bytes %lu (%.2f bits/run): starts %lu, sampled bits %lu, tags %lu, RMQ %lu | %.0f s\n",
            s, rho, nroot, ncyc, nsamp, 100.0 * nsamp / rho, X.bytes(), 8.0 * X.bytes() / rho,
            sdsl::size_in_bytes(X.RB) + sdsl::size_in_bytes(X.RBr) + sdsl::size_in_bytes(X.RBs), s > 1 ? sdsl::size_in_bytes(X.SB) + sdsl::size_in_bytes(X.SBr) : 0,
            sdsl::size_in_bytes(X.L), X.has_rmq ? sdsl::size_in_bytes(X.rmq) : 0, secs(T0));
    return 0;
}
