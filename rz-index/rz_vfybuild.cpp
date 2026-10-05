// Build the verification structures (vfy.hpp) for a byte-mapped digested reference built by rz-prep -B.
// usage: rz-vfybuild [-P old=new] <prefix.tbl> <prefix.map> <out.vfy>
//   prefix.tbl lists the documents in order (index, species, FASTA path, digest length, position in S);
//   -P rewrites a path prefix (for running where the FASTA files are mounted elsewhere).
// Checks that the per-run digests reproduce each document's digest exactly.
#include "digest.hpp"
#include "vfy.hpp"
#include <divsufsort.h>
#include <zlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <deque>
#include <chrono>
using std::string; typedef uint64_t u64;

static string read_genome(const string &path) {            // as rz-prep
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) { std::cerr << "cannot open " << path << "\n"; exit(1); }
    string g; std::vector<char> buf(1 << 20); bool header = false, bol = true, any = false; int got;
    while ((got = gzread(f, buf.data(), buf.size())) > 0)
        for (int i = 0; i < got; ++i) {
            char c = buf[i];
            if (bol && c == '>') { header = true; if (any) g.push_back('X'); }
            if (c == '\n') { header = false; bol = true; continue; }
            bol = false;
            if (header || c == '\r' || c == ' ' || c == '\t') continue;
            switch (c) { case 'A': case 'a': g.push_back('A'); break; case 'C': case 'c': g.push_back('C'); break;
                         case 'G': case 'g': g.push_back('G'); break; case 'T': case 't': g.push_back('T'); break; default: g.push_back('X'); }
            any = true;
        }
    gzclose(f); return g;
}

// digest of one A/C/G/T run as rz::digest_codes does it, with the first k-mer position of each symbol
static void digest_run_pos(const char *s, u64 len, int k, int w, std::vector<uint32_t> &codes, std::vector<uint32_t> &pos) {
    codes.clear(); pos.clear();
    if (len < (u64)k) return;
    u64 nk = len - k + 1, W = (u64)(w - k + 1); if (W > nk) W = nk;
    uint32_t mask = (uint32_t)((1ull << (2 * k)) - 1), fw = 0;
    std::vector<u64> h(nk); std::vector<uint32_t> f(nk);
    for (u64 t = 0; t < len; ++t) {
        fw = ((fw << 2) | (uint32_t)rz::dg_b2(s[t])) & mask;
        if (t + 1 >= (u64)k) { u64 p = t + 1 - k; uint32_t rc = rz::dg_rc_code(fw, k); f[p] = fw; h[p] = rz::dg_mix(fw < rc ? fw : rc); }
    }
    std::vector<char> sel(nk, 0); std::deque<u64> dq;
    for (u64 p = 0; p < nk; ++p) {
        while (!dq.empty() && h[dq.back()] > h[p]) dq.pop_back();
        dq.push_back(p);
        if (p + 1 >= W) { u64 s0 = p + 1 - W; while (dq.front() < s0) dq.pop_front(); u64 mn = h[dq.front()];
            for (u64 q : dq) { if (h[q] != mn) break; sel[q] = 1; } }
    }
    for (u64 p = 0; p < nk; ++p) if (sel[p]) { if (codes.empty() || codes.back() != f[p]) { codes.push_back(f[p]); pos.push_back((uint32_t)p); } }
}

int main(int argc, char **argv) {
    string from, to; int a = 1;
    if (argc > 2 && string(argv[1]) == "-P") { string p = argv[2]; size_t e = p.find('='); from = p.substr(0, e); to = p.substr(e + 1); a = 3; }
    if (argc - a != 3) { fprintf(stderr, "usage: rz-vfybuild [-P old=new] prefix.tbl prefix.map out.vfy\n"); return 1; }
    auto T0 = std::chrono::steady_clock::now();
    rz::dg_bytemap &M = rz::bytemap();
    if (!M.load(argv[a + 1])) { fprintf(stderr, "cannot load map\n"); return 1; }
    const int k = M.k, w = M.w;
    rz::vfy_index V;
    string T;                                          // all runs, concatenated
    std::vector<uint32_t> codes, pos, allpos;           // allpos: symbol positions, all runs
    std::vector<std::pair<u64, u64>> refrun;           // per document: longest run (tstart, len)
    std::ifstream tbl(argv[a]); string line; u64 ndoc = 0, bad = 0;
    while (std::getline(tbl, line)) {
        std::istringstream ss(line); string idx, sp, path, dl, ps;
        std::getline(ss, idx, '\t'); std::getline(ss, sp, '\t'); std::getline(ss, path, '\t'); std::getline(ss, dl, '\t'); std::getline(ss, ps, '\t');
        if (!from.empty() && path.compare(0, from.size(), from) == 0) path = to + path.substr(from.size());
        string g = read_genome(path);
        V.dfirst.push_back(V.tstart.size());
        string dg; u64 i = 0, n = g.size(); std::pair<u64, u64> best = {0, 0};
        while (i < n) {
            while (i < n && rz::dg_b2(g[i]) < 0) ++i;
            u64 j = i; while (j < n && rz::dg_b2(g[j]) >= 0) ++j;
            if (j - i >= (u64)k) {
                digest_run_pos(g.data() + i, j - i, k, w, codes, pos);
                if (!dg.empty()) dg.push_back('X');
                V.tstart.push_back(T.size()); V.rlen.push_back(j - i); V.dstart.push_back(dg.size()); V.dlen.push_back(codes.size());
                V.sfirst.push_back(allpos.size());
                for (u64 q = 0; q < codes.size(); ++q) { dg.push_back((char)M.code2byte[codes[q]]); allpos.push_back(pos[q]); }
                if (j - i > best.second) best = {T.size(), j - i};
                T.append(g, i, j - i);
            }
            i = j;
        }
        V.dlg.push_back(dg.size());
        if (dg != rz::digest_mapped(g) || dg.size() != std::stoull(dl)) ++bad;
        refrun.push_back(best); ++ndoc;
    }
    V.dfirst.push_back(V.tstart.size()); V.sfirst.push_back(allpos.size());
    if (bad) { fprintf(stderr, "ERROR: %lu documents' run digests differ from their digests\n", bad); return 1; }
    fprintf(stderr, "documents %lu, runs %zu, bases %zu, digest symbols %zu\n", ndoc, V.tstart.size(), T.size(), allpos.size());
    // symbol positions: 4-bit gaps with escapes, a sample every 32 symbols of each run
    V.gaps.assign((allpos.size() + 1) / 2, 0);
    for (u64 j = 0; j < V.tstart.size(); ++j) {
        V.sblock.push_back(V.samp_pos.size());
        u64 s0 = V.sfirst[j], s1 = V.sfirst[j + 1];
        for (u64 i = s0; i < s1; ++i) {
            uint32_t g = i == s0 ? allpos[i] : allpos[i] - allpos[i - 1];
            if (g >= 15) { V.gaps[i >> 1] |= 15 << ((i & 1) * 4); V.esc.push_back(g); } else V.gaps[i >> 1] |= g << ((i & 1) * 4);
            if ((i - s0) % 32 == 0) { V.samp_pos.push_back(allpos[i]); V.samp_esc.push_back(V.esc.size()); }
        }
    }
    V.sblock.push_back(V.samp_pos.size());
    // reference: the longest run of each document
    for (auto &r : refrun) V.R.append(T, r.first, r.second);
    std::vector<int32_t> SA(V.R.size());
    divsufsort((const unsigned char *)V.R.data(), SA.data(), (int32_t)V.R.size());
    const char *Rp = V.R.data(); u64 Rn = V.R.size();
    for (u64 j = 0; j < V.tstart.size(); ++j) {        // greedy RLZ parse of each run
        u64 st = V.tstart[j], en = st + V.rlen[j], i = st;
        while (i < en) {
            u64 lo = 0, hi = Rn, d = 0, blo = 0;
            while (i + d < en) {
                char c = T[i + d];
                u64 l = lo, h = hi;                    // first suffix in [lo,hi) with R[SA+d] >= c
                while (l < h) { u64 m = (l + h) / 2; u64 p = SA[m] + d; if (p < Rn && Rp[p] < c || p >= Rn) l = m + 1; else h = m; }
                u64 nl = l; h = hi;
                while (l < h) { u64 m = (l + h) / 2; u64 p = SA[m] + d; if (p < Rn && Rp[p] <= c || p >= Rn) l = m + 1; else h = m; }
                if (nl >= l) break;
                lo = nl; hi = l; blo = lo; ++d;
            }
            if (d == 0) { fprintf(stderr, "character not in reference\n"); return 1; }
            V.pstart.push_back(i); V.psrc.push_back(SA[blo]); i += d;
        }
    }
    V.pstart.push_back(T.size());
    // check random access on a sample of runs
    string out; u64 chk = 0;
    for (u64 j = 0; j < V.tstart.size(); j += 97) { V.run(j, out); if (out != T.substr(V.tstart[j], V.rlen[j])) ++chk;
        for (u64 q = 0; q < V.dlen[j]; ++q) if (V.sym_pos(j, q) != allpos[V.sfirst[j] + q]) { ++chk; break; } }
    if (chk) { fprintf(stderr, "ERROR: %lu runs decode wrongly\n", chk); return 1; }
    V.save(argv[a + 2]);
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    fprintf(stderr, "reference %zu bases, phrases %zu (%.1f per run, avg %.1f bases), escapes %zu\n", V.R.size(), V.pstart.size() - 1,
            (double)(V.pstart.size() - 1) / V.tstart.size(), (double)T.size() / (V.pstart.size() - 1), V.esc.size());
    fprintf(stderr, "bytes %lu: RLZ %lu (reference %zu + phrases %lu), symbol positions %lu, runs %lu | %.1f s\n", V.bytes(),
            V.R.size() + 8 * V.pstart.size(), V.R.size(), 8 * V.pstart.size(),
            V.gaps.size() + 2 * V.esc.size() + 2 * V.samp_pos.size() + 4 * V.samp_esc.size(), 4 * 4 * V.tstart.size() + 16 * V.tstart.size(), sec);
    return 0;
}
