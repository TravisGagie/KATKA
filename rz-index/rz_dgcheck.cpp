// How much do minimizer digests loosen the leftmost/rightmost answers?
//
// usage: rz-dgcheck <digest.rz> <digest.tbl> <base.tbl> <base.S> k w <patterns1> <answers1> [<patterns2> <answers2> ...]
//   patterns: base (A/C/G/T) patterns in Pizza&Chili format; answers: the exact (leftmost, rightmost)
//   positions in the undigested S (rz-bench -o on the undigested index).  Each pattern is digested with
//   windows lying entirely inside it, queried on the digested rz-index, and the genomes (and species)
//   holding the digest answers are compared with those holding the exact answers.
#include "rz_index.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <fcntl.h>
#include <unistd.h>
using std::string; using std::vector; typedef uint64_t u64;

struct table { vector<u64> start, len; vector<string> sp, name; };
static table read_tbl(const string &f) {
    table t; std::ifstream in(f); string line;
    while (std::getline(in, line)) {
        std::istringstream ss(line); string idx, sp, file; u64 len, st;
        std::getline(ss, idx, '\t'); std::getline(ss, sp, '\t'); std::getline(ss, file, '\t'); ss >> len >> st;
        t.start.push_back(st); t.len.push_back(len); t.sp.push_back(sp); t.name.push_back(file);
    }
    return t;
}
static u64 genome(const table &t, u64 pos) { return std::upper_bound(t.start.begin(), t.start.end(), pos) - t.start.begin() - 1; }
static vector<string> read_patterns(const string &f) {
    std::ifstream in(f, std::ios::binary); string header; std::getline(in, header);
    auto get = [&](const string &key) { size_t p = header.find(key + "="); return std::stoull(header.substr(p + key.size() + 1)); };
    u64 N = get("number"), m = get("length");
    vector<string> P(N, string(m, 0));
    for (u64 i = 0; i < N; ++i) in.read(&P[i][0], m);
    return P;
}
static string rc_dna(const string &s) { string r(s.rbegin(), s.rend()); for (auto &c : r) c = (char)rz::rc_symbol((unsigned char)c); return r; }
int main(int argc, char **argv) {
    if (argc < 9 || (argc - 7) % 2) { fprintf(stderr, "usage: rz-dgcheck digest.rz digest.tbl base.tbl base.S k w pat ans [pat ans ...]\n"); return 1; }
    table D = read_tbl(argv[2]), B = read_tbl(argv[3]);
    int fd = open(argv[4], O_RDONLY);
    int k = atoi(argv[5]), w = atoi(argv[6]);
    // genomes are matched by name, species by name; species order may differ between the two runs
    std::map<string, u64> bidx; for (u64 i = 0; i < B.name.size(); ++i) bidx[B.name[i]] = i;
    vector<u64> d2b(D.name.size()); for (u64 i = 0; i < D.name.size(); ++i) d2b[i] = bidx.at(D.name[i]);
    vector<string> bsp, dsp;   // species orders
    for (auto &s : B.sp) if (bsp.empty() || bsp.back() != s) bsp.push_back(s);
    for (auto &s : D.sp) if (dsp.empty() || dsp.back() != s) dsp.push_back(s);
    std::map<string, int> brank, drank; for (int i = 0; i < (int)bsp.size(); ++i) brank[bsp[i]] = i; for (int i = 0; i < (int)dsp.size(); ++i) drank[dsp[i]] = i;
    // species whose relative order differs: a leftmost/rightmost species in base order is only
    // guaranteed to be the one in digest order if no species jumps across it
    auto order_safe = [&](const string &a, bool left) {
        for (auto &s : bsp) if (s != a) {
            bool bb = brank[s] < brank[a], db = drank[s] < drank[a];
            if (bb != db && (left ? db : !db)) return false;
        }
        return true;
    };
    auto contains = [&](u64 bg, const string &p, const string &rp) {
        string g(B.len[bg], 0); pread(fd, &g[0], B.len[bg], B.start[bg]);
        return g.find(p) != string::npos || g.find(rp) != string::npos;
    };
    rz::index Z; { std::ifstream in(argv[1], std::ios::binary); Z.load(in); }
    printf("bases\tn\tdglen\tbad\tfpL\tfpR\tfpAny\tspN\tspL\tspR\tspAny\tone\toneLost\n");
    for (int a = 7; a < argc; a += 2) {
        vector<string> P = read_patterns(argv[a]);
        vector<u64> ans; { std::ifstream in(argv[a + 1], std::ios::binary); u64 x; while (in.read((char *)&x, 8)) ans.push_back(x); }
        u64 n = 0, bad = 0, fL = 0, fR = 0, fA = 0, sN = 0, sL = 0, sR = 0, sA = 0, one = 0, oneL = 0; double dl = 0;
        for (u64 i = 0; i < P.size(); ++i) {
            if (P[i].find_first_not_of("ACGT") != string::npos || ans[2 * i] == rz::NONE) continue;
            string d = rz::digest(P[i], k, w), rp = rc_dna(P[i]);
            auto r = Z.query(d);
            ++n; dl += d.size();
            if (r.first == rz::NONE) { ++bad; continue; }
            u64 gl = d2b[genome(D, r.first)], gr = d2b[genome(D, r.second)];     // base genome ids
            u64 tl = genome(B, ans[2 * i]), tr = genome(B, ans[2 * i + 1]);
            bool okL = gl == tl || contains(gl, P[i], rp), okR = gr == tr || contains(gr, P[i], rp);
            fL += !okL; fR += !okR; fA += !okL || !okR;
            const string &TL = B.sp[tl], &TR = B.sp[tr];
            if (order_safe(TL, true) && order_safe(TR, false)) {
                ++sN;
                const string &GL = B.sp[gl], &GR = B.sp[gr];
                if (drank[GL] > drank[TL] || drank[GR] < drank[TR]) { ++bad; continue; }   // false negative: never
                bool a1 = GL != TL, a2 = GR != TR;
                sL += a1; sR += a2; sA += a1 || a2;
                if (TL == TR) { ++one; oneL += a1 || a2; }
            }
        }
        printf("%zu\t%lu\t%.1f\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n", P[0].size(), n, dl / n, bad, fL, fR, fA, sN, sL, sR, sA, one, oneL);
        fflush(stdout);
    }
}
