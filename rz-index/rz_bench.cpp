// Benchmark: leftmost and rightmost occurrences with the rz-index vs the r-index.
//
// usage: rz-bench [-t S] [-o answers] <index.rz | index.rix> <patterns>
//          rz-bench -c <answers1> <answers2>
//   One index per process (its type is taken from the file extension), so only one index is in
//   memory at a time.  The index size reported is its file size.
//   -t S : also check every answer by brute force (find / rfind on S held in memory)
//   -o f : write the (leftmost, rightmost) answers to f, for comparison with -c
//   -c   : count the patterns on which two answer files disagree
//   -m f : (rz-index) do the backward searches on the move structure f (rz-mvbuild); the size
//          reported then counts the move structure instead of the RLBWT
//   -x f : (rz-index) use the auxiliary structures f (rz-aux); -T t sets the scan threshold (16)
//   An index file ending in .mv times only the rz-index's two backward searches (for P and rc(P),
//   recording BWT rows at every step) on that move structure, for indexes too large to load
//   together with it.
//   -D k,w : the patterns are DNA; digest each one (minimizers, windows lying inside the pattern) before
//          querying a digested index.  m is then the length in bases; the average digest length is printed.
//   -B f : the patterns are DNA and the index is over a byte-mapped digest (rz-prep -B, map file f): digest
//          each pattern with the map (whole pattern); m is then the length in bases
//   -n N : use only the first N patterns
//   -C f : (rz-index) do the backward searches on the RLCSA f (rz-csabuild); the size reported then
//          counts the RLCSA instead of the RLBWT
//   -b f : (sr-index) do the backward searches and LF steps on the RLBWT stored in the r-index f
//          (.rix), i.e. an RLFM-index backend, instead of a move structure
//
// rz-index time = backward searches for P and rc(P) + left grid + right grid.
// r-index time  = backward search for P with toehold + phi loop over all occurrences
//                 (tracking min and max).
// Totals and breakdowns are reported as averages per pattern (µs).

#include "sr.hpp"
#include <fstream>
#include <sstream>
#include <unistd.h>

using std::string;
using std::vector;
typedef uint64_t u64;

static vector<string> read_patterns(const string &f) {
    std::ifstream in(f, std::ios::binary);
    if (!in) { std::cerr << "cannot open " << f << "\n"; exit(1); }
    string header; std::getline(in, header);
    auto get = [&](const string &key) {
        size_t p = header.find(key + "=");
        if (p == string::npos) { std::cerr << "bad header\n"; exit(1); }
        return std::stoull(header.substr(p + key.size() + 1));
    };
    u64 N = get("number"), m = get("length");
    vector<string> P(N, string(m, 0));
    for (u64 i = 0; i < N; ++i) in.read(&P[i][0], m);
    return P;
}

int main(int argc, char **argv) {
    string sfile, afile, mvfile, rlfile, auxfile, dgopt, csafile, mapfile;
    u64 maxpat = UINT64_MAX;
    u64 scanT = 16;
    bool cmp = false;
    int opt;
    const char *usage = "usage: rz-bench [-t S] [-o answers] <index.rz|index.rix> <patterns>\n"
                        "       rz-bench -c <answers1> <answers2>\n";
    while ((opt = getopt(argc, argv, "t:o:cm:b:x:T:D:n:C:B:")) != -1) {
        if (opt == 't') sfile = optarg;
        else if (opt == 'o') afile = optarg;
        else if (opt == 'c') cmp = true;
        else if (opt == 'm') mvfile = optarg;
        else if (opt == 'b') rlfile = optarg;
        else if (opt == 'x') auxfile = optarg;
        else if (opt == 'T') scanT = std::stoull(optarg);
        else if (opt == 'D') dgopt = optarg;
        else if (opt == 'C') csafile = optarg;
        else if (opt == 'B') mapfile = optarg;
        else if (opt == 'n') maxpat = std::stoull(optarg);
        else { std::cerr << usage; return 1; }
    }
    if (argc - optind != 2) { std::cerr << usage; return 1; }
    if (cmp) {
        auto rd = [](const char *f) {
            std::ifstream in(f, std::ios::binary); vector<u64> v;
            u64 x; while (in.read((char *)&x, 8)) v.push_back(x);
            return v;
        };
        vector<u64> a = rd(argv[optind]), b = rd(argv[optind + 1]);
        if (a.size() != b.size()) { printf("answer files differ in length (%zu vs %zu)\n", a.size() / 2, b.size() / 2); return 2; }
        u64 dis = 0;
        for (u64 i = 0; i + 1 < a.size(); i += 2) if (a[i] != b[i] || a[i + 1] != b[i + 1]) ++dis;
        printf("patterns=%zu disagreements=%lu\n", a.size() / 2, dis);
        return dis ? 2 : 0;
    }
    string ifile = argv[optind];
    bool isrz = ifile.size() > 3 && ifile.substr(ifile.size() - 3) == ".rz";
    bool issr = ifile.size() > 4 && ifile.substr(ifile.size() - 4) == ".sri";
    bool ismv = ifile.size() > 3 && ifile.substr(ifile.size() - 3) == ".mv";
    if (issr && mvfile.empty() && rlfile.empty()) { std::cerr << "an .sri index needs a backend: -m index.mv or -b index.rix\n"; return 1; }
    vector<string> P = read_patterns(argv[optind + 1]);
    if (P.size() > maxpat) P.resize(maxpat);
    u64 NP = P.size(), m = NP ? P[0].size() : 0;
    double dglen = 0;
    if (!mapfile.empty()) {
        if (!rz::bytemap().load(mapfile)) { std::cerr << "cannot load map " << mapfile << "\n"; return 1; }
        dgopt = "map" + std::to_string(rz::bytemap().k) + "," + std::to_string(rz::bytemap().w);
        for (auto &p : P) {
            if (p.find_first_not_of("ACGT") != string::npos) { std::cerr << "-B: pattern with a non-ACGT base\n"; return 1; }
            p = rz::digest_mapped(p); dglen += p.size();
        }
        dglen /= NP;
    } else if (!dgopt.empty()) {
        int k = std::stoi(dgopt), w = std::stoi(dgopt.substr(dgopt.find(',') + 1));
        for (auto &p : P) {
            if (p.find_first_not_of("ACGT") != string::npos) { std::cerr << "-D: pattern with a non-ACGT base\n"; return 1; }
            if (p.size() < (size_t)w) { std::cerr << "-D: pattern shorter than the window\n"; return 1; }
            p = rz::digest(p, k, w); dglen += p.size();
        }
        dglen /= NP;
    }
    u64 bytes; { std::ifstream f(ifile, std::ios::binary | std::ios::ate); bytes = f.tellg(); }

    vector<std::pair<u64, u64>> ans(NP);
    rz::timing T;
    u64 lq = 0, rq = 0, occ = 0, maxocc = 0, srs = 0, srlf = 0, srinv = 0;
    double t0, tt;
    bool X_ef = false;
    if (ismv) {
        rz::move_bwt M;
        { std::ifstream in(ifile, std::ios::binary); M.load(in); }
        std::vector<u64> fs, fe, xs, xe;
        u64 sink = 0;
        t0 = rz::now();
        for (u64 i = 0; i < NP; ++i) {
            const std::string &Q = P[i];
            const u64 mq = Q.size();
            fs.resize(mq + 1); fe.resize(mq + 1); xs.resize(mq + 1); xe.resize(mq + 1);
            rz::move_bwt::interval I = M.full();
            bool ok = true;
            for (u64 j = mq; j-- > 0;) {
                if (!M.step(I, (unsigned char)Q[j])) { ok = false; break; }
                fs[j] = M.at(I.ks, I.os); fe[j] = M.at(I.ke, I.oe) + 1;
            }
            if (!ok) continue;
            I = M.full();
            bool alive = true;
            for (u64 j = 1; j <= mq; ++j) {
                if (alive && M.step(I, rz::comp((unsigned char)Q[j - 1]))) { xs[j] = M.at(I.ks, I.os); xe[j] = M.at(I.ke, I.oe) + 1; }
                else { alive = false; xs[j] = xe[j] = 0; }
            }
            sink += fs[0] + xs[mq];
        }
        tt = rz::now() - t0;
        T.bs = tt;
        if (sink == 42) printf(" ");
    } else if (isrz) {
        rz::index Z;
        { std::ifstream in(ifile, std::ios::binary); Z.load(in); }
        rz::aux_index A;
        if (!auxfile.empty()) {
            { std::ifstream in(auxfile, std::ios::binary); A.load(in); }
            Z.aux = &A; Z.scanT = scanT;
            // size: the unions replace the four Bf/Br bitvectors; the scan arrays are extra
            struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
                std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } } c1, c2;
            std::ostream o1(&c1), o2(&c2);
            Z.L.Bf.serialize(o2); Z.R.Bf.serialize(o2); Z.L.Br.serialize(o2); Z.R.Br.serialize(o2);
            A.serialize(o1);
            bytes = bytes - c2.c + c1.c;
        }
        rz::rlcsa_bwt X;
        if (!csafile.empty()) {
            X.load(csafile); Z.csa = &X; X_ef = X.ef;
            struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
                std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } } cnt;
            std::ostream os(&cnt);
            Z.bwt.serialize(os);
            bytes = bytes - cnt.c + X.bytes();
        }
        if (!mvfile.empty()) {
            rz::move_bwt M;
            { std::ifstream in(mvfile, std::ios::binary); M.load(in); }
            struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
                std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } } cnt;
            std::ostream os(&cnt);
            Z.bwt.serialize(os);
            bytes = bytes - cnt.c + M.bytes();
            t0 = rz::now();
            for (u64 i = 0; i < NP; ++i) { ans[i] = Z.query(P[i], M, &T); lq += Z.last_left_queries; rq += Z.last_right_queries; }
            tt = rz::now() - t0;
        } else {
            t0 = rz::now();
            for (u64 i = 0; i < NP; ++i) { ans[i] = Z.query(P[i], &T); lq += Z.last_left_queries; rq += Z.last_right_queries; }
            tt = rz::now() - t0;
        }
    } else if (issr && !rlfile.empty()) {
        rz::rlbwt B;
        { std::ifstream in(rlfile, std::ios::binary); u64 nn; in.read((char *)&nn, 8); B.load(in); }
        if (!csafile.empty()) {   // RLCSA for backward search and LF ranks; RLBWT run starts and heads kept
            rz::rlcsa_bwt CS; CS.load(csafile); X_ef = CS.ef;
            rz::csa_nav NV(B, CS);
            rz::srindex_csa X;
            { std::ifstream in(ifile, std::ios::binary); X.load(in); }
            X.mv = &NV;
            struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
                std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } } cnt;
            std::ostream os(&cnt);
            B.starts.serialize(os); B.heads.serialize(os);
            bytes += cnt.c + CS.bytes();
            t0 = rz::now();
            for (u64 i = 0; i < NP; ++i) { u64 o; ans[i] = X.query(P[i], o, &T); occ += o; maxocc = std::max(maxocc, o); }
            tt = rz::now() - t0;
            srs = X.s; srlf = X.lf_steps; srinv = X.invalid_phis;
        } else {
        rz::rl_nav NV(B);
        rz::srindex_rl X;
        { std::ifstream in(ifile, std::ios::binary); X.load(in); }
        X.mv = &NV;
        struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
            std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } } cnt;
        std::ostream os(&cnt);
        B.serialize(os);
        bytes += cnt.c;
        t0 = rz::now();
        for (u64 i = 0; i < NP; ++i) { u64 o; ans[i] = X.query(P[i], o, &T); occ += o; maxocc = std::max(maxocc, o); }
        tt = rz::now() - t0;
        srs = X.s; srlf = X.lf_steps; srinv = X.invalid_phis;
        }
    } else if (issr) {
        rz::move_bwt M;
        { std::ifstream in(mvfile, std::ios::binary); M.load(in); }
        rz::srindex X;
        { std::ifstream in(ifile, std::ios::binary); X.load(in); }
        X.mv = &M;
        bytes += M.bytes();
        t0 = rz::now();
        for (u64 i = 0; i < NP; ++i) { u64 o; ans[i] = X.query(P[i], o, &T); occ += o; maxocc = std::max(maxocc, o); }
        tt = rz::now() - t0;
        srs = X.s; srlf = X.lf_steps; srinv = X.invalid_phis;
    } else {
        rz::rindex R;
        { std::ifstream in(ifile, std::ios::binary); R.load(in); }
        rz::rlcsa_bwt X;
        if (!csafile.empty()) {   // size: the RLCSA and the RLBWT's run starts replace the RLBWT
            X.load(csafile); R.csa = &X; X_ef = X.ef;
            struct counter : std::streambuf { u64 c = 0; int overflow(int ch) override { ++c; return ch; }
                std::streamsize xsputn(const char *, std::streamsize k) override { c += k; return k; } } c1, c2;
            std::ostream o1(&c1), o2(&c2);
            R.bwt.serialize(o1); R.bwt.starts.serialize(o2);
            bytes = bytes - c1.c + X.bytes() + c2.c;
        }
        t0 = rz::now();
        for (u64 i = 0; i < NP; ++i) { u64 o; ans[i] = R.query(P[i], o, &T); occ += o; maxocc = std::max(maxocc, o); }
        tt = rz::now() - t0;
    }
    if (!afile.empty()) {
        std::ofstream o(afile, std::ios::binary);
        for (auto &a : ans) { o.write((char *)&a.first, 8); o.write((char *)&a.second, 8); }
    }
    u64 bad = 0;
    if (!sfile.empty()) {
        std::ifstream in(sfile, std::ios::binary);
        std::stringstream ss; ss << in.rdbuf(); string S = ss.str();
        for (u64 i = 0; i < NP; ++i) {
            size_t a = S.find(P[i]), b = S.rfind(P[i]);
            std::pair<u64, u64> e = {a == string::npos ? rz::NONE : a, b == string::npos ? rz::NONE : b};
            if (e != ans[i]) ++bad;
        }
    }
    double us = 1e6 / NP;
    if (ismv)
        printf("m=%lu patterns=%lu | move structure, rz backward searches only: bytes=%lu bs=%.2f us", m, NP, bytes, tt * us);
    else if (isrz)
        printf("m=%lu patterns=%lu | rz%s: bytes=%lu total=%.2f bs=%.2f left=%.2f right=%.2f us "
               "(grid queries/pattern: left %.1f right %.1f)", m, NP,
               (std::string(mvfile.empty() ? "" : "+mv") + (csafile.empty() ? "" : (X_ef ? "+csaEF" : "+csa")) + (auxfile.empty() ? "" : "+aux")).c_str(), bytes, tt * us, T.bs * us, T.left * us,
               T.right * us, (double)lq / NP, (double)rq / NP);
    else if (issr)
        printf("m=%lu patterns=%lu | sr-index(s=%lu)+%s%s: bytes=%lu total=%.2f bs=%.2f locate=%.2f us avg_occ=%.1f max_occ=%lu "
               "(LF steps/pattern %.1f, invalid phi/pattern %.2f)", m, NP, srs, rlfile.empty() ? "mv" : "rl", csafile.empty() ? "" : (X_ef ? "+csaEF" : "+csa"), bytes, tt * us, T.bs * us, T.left * us,
               (double)occ / NP, maxocc, (double)srlf / NP, (double)srinv / NP);
    else
        printf("m=%lu patterns=%lu | r-index%s: bytes=%lu total=%.2f bs=%.2f phi=%.2f us avg_occ=%.1f max_occ=%lu",
               m, NP, csafile.empty() ? "" : (X_ef ? "+csaEF" : "+csa"), bytes, tt * us, T.bs * us, T.left * us, (double)occ / NP, maxocc);
    if (!dgopt.empty()) printf(" | digest %s avg_len=%.2f", dgopt.c_str(), dglen);
    if (!sfile.empty()) printf(" | brute_force_errors=%lu", bad);
    printf("\n");
    return bad ? 2 : 0;
}
