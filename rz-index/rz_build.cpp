// Build the rz-index (and optionally the baseline r-index) from S, a BWT, and the
// left-to-right / right-to-left parses.
//
// usage: rz-build [-s step] [-a] <S> <bwt> <left.ends> <right.ends> <S.tbl> <out-prefix>
//   <bwt>        one byte per row.  Either
//                  * BWT(S$) of S as a single string (e.g. Big-BWT; one terminator), or
//                  * a multi-string BWT / eBWT with one terminator per dataset, where the
//                    datasets are the species (maximal runs of consecutive genomes with the
//                    same species in <S.tbl>), i.e. the files rz-prep -d writes (e.g. pfp-merge).
//                Terminator bytes may be 0x00, 0x01 or 0x02.  The number of terminators tells
//                which case applies.  The order of the datasets' terminator rows does not matter.
//   <left.ends>  phrase ends of a left-referencing parse of S      (rz-lz77)
//   <right.ends> phrase ends of a left-referencing parse of S^rev  (rz-lz77 -r)
//   <S.tbl>      genome table written by rz-prep
//   -s step      store grid values every `step` wavelet-matrix levels (default: bottom only)
//   -a           also build the baseline r-index on the same BWT: <out>.rix
//   -A           build only the baseline r-index (<out>.rix); the parse arguments are ignored
//                (pass "-").  Building the two indexes in separate runs lowers peak memory.
//   -G           no grids: <out>.rz holds only the RLBWT and the strand starts, which is all that
//                listing with tag arrays or sr-indexes needs (not the rz-index's own LCA queries);
//                the parse arguments are ignored (pass "-"), so no LZ77 parses are needed
// writes <out>.rz
//
// No suffix array is built.  The BWT rows of the O(z) positions we need (suffixes after
// phrase boundaries and mirrored positions) are found by one LF pass over each string,
// using run-length ("move"-style) LF arrays of O(r) words.  The same pass yields the SA
// samples at BWT run boundaries that the r-index needs.

#include "rz_index.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <unistd.h>
#include <numeric>

using namespace rz;
using std::string;
using std::vector;

static vector<u64> read_u64s(const string &f) {
    FILE *fp = fopen(f.c_str(), "rb");
    if (!fp) { std::cerr << "cannot open " << f << "\n"; exit(1); }
    fseeko(fp, 0, SEEK_END);
    u64 sz = ftello(fp) / 8;
    fseeko(fp, 0, SEEK_SET);
    vector<u64> v(sz);
    if (fread(v.data(), 8, sz, fp) != sz) { std::cerr << "short read " << f << "\n"; exit(1); }
    fclose(fp);
    return v;
}

int main(int argc, char **argv) {
    u64 step = 0;
    bool baseline = false, sides = true, nogrids = false;
    int opt;
    const char *usage = "usage: rz-build [-s step] [-a|-A|-G] <S> <bwt> <left.ends> <right.ends> <S.tbl> <out-prefix>\n";
    while ((opt = getopt(argc, argv, "s:aAG")) != -1) {
        if (opt == 's') step = std::stoull(optarg);
        else if (opt == 'a') baseline = true;
        else if (opt == 'A') { baseline = true; sides = false; }
        else if (opt == 'G') nogrids = true;
        else { std::cerr << usage; return 1; }
    }
    if (argc - optind != 6) { std::cerr << usage; return 1; }
    string Sfile = argv[optind], bwtfile = argv[optind + 1], lfile = argv[optind + 2],
           rfile = argv[optind + 3], tfile = argv[optind + 4], out = argv[optind + 5];
    double t0 = now();

    rz::index idx;
    // ---- text length, strands, species ------------------------------------------------
    u64 N;
    FILE *Sfp = fopen(Sfile.c_str(), "rb");
    if (!Sfp) { std::cerr << "cannot open " << Sfile << "\n"; return 1; }
    fseeko(Sfp, 0, SEEK_END); N = ftello(Sfp);
    idx.N = N;
    vector<u64> sstart, slen;           // strands: G_0, rc(G_0), G_1, rc(G_1), ...
    vector<u64> spstart;                // S offsets where a new species begins (+ N at the end)
    {
        std::ifstream t(tfile);
        string line, last; u64 pos = 0;
        while (std::getline(t, line)) {
            std::istringstream ss(line);
            string a, sp, path; u64 len, st;
            std::getline(ss, a, '\t'); std::getline(ss, sp, '\t'); std::getline(ss, path, '\t');
            ss >> len >> st;
            if (st != pos) { std::cerr << "table/offset mismatch\n"; return 1; }
            if (spstart.empty() || sp != last) { spstart.push_back(pos); last = sp; }
            sstart.push_back(pos); slen.push_back(len); pos += len + 1;
            sstart.push_back(pos); slen.push_back(len); pos += len + 1;
        }
        if (pos != N) { std::cerr << "table does not match |S| (" << pos << " vs " << N << ")\n"; return 1; }
    }
    idx.nstr = sstart.size();
    idx.Bs.build(sstart, N);
    auto mirror = [&](u64 q) -> u64 {
        u64 t = std::upper_bound(sstart.begin(), sstart.end(), q) - sstart.begin() - 1;
        u64 o = q - sstart[t];
        if (o >= slen[t]) return q;                       // separator X
        return sstart[t ^ 1] + slen[t] - 1 - o;
    };

    // ---- RLBWT and datasets -----------------------------------------------------------------
    std::cerr << "[1/6] run-length encoding the BWT\n";
    idx.bwt.build(bwtfile);
    const rlbwt &B = idx.bwt;
    u64 R = B.R, n = B.n;
    u64 K = B.C[2] - B.C[1];                             // number of terminators = strings in the BWT
    vector<u64> doff;                                     // dataset starts in S
    if (K == 1) doff = {0};
    else if (K == spstart.size()) doff = spstart;
    else {
        std::cerr << "the BWT has " << K << " terminators, but S has 1 text and " << spstart.size()
                  << " species datasets\n";
        return 1;
    }
    doff.push_back(N);
    if (n != N + K) { std::cerr << "BWT length " << n << " != |S| + " << K << "\n"; return 1; }
    std::cerr << "      n=" << n << " r=" << R << " strings=" << K << "\n";
    if (nogrids) {                              // RLBWT and strands only
        std::ofstream o(out + ".rz", std::ios::binary);
        o.write((char *)&idx.N, 8); o.write((char *)&idx.nstr, 8);
        u64 bytes = 16 + idx.bwt.serialize(o) + idx.Bs.serialize(o);
        o.close(); fclose(Sfp);
        std::cout << "N=" << N << " genomes=" << idx.nstr / 2 << " r=" << R << " rz-index without grids, bytes=" << bytes
                  << "\nbuild time " << now() - t0 << "s\n";
        return 0;
    }

    // ---- parses -----------------------------------------------------------------------------
    std::cerr << "[2/6] reading parses\n";
    vector<u64> eL, bR;
    if (sides) {
        eL = read_u64s(lfile);
        bR = read_u64s(rfile);                  // phrase ends in S^rev -> phrase starts in S
        for (auto &x : bR) x = N - 1 - x;
        std::reverse(bR.begin(), bR.end());
    }
    u64 zL = eL.size(), zR = bR.size();
    if (sides && (!zL || eL.back() != N - 1 || !zR || bR[0] != 0)) { std::cerr << "parse files do not cover S\n"; return 1; }
    if (2 * zL + 2 * zR >= (1ULL << 32)) { std::cerr << "too many phrases\n"; return 1; }
    idx.L.z = zL; idx.R.z = zR;

    // S positions whose BWT rows we need.  Position q < N is the rotation of its dataset
    // starting at q; position N is the terminator rotation of the last dataset.
    // The needs are generated one dataset at a time (a mirror stays inside its dataset), so only
    // the rows themselves, at ceil(log n) bits each, are kept for the whole text.
    struct __attribute__((packed)) need { u64 pos; uint32_t slot; };
    int_vector<> rows(sides ? 2 * zL + 2 * zR : 0, 0, bits::hi(N + 64) + 1);
    if (sides) rows[2 * zL + zR + 0] = 0;       // b = 0: empty prefix -> row 0 (a terminator row)

    auto needs_of = [&](u64 j, vector<need> &v) {           // needed positions q of dataset j, descending
        v.clear();
        if (!sides) return;
        u64 lo = doff[j], hi = doff[j + 1] + (j + 1 == doff.size() - 1 ? 1 : 0);
        auto in = [&](u64 q) { return q >= lo && q < hi; };
        u64 a = std::lower_bound(eL.begin(), eL.end(), lo == 0 ? 0 : lo - 1) - eL.begin();
        for (u64 k = a; k < zL && eL[k] < hi; ++k) {
            if (in(eL[k] + 1)) v.push_back({eL[k] + 1, (uint32_t)k});
            u64 m = mirror(eL[k]);
            if (in(m)) v.push_back({m, (uint32_t)(zL + k)});
        }
        a = std::lower_bound(bR.begin(), bR.end(), lo) - bR.begin();
        for (u64 k = a; k < zR && bR[k] <= hi; ++k) {
            if (in(bR[k])) v.push_back({bR[k], (uint32_t)(2 * zL + k)});
            if (bR[k] > 0) { u64 m = mirror(bR[k] - 1); if (in(m)) v.push_back({m, (uint32_t)(2 * zL + zR + k)}); }
        }
        std::sort(v.begin(), v.end(), [](const need &x, const need &y) { return x.pos > y.pos; });
    };

    // ---- LF arrays (move-style) ---------------------------------------------------------------
    std::cerr << "[3/6] LF arrays\n";
    struct __attribute__((packed)) mrun {                 // 14 bytes per run; rows < 2^40
        uint32_t s0, l0, lfr; uint8_t s1, l1;
        u64 start() const { return s0 | ((u64)s1 << 32); }
        u64 lfs() const { return l0 | ((u64)l1 << 32); }
        void set_start(u64 v) { s0 = (uint32_t)v; s1 = (uint8_t)(v >> 32); }
        void set_lfs(u64 v) { l0 = (uint32_t)v; l1 = (uint8_t)(v >> 32); }
    };
    if (n >= (1ULL << 40)) { std::cerr << "BWT longer than 2^40\n"; return 1; }
    if (R >= (1ULL << 32) - 1) { std::cerr << "too many runs\n"; return 1; }
    vector<mrun> mv(R + 1);
    vector<unsigned char> rhead(R);
    {
        u64 before[256] = {0};
        rlbwt::runs(bwtfile, [&](u64 k, u64 s, u64 len, unsigned char c) {
            mv[k].set_start(s); rhead[k] = c;
            mv[k].set_lfs(B.C[c] + before[c]);
            before[c] += len;
        });
        mv[R].set_start(n);
    }
    auto run_containing = [&](u64 row) {
        return (u64)(std::upper_bound(mv.begin(), mv.end(), row, [](u64 v, const mrun &x) { return v < x.start(); }) - mv.begin() - 1);
    };
    for (u64 k = 0; k < R; ++k) mv[k].lfr = (uint32_t)run_containing(mv[k].lfs());
    auto LF = [&](u64 &row, u64 &k) {
        const mrun &cur = mv[k];
        row = cur.lfs() + (row - cur.start());
        k = cur.lfr;
        while (mv[k + 1].start() <= row) ++k;
    };

    // ---- which terminator row belongs to which dataset ------------------------------------------
    // Row x < K is the rotation starting at some dataset's terminator; reading the BWT backwards
    // from x spells that dataset's text backwards, which we match against the datasets' tails.
    std::cerr << "[4/6] matching terminator rows to datasets\n";
    vector<u64> drow(K, NONE);
    if (K == 1) drow[0] = 0;
    else {
        std::map<string, u64> tails;
        u64 LEN = 64;
        while (true) {
            tails.clear();
            bool uniq = true;
            for (u64 j = 0; j < K; ++j) {
                u64 len = std::min(LEN, doff[j + 1] - doff[j]);
                string s(len, 0);
                fseeko(Sfp, doff[j + 1] - len, SEEK_SET);
                if (fread(&s[0], 1, len, Sfp) != len) return 1;
                std::reverse(s.begin(), s.end());
                if (!tails.emplace(s, j).second) uniq = false;
            }
            if (uniq) break;
            if (LEN > (1ULL << 26)) { std::cerr << "two datasets have identical tails\n"; return 1; }
            LEN *= 4;
        }
        for (u64 x = 0; x < K; ++x) {
            string s;
            u64 row = x, k = run_containing(x);
            while (s.size() < LEN) {
                unsigned char c = rhead[k];
                if (c == 1) break;
                s.push_back(c);
                LF(row, k);
            }
            // the dataset whose (reversed) tail equals s, or which s exhausts (short dataset)
            auto it = tails.find(s);
            if (it == tails.end()) { std::cerr << "terminator row " << x << " matches no dataset\n"; return 1; }
            if (drow[it->second] != NONE) { std::cerr << "two terminator rows match dataset " << it->second << "\n"; return 1; }
            drow[it->second] = x;
        }
    }

    // ---- one LF pass per string: rows of needed positions, SA samples ------------------------------
    std::cerr << "[5/6] LF pass over the text\n";
    int_vector<> ssa;                      // SA in D coordinates at run starts (loaded after the pass)
    vector<std::pair<u64, u64>> term;      // (row, D position) of rows whose BWT char is a terminator
    // SA samples at run starts / ends are streamed to temporary files during the LF pass and
    // loaded only after the run arrays are freed, to keep them out of the pass's peak memory.
    struct __attribute__((packed)) samp { uint32_t k; u64 v; };
    FILE *fss = nullptr, *fes = nullptr;
    string tss = out + ".ssa.tmp", tes = out + ".esa.tmp";
    if (baseline) {
        fss = fopen(tss.c_str(), "wb"); fes = fopen(tes.c_str(), "wb");
        if (!fss || !fes) { std::cerr << "cannot write temporary sample files\n"; return 1; }
        setvbuf(fss, nullptr, _IOFBF, 1 << 24); setvbuf(fes, nullptr, _IOFBF, 1 << 24);
    }
    auto load_samples = [&](const string &f, int_vector<> &a) {
        a = int_vector<>(R, 0, bits::hi(n) + 1);
        FILE *fp = fopen(f.c_str(), "rb");
        vector<samp> buf(1 << 20);
        size_t got, tot = 0;
        while ((got = fread(buf.data(), sizeof(samp), buf.size(), fp)) > 0) {
            for (size_t i = 0; i < got; ++i) a[buf[i].k] = buf[i].v;
            tot += got;
        }
        fclose(fp); unlink(f.c_str());
        if (tot != R) { std::cerr << "internal error: " << tot << " samples in " << f << ", expected " << R << "\n"; exit(1); }
    };
    vector<need> nd;
    u64 filled = 0;
    {
        u64 steps = 0;
        double tl = now();
        for (u64 j = K; j-- > 0;) {
            u64 len = doff[j + 1] - doff[j], dpos = doff[j] + j;    // D offset of dataset j
            u64 row = drow[j], k = run_containing(row), p = len, ptr = 0;
            needs_of(j, nd);
            while (true) {
                u64 q = doff[j] + p;                                  // S position (q == N only for p == len, j == K-1)
                if (p < len || j == K - 1)
                    while (ptr < nd.size() && nd[ptr].pos == q) rows[nd[ptr++].slot] = row;
                if (baseline) {
                    if (row == mv[k].start()) { samp e{(uint32_t)k, dpos + p}; fwrite(&e, sizeof e, 1, fss); }
                    if (row + 1 == mv[k + 1].start()) { samp e{(uint32_t)k, dpos + p}; fwrite(&e, sizeof e, 1, fes); }
                    if (p == 0) term.push_back({row, dpos});
                }
                if (p == 0) break;
                LF(row, k);
                --p;
                if ((++steps & ((1ULL << 30) - 1)) == 0)
                    std::cerr << "\r  " << steps << " / " << N << " (" << steps / (now() - tl) / 1e6 << " M steps/s)" << std::flush;
            }
            if (ptr != nd.size()) { std::cerr << "internal error: missing rows in dataset " << j << "\n"; return 1; }
            filled += nd.size();
            if (rhead[k] != 1) { std::cerr << "LF pass over dataset " << j << " did not end at a terminator: BWT inconsistent\n"; return 1; }
        }
        std::cerr << "\r  " << N << " / " << N << " (" << N / (now() - tl + 1e-9) / 1e6 << " M steps/s)\n";
        if (sides && filled + 1 != rows.size()) { std::cerr << "internal error: " << filled << " rows found, " << rows.size() - 1 << " needed\n"; return 1; }
    }
    vector<need>().swap(nd);
    if (baseline) { fclose(fss); fclose(fes); }

    auto sz = [](auto &x) { std::ostringstream s; return (u64)x.serialize(s); };
    if (baseline) {
        std::cerr << "[6/6] baseline r-index\n";
        rz::rindex ri;
        ri.n = n;
        u64 w = bits::hi(n) + 1;
        // phi samples: run starts x >= 1 (SA[x-1] = SA at the end of run k-1) and terminator rows
        std::sort(term.begin(), term.end());
        vector<std::pair<u64, u64>> xk;          // extra (key, index of SA[x-1]) for terminator rows
        vector<u64> ext;
        for (u64 t = 0; t < term.size(); ++t) {
            u64 x = term[t].first;
            if (x == 0) continue;
            u64 k = run_containing(x);
            if (mv[k].start() == x) continue;                              // already a run start
            if (t == 0 || term[t - 1].first != x - 1) { std::cerr << "internal error: terminator run\n"; return 1; }
            xk.push_back({term[t].second, R + ext.size()});
            ext.push_back(term[t - 1].second);
        }
        vector<mrun>().swap(mv); vector<unsigned char>().swap(rhead);
        load_samples(tss, ssa);
        u64 M = (R - 1) + xk.size();
        auto key = [&](u64 i) -> u64 { return i < R - 1 ? (u64)ssa[i + 1] : xk[i - (R - 1)].first; };
        vector<uint32_t> perm(M);
        std::iota(perm.begin(), perm.end(), 0);
        std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b) { return key(a) < key(b); });
        {
            sd_vector_builder kb(n, M);
            ri.vals = int_vector<>(M, 0, bits::hi(R + ext.size()) + 1);
            for (u64 i = 0; i < M; ++i) {
                u64 q = perm[i];
                kb.set(key(q));
                ri.vals[i] = q < R - 1 ? q : xk[q - (R - 1)].second;
            }
            vector<uint32_t>().swap(perm);
            int_vector<>().swap(ssa);
            ri.keys.build_from(kb);
        }
        load_samples(tes, ri.esa);
        ri.extra = int_vector<>(ext.size(), 0, w);
        for (u64 i = 0; i < ext.size(); ++i) ri.extra[i] = ext[i];
        vector<u64> ds(K);
        for (u64 j = 0; j < K; ++j) ds[j] = doff[j] + j;
        ri.dstart.build(ds, n);
        // same format as rindex::serialize, but writing the RLBWT already in memory
        std::ofstream ro(out + ".rix", std::ios::binary);
        ro.write((char *)&ri.n, 8);
        u64 bb = idx.bwt.serialize(ro);
        u64 rb = 8 + bb + ri.esa.serialize(ro) + ri.keys.serialize(ro) + ri.vals.serialize(ro)
                 + ri.extra.serialize(ro) + ri.dstart.serialize(ro);
        std::cout << "r-index bytes=" << rb << " rlbwt=" << bb << " samples=" << rb - bb << "\n";
    }
    vector<mrun>().swap(mv); vector<unsigned char>().swap(rhead);

    if (sides) {
    // ---- sides ---------------------------------------------------------------------------------
    std::cerr << "[7/7] bitvectors and grids\n";
    auto build_side = [&](side &sd, u64 z, u64 oS, u64 oM, const vector<u64> &pos, bool maxside) {
        vector<u64> a(z); for (u64 k = 0; k < z; ++k) a[k] = rows[oS + k];
        std::sort(a.begin(), a.end());
        sd.Bf.build(a, n);
        vector<u64> b(z); for (u64 k = 0; k < z; ++k) b[k] = rows[oM + k];
        std::sort(b.begin(), b.end());
        for (u64 i = 1; i < z; ++i) if (b[i] == b[i - 1] || a[i] == a[i - 1]) { std::cerr << "internal error: duplicate row\n"; exit(1); }
        sd.Br.build(b, n);
        sd.Bp.build(pos, N);
        vector<u64> Y(z), Kv(z);
        for (u64 k = 0; k < z; ++k) {
            u64 x = sd.Br.rank(rows[oM + k]), y = sd.Bf.rank(rows[oS + k]);
            Y[x] = y; Kv[x] = maxside ? z - 1 - k : k;
        }
        sd.grid.build(Y, Kv, step);
    };
    build_side(idx.L, zL, 0, zL, eL, false);
    build_side(idx.R, zR, 2 * zL, 2 * zL + zR, bR, true);

    std::ofstream o(out + ".rz", std::ios::binary);
    u64 bytes = idx.serialize(o);
    o.close();
    std::cout << "N=" << N << " genomes=" << idx.nstr / 2 << " strings_in_bwt=" << K << " r=" << R
              << " zL=" << zL << " zR=" << zR << "\n";
    std::cout << "rz-index bytes=" << bytes << " rlbwt=" << sz(idx.bwt) << " strands=" << sz(idx.Bs)
              << " left=" << sz(idx.L) << " (grid " << sz(idx.L.grid) << ") right=" << sz(idx.R)
              << " (grid " << sz(idx.R.grid) << ")\n";
    }
    fclose(Sfp);
    std::cout << "build time " << now() - t0 << "s\n";
    return 0;
}
