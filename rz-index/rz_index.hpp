// rz-index: leftmost and rightmost occurrences of a pattern in a collection of
// DNA strings closed under reverse complement, laid out as
//     S = G_1 X rc(G_1) X G_2 X rc(G_2) X ... G_k X rc(G_k) X
// (X is self-complementary and never occurs in patterns).
//
// Components
//   bwt            : run-length BWT of S$ (the same RLBWT an r-index for S uses)
//   Bs             : sparse bitvector marking string starts in S (string t = G_{t/2} or its rc)
//   left  grid     : from a left-referencing LZ-parse of S (every phrase has an earlier source).
//                    For each phrase end e:  y = rank of row(S[e+1..]) in BfL,
//                    x = rank of row(mirror(e)) in BrL, value = k (text order of e).
//                    The leftmost occurrence of P contains a phrase end.
//   right grid     : from a left-referencing parse of S^rev, i.e. a right-referencing parse of S.
//                    For each phrase start b: y = rank of row(S[b..]) in BfR,
//                    x = rank of row(mirror(b-1)) in BrR (row 0 if b = 0), value = zR-1-k.
//                    The rightmost occurrence of P contains a phrase start.
//   mirror(q)      : if q is in string t at offset o, the position of offset |t|-1-o in the
//                    partner string t^1, so that S[..q] ends with Q iff the suffix at mirror(q)
//                    starts with rc(Q) (for Q without X);  if S[q] = X, mirror(q) = q.
//
// Query:  backward search P (intervals of P[i..m)) and rc(P) (intervals of rc(P[0..i))),
//         then, for each split, map intervals through the sparse bitvectors and run a
//         2D range-min on each grid.
//
// All positions 0-based; the terminator $ is 0x01 inside the RLBWT.

#pragma once
#include <sdsl/sd_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/construct.hpp>
#include <vector>
#include "mv.hpp"
#include "digest.hpp"
#include "rlcsa.hpp"
#include <string>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <memory>

namespace rz {

using namespace sdsl;
typedef uint64_t u64;
static const u64 NONE = UINT64_MAX;

inline unsigned char comp(unsigned char c) {
    return rc_symbol(c);                       // DNA and minimizer-digest symbols (digest.hpp)
}

// ---------------------------------------------------------------------------
// Sparse bitvector (Elias-Fano).  Build in place; supports hold pointers.
// ---------------------------------------------------------------------------
struct sparse_bv {
    sd_vector<> bv;
    rank_support_sd<1> rk;
    select_support_sd<1> sl;
    sparse_bv() {}
    sparse_bv(const sparse_bv &) = delete;
    sparse_bv &operator=(const sparse_bv &) = delete;
    template <class V> void build(const V &ones, u64 len) {   // sorted, distinct
        sd_vector_builder b(len, ones.size());
        for (u64 x : ones) b.set(x);
        bv = sd_vector<>(b);
        util::init_support(rk, &bv);
        util::init_support(sl, &bv);
    }
    void build_from(sd_vector_builder &b) {
        bv = sd_vector<>(b);
        util::init_support(rk, &bv);
        util::init_support(sl, &bv);
    }
    inline u64 rank(u64 i) const { return rk(i); }       // # ones in [0,i)
    inline u64 select(u64 k) const { return sl(k); }     // k-th one, k >= 1
    u64 ones() const { return rk(bv.size()); }
    u64 serialize(std::ostream &out) const { return bv.serialize(out) + rk.serialize(out) + sl.serialize(out); }
    void load(std::istream &in) { bv.load(in); rk.load(in, &bv); sl.load(in, &bv); }
};

// ---------------------------------------------------------------------------
// Run-length BWT with rank (Mäkinen-Navarro style), built by streaming a BWT file.
//   starts : sd_vector over [0,n) marking run starts           (R ones)
//   heads  : Huffman-shaped wavelet tree over the R run heads
//   cum[c] : sd_vector over [0, count_c) marking the last position of each c-run
//            in the concatenation of all c characters             (R_c ones)
// ---------------------------------------------------------------------------
struct rlbwt {
    u64 n = 0, R = 0;
    u64 C[257];
    sparse_bv starts;
    wt_huff<> heads;
    std::vector<sparse_bv> cum = std::vector<sparse_bv>(256);

    rlbwt() { std::fill(C, C + 257, 0); }
    rlbwt(const rlbwt &) = delete;
    rlbwt &operator=(const rlbwt &) = delete;

    // BWT file, one byte per row; terminator bytes 0x00-0x02 (one per string) are mapped to 0x01
    template <class F> static void scan(const std::string &file, F f) {
        FILE *fp = fopen(file.c_str(), "rb");
        if (!fp) { std::cerr << "cannot open " << file << "\n"; exit(1); }
        std::vector<unsigned char> buf(1 << 24);
        size_t got; u64 pos = 0;
        while ((got = fread(buf.data(), 1, buf.size(), fp)) > 0) {
            for (size_t i = 0; i < got; ++i) { unsigned char c = buf[i] > 2 ? buf[i] : 1; f(pos++, c); }
        }
        fclose(fp);
    }

    // calls g(run_index, start, length, head) for every run of the BWT file
    template <class G> static void runs(const std::string &file, G g) {
        u64 k = 0, s = 0; int last = -1;
        u64 len = 0;
        scan(file, [&](u64 pos, unsigned char c) {
            if ((int)c != last) {
                if (last >= 0) g(k++, s, len, (unsigned char)last);
                last = c; s = pos; len = 0;
            }
            ++len;
        });
        if (last >= 0) g(k++, s, len, (unsigned char)last);
    }

    void build(const std::string &file) {
        u64 cnt[256] = {0}, rc[256] = {0};
        runs(file, [&](u64 k, u64 s, u64 len, unsigned char c) { cnt[c] += len; rc[c]++; R = k + 1; n = s + len; });
        C[0] = 0;
        for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + cnt[c];
        sd_vector_builder sb(n, R);
        std::vector<std::unique_ptr<sd_vector_builder>> cb(256);
        for (int c = 0; c < 256; ++c) cb[c].reset(new sd_vector_builder(cnt[c], rc[c]));
        int_vector<8> hv(R);
        u64 acc[256] = {0};
        runs(file, [&](u64 k, u64 s, u64 len, unsigned char c) {
            sb.set(s);
            acc[c] += len;
            cb[c]->set(acc[c] - 1);
            hv[k] = c;
        });
        starts.build_from(sb);
        for (int c = 0; c < 256; ++c) cum[c].build_from(*cb[c]);
        construct_im(heads, hv);
    }

    // # of c in [0, i)
    inline u64 rank(u64 i, unsigned char c) const {
        if (i == 0 || C[c + 1] == C[c]) return 0;
        if (i >= n) return C[c + 1] - C[c];
        u64 q = starts.rank(i) - 1;                 // run containing i-1
        u64 s = starts.select(q + 1);
        auto is = heads.inverse_select(q);          // (rank of heads[q] in [0,q), heads[q])
        u64 t = is.second == c ? is.first : heads.rank(q, c);
        u64 base = t == 0 ? 0 : cum[c].select(t) + 1;
        if (is.second == c) base += i - s;
        return base;
    }
    inline u64 run_of(u64 i) const { return starts.rank(i + 1) - 1; }
    inline unsigned char at(u64 i) const { return heads[run_of(i)]; }
    // position of the j-th (0-based) occurrence of c
    inline u64 select(u64 j, unsigned char c) const {
        u64 tt = cum[c].rank(j);                              // index of the c-run containing it
        u64 off = j - (tt ? cum[c].select(tt) + 1 : 0);
        u64 q = heads.select(tt + 1, c);
        return starts.select(q + 1) + off;
    }
    // [sp, ep) -> interval of cX
    inline void extend(u64 &sp, u64 &ep, unsigned char c) const {
        sp = C[c] + rank(sp, c);
        ep = C[c] + rank(ep, c);
    }
    u64 serialize(std::ostream &out) const {
        u64 w = 0;
        out.write((char *)&n, 8); out.write((char *)&R, 8); out.write((char *)C, 257 * 8); w += 16 + 257 * 8;
        w += starts.serialize(out);
        w += heads.serialize(out);
        for (int c = 0; c < 256; ++c) w += cum[c].serialize(out);
        return w;
    }
    void load(std::istream &in) {
        in.read((char *)&n, 8); in.read((char *)&R, 8); in.read((char *)C, 257 * 8);
        starts.load(in);
        heads.load(in);
        for (int c = 0; c < 256; ++c) cum[c].load(in);
    }
};

// ---------------------------------------------------------------------------
// 2D range-min over a permutation grid: wavelet matrix + succinct RMQ per level.
// Values are recovered by descending to the next level that stores them explicitly.
// ---------------------------------------------------------------------------
struct grid_rmq {
    u64 z = 0, L = 0, step = 0;
    std::vector<bit_vector> bv;
    std::vector<rank_support_v5<1>> rk;
    std::vector<u64> Z;
    std::vector<rmq_succinct_sct<true>> rmq;
    std::vector<int_vector<>> val;

    grid_rmq() {}
    grid_rmq(const grid_rmq &) = delete;
    grid_rmq &operator=(const grid_rmq &) = delete;

    void build(std::vector<u64> &Y, std::vector<u64> &K, u64 step_) {
        z = Y.size();
        L = 1;
        while ((1ULL << L) < z) ++L;
        step = step_ == 0 ? L : step_;
        bv.resize(L); rk.resize(L); Z.resize(L); rmq.resize(L); val.resize(L + 1);
        u64 w = std::max<u64>(1, bits::hi(std::max<u64>(1, z - 1)) + 1);
        std::vector<u64> y2(z), k2(z);
        for (u64 l = 0; l <= L; ++l) {
            int_vector<> kv(z, 0, w);
            for (u64 i = 0; i < z; ++i) kv[i] = K[i];
            if (l == L || l % step == 0) val[l] = kv;
            if (l == L) break;
            rmq[l] = rmq_succinct_sct<true>(&kv);
            bv[l] = bit_vector(z, 0);
            u64 bit = L - 1 - l, zeros = 0;
            for (u64 i = 0; i < z; ++i)
                if ((Y[i] >> bit) & 1) bv[l][i] = 1; else ++zeros;
            Z[l] = zeros;
            util::init_support(rk[l], &bv[l]);
            u64 a = 0, b = zeros;
            for (u64 i = 0; i < z; ++i) {
                if ((Y[i] >> bit) & 1) { y2[b] = Y[i]; k2[b] = K[i]; ++b; }
                else { y2[a] = Y[i]; k2[a] = K[i]; ++a; }
            }
            Y.swap(y2); K.swap(k2);
        }
    }
    inline u64 value(u64 l, u64 p) const {
        while (val[l].size() == 0) {
            u64 r1 = rk[l](p);
            p = bv[l][p] ? Z[l] + r1 : p - r1;
            ++l;
        }
        return val[l][p];
    }
    // min value over x in [x1,x2), y in [y1,y2); NONE if empty
    u64 query(u64 x1, u64 x2, u64 y1, u64 y2) const {
        if (x1 >= x2 || y1 >= y2) return NONE;
        u64 best = NONE;
        rec(0, 0, x1, x2, y1, y2, best);
        return best;
    }
    void rec(u64 l, u64 prefix, u64 s, u64 e, u64 y1, u64 y2, u64 &best) const {
        if (s >= e) return;
        u64 lo = prefix << (L - l), hi = (prefix + 1) << (L - l);
        if (hi <= y1 || lo >= y2) return;
        if (y1 <= lo && hi <= y2) {
            u64 p = (l == L || e - s == 1) ? s : rmq[l](s, e - 1);
            u64 v = value(l, p);
            if (v < best) best = v;
            return;
        }
        u64 rs = rk[l](s), re = rk[l](e);
        rec(l + 1, prefix << 1, s - rs, e - re, y1, y2, best);
        rec(l + 1, (prefix << 1) | 1, Z[l] + rs, Z[l] + re, y1, y2, best);
    }
    u64 serialize(std::ostream &out) const {
        u64 w = 24;
        out.write((char *)&z, 8); out.write((char *)&L, 8); out.write((char *)&step, 8);
        for (u64 l = 0; l < L; ++l) {
            w += bv[l].serialize(out) + rk[l].serialize(out);
            out.write((char *)&Z[l], 8); w += 8;
            w += rmq[l].serialize(out);
        }
        for (u64 l = 0; l <= L; ++l) w += val[l].serialize(out);
        return w;
    }
    void load(std::istream &in) {
        in.read((char *)&z, 8); in.read((char *)&L, 8); in.read((char *)&step, 8);
        bv.resize(L); rk.resize(L); Z.resize(L); rmq.resize(L); val.resize(L + 1);
        for (u64 l = 0; l < L; ++l) {
            bv[l].load(in); rk[l].load(in, &bv[l]);
            in.read((char *)&Z[l], 8);
            rmq[l].load(in);
        }
        for (u64 l = 0; l <= L; ++l) val[l].load(in);
    }
};

// one side (left or right) of the index
struct side {
    u64 z = 0;
    sparse_bv Bf;    // rows of suffixes on the "suffix" side of each boundary
    sparse_bv Br;    // rows of mirrored positions ("prefix" side)
    sparse_bv Bp;    // text positions: phrase ends (left) or phrase starts (right)
    grid_rmq grid;
    u64 serialize(std::ostream &out) const {
        out.write((char *)&z, 8);
        return 8 + Bf.serialize(out) + Br.serialize(out) + Bp.serialize(out) + grid.serialize(out);
    }
    void load(std::istream &in) {
        in.read((char *)&z, 8);
        Bf.load(in); Br.load(in); Bp.load(in); grid.load(in);
    }
};

struct timing { double bs = 0, left = 0, right = 0; };

static inline double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Auxiliary structures for faster grid queries (built from an existing index by rz-aux):
//  * union bitvectors: the left and right sides rank the same BWT intervals, so each pair of
//    sparse bitvectors (L.Bf + R.Bf, L.Br + R.Br) is replaced by one Elias-Fano union plus two
//    flag bitvectors, halving the Elias-Fano ranks per split;
//  * scan arrays per grid: y by x, x by y and value by x, so that a query whose narrower side
//    has at most T points is answered by scanning those points instead of walking the wavelet
//    matrix.
struct union_bv {
    sparse_bv U;
    bit_vector fl, fr; rank_support_v5<1> rl, rr;
    inline void rank2(u64 p, u64 &a, u64 &b) const { u64 j = U.rank(p); a = rl(j); b = rr(j); }
    u64 serialize(std::ostream &out) const { return U.serialize(out) + fl.serialize(out) + fr.serialize(out) + rl.serialize(out) + rr.serialize(out); }
    void load(std::istream &in) { U.load(in); fl.load(in); fr.load(in); rl.load(in, &fl); rr.load(in, &fr); }
};
struct scan_arrays {
    int_vector<> yx, xy, kx;                 // y of x, x of y, value of x
    u64 serialize(std::ostream &out) const { return yx.serialize(out) + xy.serialize(out) + kx.serialize(out); }
    void load(std::istream &in) { yx.load(in); xy.load(in); kx.load(in); }
    // min value over x in [x1,x2), y in [y1,y2), scanning the narrower side
    inline u64 query(u64 x1, u64 x2, u64 y1, u64 y2) const {
        u64 best = NONE;
        if (x2 - x1 <= y2 - y1) {
            for (u64 x = x1; x < x2; ++x) { u64 y = yx[x]; if (y >= y1 && y < y2) { u64 v = kx[x]; if (v < best) best = v; } }
        } else {
            for (u64 y = y1; y < y2; ++y) { u64 x = xy[y]; if (x >= x1 && x < x2) { u64 v = kx[x]; if (v < best) best = v; } }
        }
        return best;
    }
};
struct aux_index {
    union_bv F, Rr;                          // F: L.Bf + R.Bf ; Rr: L.Br + R.Br
    scan_arrays SL, SR;
    u64 serialize(std::ostream &out) const { return F.serialize(out) + Rr.serialize(out) + SL.serialize(out) + SR.serialize(out); }
    void load(std::istream &in) { F.load(in); Rr.load(in); SL.load(in); SR.load(in); }
};

struct index {
    u64 N = 0;          // |S|
    u64 nstr = 0;       // number of strings (2 per genome)
    rlbwt bwt;
    sparse_bv Bs;       // string starts in S
    side L, R;

    // scratch
    std::vector<u64> fs, fe, xs, xe;
    u64 last_left_queries = 0, last_right_queries = 0;
    const aux_index *aux = nullptr;       // optional (rz-aux)
    const rlcsa_bwt *csa = nullptr;       // optional: do backward searches on an RLCSA instead of bwt
    inline void ext(u64 &sp, u64 &ep, unsigned char c) const { if (csa) csa->extend(sp, ep, c); else bwt.extend(sp, ep, c); }
    u64 scanT = 16;                       // scan when the narrower side has at most scanT points

    u64 genome_of(u64 pos) const { return (Bs.rank(pos + 1) - 1) / 2; }

    // leftmost and rightmost starting positions of P in S (NONE if P does not occur)
    std::pair<u64, u64> query(const std::string &P, timing *T = nullptr) {
        u64 m = P.size();
        double t0 = T ? now() : 0;
        fs.resize(m + 1); fe.resize(m + 1); xs.resize(m + 1); xe.resize(m + 1);
        u64 sp = 0, ep = bwt.n;
        fs[m] = sp; fe[m] = ep;
        for (u64 i = m; i-- > 0;) {
            ext(sp, ep, (unsigned char)P[i]);
            if (sp >= ep) return {NONE, NONE};
            fs[i] = sp; fe[i] = ep;
        }
        sp = 0; ep = bwt.n;
        xs[0] = 0; xe[0] = bwt.n;
        for (u64 i = 1; i <= m; ++i) {           // rc(P[0..i)) = comp(P[i-1]) rc(P[0..i-1))
            ext(sp, ep, comp((unsigned char)P[i - 1]));
            xs[i] = sp; xe[i] = ep;
        }
        return grids(m, t0, T);
    }
    // same query when the caller has already filled fs/fe[0..m] with the BWT intervals of the suffixes
    // P[i..m) (fs[m] = 0, fe[m] = n): only the search for rc(P) and the grids remain
    std::pair<u64, u64> query_with_suffixes(const std::string &P, timing *T = nullptr) {
        u64 m = P.size();
        double t0 = T ? now() : 0;
        xs.resize(m + 1); xe.resize(m + 1);
        u64 sp = 0, ep = bwt.n;
        xs[0] = 0; xe[0] = bwt.n;
        for (u64 i = 1; i <= m; ++i) {
            ext(sp, ep, comp((unsigned char)P[i - 1]));
            xs[i] = sp; xe[i] = ep;
        }
        return grids(m, t0, T);
    }
    // same query when the caller has filled both fs/fe[0..m] and xs/xe[0..m] (xs/xe[i]: rc(P[0..i)))
    std::pair<u64, u64> query_with_both(u64 m, timing *T = nullptr) {
        double t0 = T ? now() : 0;
        return grids(m, t0, T);
    }
    // same query, with the backward searches done on a move structure
    std::pair<u64, u64> query(const std::string &P, move_bwt &mv, timing *T = nullptr) {
        u64 m = P.size();
        double t0 = T ? now() : 0;
        fs.resize(m + 1); fe.resize(m + 1); xs.resize(m + 1); xe.resize(m + 1);
        move_bwt::interval I = mv.full();
        fs[m] = 0; fe[m] = mv.n;
        for (u64 i = m; i-- > 0;) {
            if (!mv.step(I, (unsigned char)P[i])) return {NONE, NONE};
            fs[i] = mv.at(I.ks, I.os); fe[i] = mv.at(I.ke, I.oe) + 1;
        }
        I = mv.full();
        xs[0] = 0; xe[0] = mv.n;
        bool alive = true;
        for (u64 i = 1; i <= m; ++i) {
            if (alive && mv.step(I, comp((unsigned char)P[i - 1]))) { xs[i] = mv.at(I.ks, I.os); xe[i] = mv.at(I.ke, I.oe) + 1; }
            else { alive = false; xs[i] = xe[i] = 0; }
        }
        return grids(m, t0, T);
    }
    // grid part of a query, given fs/fe (suffixes of P) and xs/xe (rc of prefixes of P)
    std::pair<u64, u64> grids(u64 m, double t0, timing *T) {
        if (aux) return grids_aux(m, t0, T);
        double t1 = T ? now() : 0;
        // left: splits i = 1..m, boundary after P[i-1] is a phrase end e; start = e - i + 1
        u64 best = NONE;
        last_left_queries = last_right_queries = 0;
        for (u64 i = 1; i <= m; ++i) {
            u64 x1 = L.Br.rank(xs[i]), x2 = L.Br.rank(xe[i]);
            if (x1 >= x2) continue;
            u64 y1 = L.Bf.rank(fs[i]), y2 = L.Bf.rank(fe[i]);
            if (y1 >= y2) continue;
            ++last_left_queries;
            u64 k = L.grid.query(x1, x2, y1, y2);
            if (k == NONE) continue;
            u64 s = L.Bp.select(k + 1) + 1 - i;
            if (s < best) best = s;
        }
        double t2 = T ? now() : 0;
        // right: splits i = 0..m-1, P[i] is at a phrase start b; start = b - i
        u64 worst = NONE;
        for (u64 i = 0; i < m; ++i) {
            u64 x1, x2;
            if (i == 0) { x1 = 0; x2 = R.z; }
            else { x1 = R.Br.rank(xs[i]); x2 = R.Br.rank(xe[i]); }
            if (x1 >= x2) continue;
            u64 y1 = R.Bf.rank(fs[i]), y2 = R.Bf.rank(fe[i]);
            if (y1 >= y2) continue;
            ++last_right_queries;
            u64 v = R.grid.query(x1, x2, y1, y2);
            if (v == NONE) continue;
            u64 k = R.z - 1 - v;
            u64 s = R.Bp.select(k + 1) - i;
            if (worst == NONE || s > worst) worst = s;
        }
        if (T) { double t3 = now(); T->bs += t1 - t0; T->left += t2 - t1; T->right += t3 - t2; }
        return {best, worst};
    }

    // grid part with the auxiliary structures; left and right are timed together (T->left)
    std::pair<u64, u64> grids_aux(u64 m, double t0, timing *T) {
        double t1 = T ? now() : 0;
        const aux_index &A = *aux;
        u64 best = NONE, worst = NONE;
        last_left_queries = last_right_queries = 0;
        // x ranks for rc(P[0..i)) and y ranks for P[i..m), for both sides at once
        for (u64 i = 0; i <= m; ++i) {
            u64 lx1 = 0, lx2 = 0, rx1 = 0, rx2 = 0;
            if (i == 0) { rx1 = 0; rx2 = R.z; }
            else {
                A.Rr.rank2(xs[i], lx1, rx1); A.Rr.rank2(xe[i], lx2, rx2);
            }
            bool lok = i >= 1 && lx1 < lx2, rok = i < m && rx1 < rx2;
            if (!lok && !rok) continue;
            u64 ly1, ly2, ry1, ry2;
            A.F.rank2(fs[i], ly1, ry1); A.F.rank2(fe[i], ly2, ry2);
            if (lok && ly1 < ly2) {
                ++last_left_queries;
                u64 k = std::min(lx2 - lx1, ly2 - ly1) <= scanT ? A.SL.query(lx1, lx2, ly1, ly2) : L.grid.query(lx1, lx2, ly1, ly2);
                if (k != NONE) { u64 s = L.Bp.select(k + 1) + 1 - i; if (s < best) best = s; }
            }
            if (rok && ry1 < ry2) {
                ++last_right_queries;
                u64 v = std::min(rx2 - rx1, ry2 - ry1) <= scanT ? A.SR.query(rx1, rx2, ry1, ry2) : R.grid.query(rx1, rx2, ry1, ry2);
                if (v != NONE) {
                    u64 k = R.z - 1 - v;
                    u64 s = R.Bp.select(k + 1) - i;
                    if (worst == NONE || s > worst) worst = s;
                }
            }
        }
        if (T) { double t2 = now(); T->bs += t1 - t0; T->left += t2 - t1; }
        return {best, worst};
    }

    u64 serialize(std::ostream &out) const {
        out.write((char *)&N, 8); out.write((char *)&nstr, 8);
        return 16 + bwt.serialize(out) + Bs.serialize(out) + L.serialize(out) + R.serialize(out);
    }
    void load(std::istream &in) {
        in.read((char *)&N, 8); in.read((char *)&nstr, 8);
        bwt.load(in); Bs.load(in); L.load(in); R.load(in);
    }
};

// ---------------------------------------------------------------------------
// Baseline r-index (Gagie, Navarro, Prezza): backward search with toehold, then phi.
// Works on a single-string BWT or on a multi-string BWT with one terminator per string
// (e.g. pfp-merge's eBWT).  Positions are in "D" coordinates, D = T_0 $ T_1 $ ... T_{k-1} $,
// and are converted to S coordinates (no terminators) on output.
//   esa  : SA (in D) at the last row of every BWT run            (toehold)
//   keys : SA values of rows x >= 1 that start a BWT run or hold a terminator in the BWT;
//   vals : for those rows, the index of SA[x-1] in esa (x-1 ends a run), or R + t for the few
//          terminator rows that do not start a run (their SA[x-1] is in `extra`);
//          then phi(i) = SA[x-1] + (i - keys[j]) with keys[j] the predecessor of i.  Sampling the terminator rows keeps phi from running across
//          string boundaries.
// ---------------------------------------------------------------------------
struct rindex {
    u64 n = 0;                 // BWT length = |S| + number of strings
    rlbwt bwt;
    int_vector<> esa;
    sparse_bv keys;
    int_vector<> vals;
    int_vector<> extra;
    sparse_bv dstart;          // string (dataset) starts in D
    const rlcsa_bwt *csa = nullptr;   // optional: backward search on an RLCSA (bwt.starts still names runs)

    rindex() {}
    rindex(const rindex &) = delete;
    rindex &operator=(const rindex &) = delete;

    inline u64 phi(u64 i) const {
        u64 t = keys.rank(i + 1);
        u64 j = keys.select(t);
        u64 v = vals[t - 1];
        return (v < bwt.R ? esa[v] : extra[v - bwt.R]) + (i - j);
    }
    inline u64 to_S(u64 d) const { return d - (dstart.rank(d + 1) - 1); }

    // leftmost/rightmost occurrence (S coordinates); T->bs = backward search, T->left = phi loop
    std::pair<u64, u64> query(const std::string &P, u64 &occ, timing *T = nullptr) {
        double t0 = T ? now() : 0;
        u64 m = P.size(), sp = 0, ep = n;
        u64 k = esa[bwt.R - 1];                 // SA[ep-1]
        occ = 0;
        for (u64 i = m; i-- > 0;) {
            unsigned char c = P[i];
            if (csa) {
                bool cov; u64 y;
                if (!csa->step_th(sp, ep, c, cov, y)) { if (T) T->bs += now() - t0; return {NONE, NONE}; }
                k = cov ? k - 1 : esa[bwt.run_of(y)] - 1;
                continue;
            }
            u64 a = bwt.rank(sp, c), b = bwt.rank(ep, c);
            if (a == b) { if (T) T->bs += now() - t0; return {NONE, NONE}; }
            if (bwt.at(ep - 1) == c) k = k - 1;
            else {
                u64 y = bwt.select(b - 1, c);   // last c in [sp, ep): end of a run
                k = esa[bwt.run_of(y)] - 1;
            }
            sp = bwt.C[c] + a; ep = bwt.C[c] + b;
        }
        double t1 = T ? now() : 0;
        u64 mn = k, mx = k;
        for (u64 x = ep - 1; x > sp; --x) {
            k = phi(k);
            if (k < mn) mn = k;
            if (k > mx) mx = k;
        }
        occ = ep - sp;
        if (T) { double t2 = now(); T->bs += t1 - t0; T->left += t2 - t1; }
        return {to_S(mn), to_S(mx)};
    }
    u64 serialize(std::ostream &out) const {
        out.write((char *)&n, 8);
        return 8 + bwt.serialize(out) + esa.serialize(out) + keys.serialize(out) + vals.serialize(out) + extra.serialize(out) + dstart.serialize(out);
    }
    void load(std::istream &in) {
        in.read((char *)&n, 8);
        bwt.load(in); esa.load(in); keys.load(in); vals.load(in); extra.load(in); dstart.load(in);
    }
};

} // namespace rz
