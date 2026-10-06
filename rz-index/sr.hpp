// Subsampled r-index (sr-index; Cobas, Gagie & Navarro) on a move structure, for leftmost /
// rightmost occurrences in a multi-string BWT.
//
// A phi entry (key, value) = (SA[x], SA[x-1]) exists for every BWT run start x >= 1 and for every
// terminator row x that is not a run start (as in rz::rindex).  Entries are subsampled by the
// sr rule over their values (run-end samples) in text order: going left to right, the previous
// entry is kept when the current value is more than s past the last kept value (the first and
// last entries, and entries tied to terminator rows, are always kept).  A run-end sample is kept
// iff its entry is kept.  phi(i) uses the kept predecessor key p <= i and is valid iff no original
// key lies in (p, i]; "trusted" keys are followed by a kept key, the others store the distance to
// the next original key.  When phi is invalid, the remaining range is LF-navigated run by run
// (at most s levels) until a kept run-end sample is found, as in the sr-index; the toehold walks
// at most s LF steps back to a kept sample.  With s = 0 every entry is kept: a plain r-index.
#pragma once
#include "rz_index.hpp"

namespace rz {

// RLFM-index backend with the same interface as move_bwt: rank/select on the RLBWT.
struct rl_nav {
    const rlbwt *b = nullptr;
    u64 n = 0, R = 0;
    u64 bytes_ = 0;
    explicit rl_nav(const rlbwt &bw) : b(&bw), n(bw.n), R(bw.R) {}
    struct interval { u64 sp, ep; };                      // rows [sp, ep)
    inline u64 start(u64 k) const { return k == R ? n : b->starts.select(k + 1); }
    inline u64 len(u64 k) const { return start(k + 1) - start(k); }
    inline unsigned char chr(u64 k) const { return b->heads[k]; }
    inline u64 at(u64 k, u64 o) const { return start(k) + o; }
    inline interval full() const { return {0, n}; }
    inline void LF(u64 &k, u64 &o) const {
        u64 row = start(k) + o;
        unsigned char c = chr(k);
        u64 r2 = b->C[c] + b->rank(row, c);
        k = b->run_of(r2); o = r2 - start(k);
    }
    inline bool step_toehold(interval &I, unsigned char c, u64 &tk, u64 &to, u64 &td) const {
        u64 a = b->rank(I.sp, c), e = b->rank(I.ep, c);
        if (a == e) return false;
        if (b->at(I.ep - 1) == c) ++td;
        else { u64 y = b->select(e - 1, c); tk = b->run_of(y); to = y - start(tk); td = 1; }
        I.sp = b->C[c] + a; I.ep = b->C[c] + e;
        return true;
    }
    inline u64 rows_sp(const interval &I) const { return I.sp; }
    inline u64 rows_ep(const interval &I) const { return I.ep; }
    inline void first(const interval &I, u64 &k, u64 &o) const { k = b->run_of(I.sp); o = I.sp - start(k); }
    inline void init_toehold(u64 &tk, u64 &to) const { tk = R - 1; to = len(R - 1) - 1; }
};

// RLCSA backend: backward search (with toehold) and the rank in LF steps on the RLCSA; the run
// starts and heads of the RLBWT are still used to name runs and read BWT characters.
struct csa_nav : rl_nav {
    const rlcsa_bwt *x = nullptr;
    csa_nav(const rlbwt &bw, const rlcsa_bwt &cs) : rl_nav(bw), x(&cs) {}
    inline void LF(u64 &k, u64 &o) const {
        u64 row = start(k) + o;
        unsigned char c = chr(k);
        u64 r2 = x->C[c] + x->rank(c, row);
        k = b->run_of(r2); o = r2 - start(k);
    }
    inline bool step_toehold(interval &I, unsigned char c, u64 &tk, u64 &to, u64 &td) const {
        bool cov; u64 y;
        if (!x->step_th(I.sp, I.ep, c, cov, y)) return false;
        if (cov) ++td;
        else { tk = b->run_of(y); to = y - start(tk); td = 1; }
        return true;
    }
};

template <class NAV>
struct srindex_t {
    u64 n = 0, s = 0, R = 0;
    sparse_bv keys;                 // kept keys
    int_vector<> vals;              // their values
    bit_vector trusted;  rank_support_v<0> trusted_r0;
    int_vector<> area;              // for untrusted keys: next original key - key
    bit_vector ekept;    rank_support_v5<1> ekept_r;
    int_vector<> esav;              // kept run-end samples, in run order
    std::vector<u64> trow, tpos;    // terminator rows and their SA values (sorted by row)
    sparse_bv dstart;               // dataset starts in D
    NAV *mv = nullptr;              // backward-search backend, not owned
    static inline thread_local u64 lf_steps = 0, invalid_phis = 0;   // per thread (rz-classify -j)

    srindex_t() {}
    srindex_t(const srindex_t &) = delete;
    srindex_t &operator=(const srindex_t &) = delete;

    inline u64 to_S(u64 d) const { return d - (dstart.rank(d + 1) - 1); }
    inline std::pair<u64, bool> phi(u64 i) const {
        u64 t = keys.rank(i + 1);
        if (t == 0) return {0, false};          // no kept key <= i (only after the last row of a walk)
        u64 p = keys.select(t);
        u64 d = i - p;
        bool ok = trusted[t - 1] || d < area[trusted_r0(t - 1)];
        return {vals[t - 1] + d, ok};
    }
    // SA of row (k, o) if it is sampled
    inline bool sample(u64 k, u64 o, u64 &v) const {
        if (mv->chr(k) == 1) {
            u64 row = mv->at(k, o);
            u64 j = std::lower_bound(trow.begin(), trow.end(), row) - trow.begin();
            v = tpos[j];
            return true;
        }
        if (o + 1 == mv->len(k) && ekept[k]) { v = esav[ekept_r(k)]; return true; }
        return false;
    }

    // ---- locate over a range, reporting min and max -------------------------------------------
    struct run { u64 k, lo, hi; };               // rows lo..hi (offsets) of move row k; empty if hi < lo
    static inline thread_local u64 value = 0; static inline thread_local bool valid = false; static inline thread_local u64 mn = 0, mx = 0;
    static inline thread_local std::vector<u64> *collect = nullptr;   // if set, every occurrence (in D coordinates) is appended
    inline void report(u64 v) { if (v < mn) mn = v; if (v > mx) mx = v; if (collect) collect->push_back(v); }

    // runs covering the rows [(k, o), (k, o) + cnt)
    void split(u64 k, u64 o, u64 cnt, std::vector<run> &out) const {
        out.clear();
        while (cnt) {
            u64 l = mv->len(k) - o, t = std::min(l, cnt);
            out.push_back({k, o, o + t - 1});
            cnt -= t; ++k; o = 0;
        }
    }
    void runs_at(std::vector<run> &rs, u64 level) {
        for (u64 j = rs.size(); j-- > 0;) {
            run r = rs[j];
            if (r.hi + 1 == r.lo) continue;
            long long hi = (long long)r.hi, lo = (long long)r.lo;
            do {
                while (hi >= lo && valid) {
                    report(value); --hi;
                    auto pv = phi(value); value = pv.first; valid = pv.second;
                    if (!valid) ++invalid_phis;
                }
                if (hi < lo) break;
                u64 v;
                if (sample(r.k, (u64)hi, v)) { value = v + level; valid = true; }
            } while (valid);
            if (hi >= lo) navigate(r.k, (u64)lo, (u64)(hi - lo + 1), level);
        }
    }
    void navigate(u64 k, u64 o, u64 cnt, u64 level) {
        if (mv->chr(k) == 1) { std::cerr << "sr: navigating a terminator run\n"; exit(1); }
        mv->LF(k, o); ++lf_steps;
        if (level + 1 >= s) {                       // phi is valid again
            for (u64 i = 0; i < cnt; ++i) { report(value); value = phi(value).first; }
            valid = true;
            return;
        }
        std::vector<run> rs;
        split(k, o, cnt, rs);
        runs_at(rs, level + 1);
    }

    // leftmost / rightmost occurrence in S coordinates
    std::pair<u64, u64> query(const std::string &P, u64 &occ, timing *T = nullptr) {
        double t0 = T ? now() : 0;
        u64 m = P.size();
        typename NAV::interval I = mv->full();
        u64 tk, to, td = 0;
        mv->init_toehold(tk, to);
        occ = 0;
        for (u64 i = m; i-- > 0;)
            if (!mv->step_toehold(I, (unsigned char)P[i], tk, to, td)) { if (T) T->bs += now() - t0; return {NONE, NONE}; }
        double t1 = T ? now() : 0;
        std::pair<u64, u64> r = ends(I, tk, to, td, occ);
        if (T) { double t2 = now(); T->bs += t1 - t0; T->left += t2 - t1; }
        return r;
    }
    // leftmost and rightmost occurrences, given the interval I of a pattern found by step_toehold
    // (with its toehold tk, to, td)
    std::pair<u64, u64> ends(const typename NAV::interval &I, u64 tk, u64 to, u64 td, u64 &occ) {
        u64 v, j = 0;
        while (!sample(tk, to, v)) { mv->LF(tk, to); ++j; ++lf_steps; }
        u64 last = v + j - td;                       // SA[ep - 1]
        u64 sp = mv->rows_sp(I), ep = mv->rows_ep(I);
        occ = ep - sp;
        mn = mx = last; if (collect) collect->push_back(last);
        if (occ > 1) {
            auto pv = phi(last); value = pv.first; valid = pv.second;
            if (!valid) ++invalid_phis;
            std::vector<run> rs;
            // rows sp .. ep-2
            u64 k, o;
            mv->first(I, k, o);
            split(k, o, occ - 1, rs);
            runs_at(rs, 0);
        }
        return {to_S(mn), to_S(mx)};
    }

    u64 serialize(std::ostream &out) const {
        out.write((char *)&n, 8); out.write((char *)&s, 8); out.write((char *)&R, 8);
        u64 w = 24 + keys.serialize(out) + vals.serialize(out) + trusted.serialize(out) + trusted_r0.serialize(out)
                + area.serialize(out) + ekept.serialize(out) + ekept_r.serialize(out) + esav.serialize(out);
        u64 K = trow.size(); out.write((char *)&K, 8);
        out.write((char *)trow.data(), K * 8); out.write((char *)tpos.data(), K * 8);
        return w + 8 + 16 * K + dstart.serialize(out);
    }
    void load(std::istream &in) {
        in.read((char *)&n, 8); in.read((char *)&s, 8); in.read((char *)&R, 8);
        keys.load(in); vals.load(in); trusted.load(in); trusted_r0.load(in, &trusted);
        area.load(in); ekept.load(in); ekept_r.load(in, &ekept); esav.load(in);
        u64 K; in.read((char *)&K, 8); trow.resize(K); tpos.resize(K);
        in.read((char *)trow.data(), K * 8); in.read((char *)tpos.data(), K * 8);
        dstart.load(in);
    }
};

typedef srindex_t<move_bwt> srindex;
typedef srindex_t<rl_nav> srindex_rl;
typedef srindex_t<csa_nav> srindex_csa;

}  // namespace rz
