// Read classification with the rz-index or the sr-index.  For each read, the features are either
//  (default) the phrases of the Ziv-Merhav cross parse of the read against the collection, computed
//      right to left as in Cliffy (Ahmed, Boucher, Langmead, Genome Res 2025): backward search from the
//      end of the read until the interval is empty, then restart where the previous phrase started; or
//  (-M) the MEMs (exact matches maximal in both directions) of the read, found left to right by
//      forward-backward: backward search from one past the end of the previous MEM to its leftmost
//      start s, then forward search from s to the MEM's end.  A forward step on P is a backward step on
//      rc(P), which is in the same BWT because S has both strands (and digests are strand-symmetric).
// For each feature, the leftmost and rightmost documents (genomes) containing it or its reverse
// complement are reported, with the rz-index's grids or (-S) with the sr-index.
//
// usage: rz-classify [options] <index.rz> <reads.fq> <out.listings>
//   out.listings has Cliffy's listing format, so the Cliffy scoring scripts can read it:
//     >readname
//     [start,end] {leftdoc,rightdoc} [start,end] {leftdoc,rightdoc} ...
//   (read positions 0-based and inclusive; documents 0-based, in the order of S)
//   -M : MEMs instead of Ziv-Merhav phrases
//   -L bases : Boyer-Moore-Li: only MEMs of at least that many bases (with -B: digested MEMs whose span in
//          the read is at least that many bases); positions in the listings are read positions (bases)
//   -l : with -L and -S, list the documents of all occurrences of each MEM instead of the two ends
//   -H T, -A a : with -l, answer a MEM by its LCA (rz-index grids, loaded with the index) instead of listing
//          it when occ > T, or when occ >= a * (sequences in the documents between its ends: -N file with the
//          number of sequences per document, one per line); LCA answers are written <l,r>
//   -T f : with -L (and no -S), answer from the tag array f (rz-tagbuild): the documents of each MEM (-l) or
//          its LCA (the smallest and largest document); the search uses the RLBWT (or the RLCSA, with -C)
//   -V f : with -l on a digest, verify each occurrence against the DNA (rz-vfybuild file f): align at the
//          MEM's middle minimizer and extend the exact match; report the maximal verified matches of at
//          least L bases, each with the genera whose occurrences contain it
//   -U : with -L, homopolymer compression (an index built with rz-prep -U): each read's runs of equal bases
//          are collapsed before BML, L counts compressed bases, and the listings give each MEM's span in the
//          original read (bases), so its weight is its length in the read
//   -m minlen : ignore features shorter than minlen (default 1: report all, as Cliffy does)
//   -x f : the rz-index's aux structures
//   -C f : do the backward searches on the RLCSA f (rz-csabuild) instead of the RLBWT
//   -B f : the index is over a byte-mapped minimizer digest (rz-prep -B, map file f): digest each read
//          the same way first; read positions in the listings are then digest positions
//   -S f -R g : use the sr-index f (sr-build) on the RLBWT stored in the r-index g (or on the RLCSA, with
//          -C) for the searches and the ends, instead of the rz-index (which is then used only to map
//          positions to documents)
// A separator X in a digested read (from a non-ACGT character, e.g. N) is never matched: features stop at it,
// as undigested features stop at N (X also separates A/C/G/T runs inside reference genomes).
// All searches are timed: the reported time per read includes finding the features.
// rz-index:  ZM: one backward search per phrase (+ the failing step), then the rc search and the grids.
//            MEM: backward scan, forward (rc) scan, and a second backward scan from the MEM's true end
//            when it differs from where the first scan started; then the grids.
// sr-index:  ZM: one backward search with toehold per phrase.  MEM: backward scan, then forward (rc)
//            scan with toehold; rc(MEM) has the same documents as the MEM.
#include "sr.hpp"
#include "vfy.hpp"
#include "tag.hpp"
#include "kbloom.hpp"
#include <unordered_map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <fstream>
#include <unistd.h>
using std::string; typedef uint64_t u64;

struct RzEng {                                   // backward search on the rz-index's backend
    static const bool isrz = true;
    rz::index &Z;
    struct St { u64 sp, ep; };
    St full() const { return {0, Z.bwt.n}; }
    static u64 sp(const St &s) { return s.sp; }
    static u64 ep(const St &s) { return s.ep; }
    bool step(St &s, unsigned char c) const {
        u64 a = s.sp, b = s.ep; Z.ext(a, b, c);
        if (a >= b) return false;
        s.sp = a; s.ep = b; return true;
    }
};
template <class NAV> struct SrEng {              // backward search with toehold on the sr-index's backend
    static const bool isrz = false;
    rz::srindex_t<NAV> &S;
    struct St { typename NAV::interval I; u64 tk, to, td; };
    St full() const { St s; s.I = S.mv->full(); S.mv->init_toehold(s.tk, s.to); s.td = 0; return s; }
    static u64 sp(const St &s) { return s.I.sp; }
    static u64 ep(const St &s) { return s.I.ep; }
    bool step(St &s, unsigned char c) const {
        St t = s;
        if (!S.mv->step_toehold(t.I, c, t.tk, t.to, t.td)) return false;
        s = t; return true;
    }
    std::pair<u64, u64> ends(const St &s) { u64 occ; return S.ends(s.I, s.tk, s.to, s.td, occ); }
    static inline thread_local std::vector<u64> occs;   // per thread (-j)
    // all occurrences, as positions in S
    void list(const St &s, std::vector<u64> &out) {
        occs.clear(); S.collect = &occs; u64 occ; S.ends(s.I, s.tk, s.to, s.td, occ); S.collect = nullptr;
        out.clear(); for (u64 v : occs) out.push_back(S.to_S(v));
    }
};


// RZ_COUNTS=1: list each genus with the number of occurrences of the MEM in it, as doc:count
static bool COUNTS = false;
static void append_docs(const std::vector<u64> &pos, rz::index &Z, std::vector<u64> &docs, string &out, u64 &nlisted) {
    char b[48]; docs.clear();
    for (u64 v : pos) docs.push_back(Z.genome_of(v));
    std::sort(docs.begin(), docs.end());
    bool first = true;
    for (u64 q = 0; q < docs.size();) {
        u64 r = q; while (r < docs.size() && docs[r] == docs[q]) ++r;
        if (COUNTS) snprintf(b, sizeof b, first ? "%lu:%lu" : ",%lu:%lu", docs[q], r - q);
        else snprintf(b, sizeof b, first ? "%lu" : ",%lu", docs[q]);
        out += b; first = false; ++nlisted; q = r;
    }
}

struct Stats { u64 nreads = 0, nfeat = 0, totlen = 0, nfail = 0, steps = 0, nuncl = 0, nlisted = 0, nocc = 0, nver = 0; double tsearch = 0, tquery = 0, tver = 0; };
static const rz::vfy_index *VFY = nullptr;
static rz::tag_index *TAGX = nullptr;         // -T: answer from the tag array (list, or LCA = smallest and largest tag)
static bool TAGLIST = false;
static thread_local std::vector<u64> TL;   // -V: verify digested MEMs against the DNA
// answer a MEM from the tag array: its genera (-l), with occurrence counts (RZ_COUNTS), or the smallest and largest
static thread_local std::vector<std::pair<u64, u64>> TC;
// RZ_IVLOG=file: log the BWT interval [sp, ep) of every MEM answered from the tag array (two u64 each), for benchmarks
static FILE *IVLOG = nullptr; static std::mutex IVLOG_M;
static void tag_answer(u64 sp, u64 ep, string &out, u64 &nlisted) {
    if (IVLOG) { std::lock_guard<std::mutex> g(IVLOG_M); u64 v[2] = {sp, ep}; fwrite(v, 8, 2, IVLOG); }
    char b[48];
    if (TAGLIST && COUNTS) {
        TAGX->count(sp, ep, TC);
        for (u64 q = 0; q < TC.size(); ++q) { snprintf(b, sizeof b, q ? ",%lu:%lu" : "%lu:%lu", TC[q].first, TC[q].second); out += b; }
        nlisted += TC.size(); return;
    }
    TL.clear(); TAGX->list(sp, ep, TL); std::sort(TL.begin(), TL.end());
    if (TAGLIST) { for (u64 q = 0; q < TL.size(); ++q) { snprintf(b, sizeof b, q ? ",%lu" : "%lu", TL[q]); out += b; } nlisted += TL.size(); }
    else { snprintf(b, sizeof b, "%lu,%lu", TL.front(), TL.back()); out += b; }
}

// -F t: table of the search states (BWT interval, plus the toehold for the sr-index) of all strings of t
// symbols (A/C/G/T, or digest bytes with -B), as in Bowtie 2's and Cliffy's ftab, so that a backward search
// from scratch starts with one lookup instead of t steps.  Not used where every suffix's interval is needed
// (the rz-index's grids, -H/-A).
template <class St> struct GTab {
    int t = 0, sig = 4; std::vector<St> v; std::vector<uint8_t> ok; static inline thread_local u64 lookups = 0;
    template <class E> void build(E &E_, int tt, int sg) {
        t = tt; sig = sg;
        std::vector<St> a{E_.full()}; std::vector<uint8_t> f{1};
        for (int i = 0; i < t; ++i) {                   // strings of length i -> i + 1, key = c * sig^i + key(w)
            u64 m = a.size(); std::vector<St> a2(sig * m); std::vector<uint8_t> f2(sig * m, 0);
            for (int c = 0; c < sig; ++c) {
                unsigned char ch = sig == 4 ? (unsigned char)"ACGT"[c] : (unsigned char)c;
                for (u64 w = 0; w < m; ++w) if (f[w]) { St I = a[w]; if (E_.step(I, ch)) { a2[c * m + w] = I; f2[c * m + w] = 1; } }
            }
            a.swap(a2); f.swap(f2);
        }
        v.swap(a); ok.swap(f);
    }
    u64 bytes() const { return v.size() * (sizeof(St) + 1); }
    // state for the t symbols p[0..t-1], or false if they do not occur
    bool get(const unsigned char *p, St &I) {
        u64 key = 0, mul = 1;
        for (int i = t - 1; i >= 0; --i) { key += (u64)(sig == 4 ? rz::dg_b2((char)p[i]) : p[i]) * mul; mul *= sig; }
        ++lookups;
        if (!ok[key]) return false;
        I = v[key]; return true;
    }
};
template <class St> GTab<St> &gtab() { static GTab<St> g; return g; }
static int FT_t = 0; static u64 FT_bytes = 0;
// -H T / -A a (with -l): answer a MEM by its LCA (the rz-index's grids) instead of listing its occurrences
// when it has more than T occurrences, or when occ >= a * (number of sequences in the documents between its two
// ends, i.e. in the LCA's subtree; -N gives the sequences per document)
static u64 TRIM_FIX = 0; static bool TRIM_FIXON = false;
static u64 TRIM_K = 0, TRIM_W = 0; static thread_local u64 TR_mems = 0, TR_trim = 0, TR_skip = 0;   // RZ_TRIM=k,w: emulate a phrase index (minimizer phrases) by trimming MEMs
static u64 KEBAB_K = 0; static thread_local u64 KB_total = 0, KB_kept = 0;   // RZ_KEBAB=k: ideal KeBaB pseudo-MEMs before BML
static rz::kbloom KBF; static bool KBON = false;   // -K file.kbf: KeBaB with a blocked Bloom filter of the text's canonical k-mers (rz-kbbuild)
static u64 HYB_T = 0; static double HYB_A = 0; static bool HYB = false; static thread_local u64 NHYB = 0;
static std::vector<u64> NSEQ;               // -N: prefix sums of the number of sequences per document
template <class E, class SV>
std::pair<u64, u64> grid_ends(rz::index &Z, const SV &suf, const SV &xiv, u64 m) {
    Z.fs.resize(m + 1); Z.fe.resize(m + 1); Z.xs.resize(m + 1); Z.xe.resize(m + 1);
    for (u64 q = 0; q <= m; ++q) { Z.fs[q] = E::sp(suf[m - q]); Z.fe[q] = E::ep(suf[m - q]); Z.xs[q] = E::sp(xiv[q]); Z.xe[q] = E::ep(xiv[q]); }
    return Z.query_with_both(m);
}

template <class E>
void classify_read(E &E_, rz::index &Z, const string &s, bool mem, u64 minlen, string &line, Stats &st) {
    typedef typename E::St St;
    static thread_local std::vector<St> ivs, iv2, xiv;
    char buf[96];
    auto emit = [&](u64 i, u64 e, std::pair<u64, u64> r) {   // feature s[i..e)
        if (r.first == rz::NONE) { ++st.nfail; fprintf(stderr, "FAIL read %lu feature [%lu,%lu):", st.nreads, i, e); for (u64 k = i; k < e; ++k) fprintf(stderr, " %d", (unsigned char)s[k]); fprintf(stderr, "\n"); return; }
        snprintf(buf, sizeof buf, "[%lu,%lu] {%lu,%lu} ", i, e - 1, Z.genome_of(r.first), Z.genome_of(r.second));
        line += buf; ++st.nfeat; st.totlen += e - i;
    };
    if (!mem) {                                  // Ziv-Merhav cross parse, right to left
        u64 e = s.size();
        while (e > 0) {
            double t0 = rz::now();
            St c = E_.full(); u64 i = e;
            if (E::isrz) { ivs.clear(); ivs.push_back(c); }
            while (i > 0) { ++st.steps; if (s[i - 1] == 'X' || !E_.step(c, (unsigned char)s[i - 1])) break; --i; if (E::isrz) ivs.push_back(c); }
            double t1 = rz::now(); st.tsearch += t1 - t0;
            if (i == e) { --e; continue; }       // character absent from the index
            if (e - i >= minlen) {
                std::pair<u64, u64> r;
                if constexpr (E::isrz) {
                    u64 m = e - i; Z.fs.resize(m + 1); Z.fe.resize(m + 1);
                    for (u64 j = 0; j <= m; ++j) { Z.fs[j] = ivs[m - j].sp; Z.fe[j] = ivs[m - j].ep; }
                    r = Z.query_with_suffixes(s.substr(i, m));
                } else r = E_.ends(c);
                st.tquery += rz::now() - t1;
                emit(i, e, r);
            }
            e = i;
        }
        return;
    }
    u64 n = s.size(), p = 0;                     // MEMs, left to right; p = one past the previous MEM's end
    while (p < n) {
        double t0 = rz::now();
        St c = E_.full(); u64 i = p + 1;         // backward from s[p]: matches s[i..p]
        if (E::isrz) { ivs.clear(); ivs.push_back(c); }
        while (i > 0) { ++st.steps; if (s[i - 1] == 'X' || !E_.step(c, (unsigned char)s[i - 1])) break; --i; if (E::isrz) ivs.push_back(c); }
        if (i == p + 1) { ++p; st.tsearch += rz::now() - t0; continue; }   // s[p] absent from the index
        u64 b = i, e = b;                        // forward from b: rc(s[b..e))
        St x = E_.full();
        if (E::isrz) { xiv.clear(); xiv.push_back(x); }
        while (e < n) { ++st.steps; if (s[e] == 'X' || !E_.step(x, rz::comp((unsigned char)s[e]))) break; ++e; if (E::isrz) xiv.push_back(x); }
        // the MEM is s[b..e), e >= p + 1
        u64 m = e - b;
        if constexpr (E::isrz) {
            const std::vector<St> *suf = &ivs;   // suffix intervals of s[b..e)
            if (e != p + 1) {                    // second backward scan, from the MEM's true end
                St y = E_.full(); iv2.clear(); iv2.push_back(y);
                for (u64 k = e; k > b; --k) { ++st.steps; E_.step(y, (unsigned char)s[k - 1]); iv2.push_back(y); }
                suf = &iv2;
            }
            double t1 = rz::now(); st.tsearch += t1 - t0;
            if (m >= minlen) {
                Z.fs.resize(m + 1); Z.fe.resize(m + 1); Z.xs.resize(m + 1); Z.xe.resize(m + 1);
                for (u64 j = 0; j <= m; ++j) { Z.fs[j] = (*suf)[m - j].sp; Z.fe[j] = (*suf)[m - j].ep; Z.xs[j] = xiv[j].sp; Z.xe[j] = xiv[j].ep; }
                auto r = Z.query_with_both(m);
                st.tquery += rz::now() - t1;
                emit(b, e, r);
            }
        } else {
            double t1 = rz::now(); st.tsearch += t1 - t0;
            if (m >= minlen) { auto r = E_.ends(x); st.tquery += rz::now() - t1; emit(b, e, r); }
        }
        p = e;
    }
}

// ---- Boyer-Moore-Li on the digest, with the length threshold L in bases ---------------------------
// The read (DNA) is split into maximal A/C/G/T runs; each run is digested exactly as the reference
// (byte-mapped minimizers, all tied minima, equal consecutive symbols collapsed).  Window y covers
// bases [y, y+w) and selects the digest symbols lo[y]..hi[y] (both non-decreasing in y).  A DNA match
// R[x..x+L) contains exactly the windows y in [x, x+L-w], so its digest is dg[lo[x]..hi[x+L-w]]; if that
// does not occur in the reference digest, neither does R[x..x+L) (no false negatives).
//   candidate x: backward-search dg[lo[x]..hi[x+L-w]] from its right end.
//   failure at symbol j (dg[j..hi] does not occur): every x' >= x with lo[x'] <= j fails too, so jump to
//     the first y with lo[y] > j.
//   success: extend left (backward steps) and right (backward steps on the reverse complement) to a
//     maximal digested MEM D = dg[a..b]; its span in the read is the largest DNA interval whose windows
//     all select symbols in [a, b]: from the first window with lo >= a to just before the end of the first
//     window with hi > b.  Report D with that span; continue from the first x whose last window has hi > b.
// This finds every digested MEM whose span is at least L bases.  For each, the rz-index (or sr-index)
// gives the leftmost and rightmost documents, or (-l) the sr-index lists the documents of all occurrences.
// Hybrid answer (-H/-A): decide whether to answer the MEM P[aa..bb] (inclusive; the backward scan started
// at b, inclusive) by its two ends (rz-index grids) rather than by listing.  Returns true and the ends in r
// if so.  ivs holds the suffix intervals from the first backward scan, xiv the reverse-complement prefixes.
template <class E, class SV, class CH>
bool hyb_lca(E &E_, rz::index &Z, const typename E::St &xr, SV &ivs, SV &iv2, SV &xiv, u64 m, u64 b, u64 bb, u64 aa, CH ch, Stats &st, std::pair<u64, u64> &r) {
    u64 occ = E::ep(xr) - E::sp(xr);
    bool big = HYB_T && occ > HYB_T;
    if (!big && HYB_A <= 0) return false;
    const SV *suf = &ivs;
    if (bb != b) {                                   // suffix intervals of the MEM itself
        typename E::St y = E_.full(); iv2.clear(); iv2.push_back(y);
        for (u64 q = bb + 1; q > aa; --q) { ++st.steps; E_.step(y, ch(q - 1)); iv2.push_back(y); }
        suf = &iv2;
    }
    r = grid_ends<E>(Z, *suf, xiv, m);
    if (r.first == rz::NONE) return false;
    if (!big) { u64 gl = Z.genome_of(r.first), gr = Z.genome_of(r.second); if ((double)occ < HYB_A * (NSEQ[gr + 1] - NSEQ[gl])) return false; }
    ++NHYB; return true;
}

template <class E>
void classify_bml(E &E_, rz::index &Z, const string &read, u64 L, bool list, string &line, Stats &st) {
    typedef typename E::St St;
    const rz::dg_bytemap &M = rz::bytemap();
    const u64 k = M.k, w = M.w, Wk = w - k + 1;            // window: w bases = Wk k-mers
    static thread_local std::vector<St> ivs, iv2, xiv;
    static thread_local std::vector<u64> h, mlo, mhi, idx, lo, hi, pos, docs;
    static thread_local std::vector<uint32_t> f;
    static thread_local string dg;
    char buf[96];
    u64 nfeat0 = st.nfeat;
    u64 i = 0, n0 = read.size();
    while (i < n0) {
        while (i < n0 && rz::dg_b2(read[i]) < 0) ++i;
        u64 r0 = i;
        while (i < n0 && rz::dg_b2(read[i]) >= 0) ++i;
        u64 n = i - r0;
        if (n < L || n < w) continue;
        double t0 = rz::now();
        // k-mer hashes and codes
        u64 nk = n - k + 1, nwin = nk - Wk + 1;
        h.resize(nk); f.resize(nk);
        uint32_t mask = (uint32_t)((1ull << (2 * k)) - 1), fw = 0;
        for (u64 t = 0; t < n; ++t) {
            fw = ((fw << 2) | (uint32_t)rz::dg_b2(read[r0 + t])) & mask;
            if (t + 1 >= k) { u64 p = t + 1 - k; uint32_t rc = rz::dg_rc_code(fw, k); f[p] = fw; h[p] = rz::dg_mix(fw < rc ? fw : rc); }
        }
        // windows: first and last tied minimum
        mlo.resize(nwin); mhi.resize(nwin);
        std::vector<char> sel(nk, 0); std::deque<u64> dq;
        for (u64 p = 0; p < nk; ++p) {
            while (!dq.empty() && h[dq.back()] > h[p]) dq.pop_back();
            dq.push_back(p);
            if (p + 1 >= Wk) {
                u64 y = p + 1 - Wk;
                while (dq.front() < y) dq.pop_front();
                u64 mn = h[dq.front()], last = dq.front();
                for (u64 q : dq) { if (h[q] != mn) break; sel[q] = 1; last = q; }
                mlo[y] = dq.front(); mhi[y] = last;
            }
        }
        dg.clear(); idx.assign(nk, 0);
        for (u64 p = 0; p < nk; ++p) if (sel[p]) {
            char c = (char)M.code2byte[f[p]];
            if (dg.empty() || dg.back() != c) dg.push_back(c);
            idx[p] = dg.size() - 1;
        }
        lo.resize(nwin); hi.resize(nwin);
        for (u64 y = 0; y < nwin; ++y) { lo[y] = idx[mlo[y]]; hi[y] = idx[mhi[y]]; }
        static thread_local std::vector<u64> rfirst, rlast;
        if (VFY) { rfirst.assign(dg.size(), ~0ull); rlast.assign(dg.size(), 0);
            for (u64 p = 0; p < nk; ++p) if (sel[p]) { u64 q = idx[p]; if (rfirst[q] == ~0ull) rfirst[q] = p; rlast[q] = p; } }
        struct vint { u64 s, e, d; };
        std::vector<vint> vi;                        // verified matches (read-run coordinates, document)
        st.tsearch += rz::now() - t0;
        u64 x = 0, nd = dg.size();
        static const bool check = getenv("RZ_BML_CHECK") != nullptr;
        std::vector<std::pair<u64, u64>> found;
        struct feat { u64 l, r; string d; };
        std::vector<feat> fr; fr.clear();
        while (x + L <= n) {
            double t1 = rz::now();
            u64 a = lo[x], b = hi[x + L - w];
            St c = E_.full(); u64 j = b + 1;
            if (E::isrz || HYB) { ivs.clear(); ivs.push_back(c); }
            auto &G = gtab<St>(); bool useG = G.t && !HYB && !VFY && (!E::isrz || TAGX);
            if (useG && b + 1 >= a + G.t) {
                St c2 = c; if (G.get((const unsigned char *)dg.data() + b + 1 - G.t, c2)) { c = c2; j = b + 1 - G.t; ++st.steps; }
            }
            bool fail = false;
            while (j > a) { ++st.steps; if (!E_.step(c, (unsigned char)dg[j - 1])) { fail = true; break; } --j; if (E::isrz || HYB) ivs.push_back(c); }
            if (fail) {                              // dg[j-1..b] does not occur
                u64 nx = std::upper_bound(lo.begin(), lo.end(), j - 1) - lo.begin();
                st.tsearch += rz::now() - t1;
                if (nx >= nwin) break;
                x = nx; continue;
            }
            while (j > 0) { ++st.steps; if (!E_.step(c, (unsigned char)dg[j - 1])) break; --j; if (E::isrz || HYB) ivs.push_back(c); }
            u64 aa = j, e = aa;                      // forward (rc) scan from aa
            St xr = E_.full();
            if (E::isrz || HYB) { xiv.clear(); xiv.push_back(xr); }
            if (useG && e + G.t <= nd) {
                unsigned char rcb[8]; for (int q = 0; q < G.t; ++q) rcb[q] = rz::comp((unsigned char)dg[e + G.t - 1 - q]);
                St x2 = xr; if (G.get(rcb, x2)) { xr = x2; e += G.t; ++st.steps; }
            }
            while (e < nd) { ++st.steps; if (!E_.step(xr, rz::comp((unsigned char)dg[e]))) break; ++e; if (E::isrz || HYB) xiv.push_back(xr); }
            u64 bb = e - 1, m = e - aa;              // D = dg[aa..bb]
            u64 left = std::lower_bound(lo.begin(), lo.end(), aa) - lo.begin();
            u64 yR = std::upper_bound(hi.begin(), hi.end(), bb) - hi.begin();
            u64 right = yR >= nwin ? n - 1 : yR + w - 2;
            std::pair<u64, u64> r;
            double t2;
            if constexpr (E::isrz) {
              if (TAGX) {                            // answer from the tag array (rc(D) has the same documents)
                t2 = rz::now(); st.tsearch += t2 - t1;
                string d = "{"; { u64 nl = 0; tag_answer(E::sp(xr), E::ep(xr), d, nl); st.nlisted += nl; }
                d += "} "; fr.push_back({left, right, d});
              } else {
                const std::vector<St> *suf = &ivs;
                if (bb != b) {
                    St y = E_.full(); iv2.clear(); iv2.push_back(y);
                    for (u64 q = e; q > aa; --q) { ++st.steps; E_.step(y, (unsigned char)dg[q - 1]); iv2.push_back(y); }
                    suf = &iv2;
                }
                t2 = rz::now(); st.tsearch += t2 - t1;
                Z.fs.resize(m + 1); Z.fe.resize(m + 1); Z.xs.resize(m + 1); Z.xe.resize(m + 1);
                for (u64 q = 0; q <= m; ++q) { Z.fs[q] = (*suf)[m - q].sp; Z.fe[q] = (*suf)[m - q].ep; Z.xs[q] = xiv[q].sp; Z.xe[q] = xiv[q].ep; }
                r = Z.query_with_both(m);
                if (r.first == rz::NONE) ++st.nfail;
                else { snprintf(buf, sizeof buf, "{%lu,%lu} ", Z.genome_of(r.first), Z.genome_of(r.second)); fr.push_back({left, right, buf}); }
              }
            } else {
                t2 = rz::now(); st.tsearch += t2 - t1;
                if (!list) {
                    r = E_.ends(xr);
                    snprintf(buf, sizeof buf, "{%lu,%lu} ", Z.genome_of(r.first), Z.genome_of(r.second)); fr.push_back({left, right, buf});
                } else if (HYB && hyb_lca<E>(E_, Z, xr, ivs, iv2, xiv, m, b, bb, aa, [&](u64 q) { return (unsigned char)dg[q]; }, st, r)) {
                    snprintf(buf, sizeof buf, "<%lu,%lu> ", Z.genome_of(r.first), Z.genome_of(r.second)); fr.push_back({left, right, buf});
                } else {
                    E_.list(xr, pos); docs.clear();
                    if (VFY) {                       // check every occurrence against the DNA
                        double tv = rz::now();
                        static thread_local std::unordered_map<uint32_t, string> cache; static thread_local u64 cread = ~0ull;
                        if (cread != st.nreads) { cache.clear(); cread = st.nreads; }
                        const rz::vfy_index &V = *VFY; const u64 kk = rz::bytemap().k, t = m / 2;
                        const char *RR = read.data() + r0;
                        for (u64 v : pos) {
                            ++st.nocc;
                            u64 tt = Z.Bs.rank(v + 1) - 1, s0 = Z.Bs.select(tt + 1), doc = tt / 2, o = v - s0;
                            // the listed occurrences are of rc(D) (the state of the reverse-complement scan): an
                            // occurrence of rc(D) in one strand is an occurrence of D in the other
                            bool rev = !(tt & 1);
                            u64 of = rev ? o : V.dlg[doc] - o - m;     // forward symbols [of, of+m) hold D (or rc(D) if rev)
                            u64 j = std::upper_bound(V.dstart.begin() + V.dfirst[doc], V.dstart.begin() + V.dfirst[doc + 1], (uint32_t)of) - V.dstart.begin() - 1;
                            u64 qf = of - V.dstart[j];
                            if (qf + m > V.dlen[j]) continue;          // should not happen
                            auto it = cache.find(j);
                            if (it == cache.end()) { it = cache.emplace(j, string()).first; V.run(j, it->second); }
                            const string &F = it->second; u64 rn = F.size();
                            u64 rp, fp;                                    // read and reference anchor positions
                            if (!rev) { rp = rfirst[aa + t]; fp = V.sym_pos(j, qf + t); }
                            else { rp = rlast[aa + t]; fp = rn - V.sym_pos(j, qf + m - 1 - t) - kk; }
                            auto dnacomp = [](char c) -> char { switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; } };
                            auto ref = [&](u64 x) -> char { return rev ? dnacomp(F[rn - 1 - x]) : F[x]; };
                            static thread_local int dbg = getenv("RZ_VDBG") ? 20 : 0;
                            if (dbg > 0) { --dbg; string a(RR + rp, 12), b; for (u64 x = 0; x < 12 && fp + x < rn; ++x) b.push_back(ref(fp + x));
                                fprintf(stderr, "VDBG doc %lu rev %d o %lu Lg %u of %lu run %lu qf %lu dlen %u m %lu t %lu rp %lu fp %lu rn %lu read %s ref %s\n",
                                        doc, (int)rev, o, V.dlg[doc], of, j, qf, V.dlen[j], m, t, rp, fp, rn, a.c_str(), b.c_str()); }
                            u64 r = 0, l = 0;
                            while (rp + r < n && fp + r < rn && RR[rp + r] == ref(fp + r)) ++r;
                            if (r == 0) continue;
                            while (l < rp && l < fp && RR[rp - l - 1] == ref(fp - l - 1)) ++l;
                            if (l + r >= L) { vi.push_back({rp - l, rp + r - 1, doc}); ++st.nver; }
                        }
                        st.tver += rz::now() - tv;
                    }
                    string d = "{"; { u64 nl = 0; append_docs(pos, Z, docs, d, nl); st.nlisted += nl; }
                    d += "} "; fr.push_back({left, right, d});
                }
            }
            st.tquery += rz::now() - t2;
            if (yR >= nwin) break;
            u64 nx = yR >= L - w ? yR - (L - w) : 0;
            x = nx > x ? nx : x + 1;
        }
        if (VFY) {                                   // MEMs = maximal verified intervals; genera = documents whose match contains them
            std::vector<std::pair<u64, u64>> iv;
            for (auto &x : vi) iv.push_back({x.s, x.e});
            std::sort(iv.begin(), iv.end()); iv.erase(std::unique(iv.begin(), iv.end()), iv.end());
            for (auto &p : iv) {
                bool in = false;
                for (auto &q : iv) if (q != p && q.first <= p.first && p.second <= q.second) { in = true; break; }
                if (in) continue;
                docs.clear();
                for (auto &x : vi) if (x.s <= p.first && p.second <= x.e) docs.push_back(x.d);
                std::sort(docs.begin(), docs.end()); docs.erase(std::unique(docs.begin(), docs.end()), docs.end());
                snprintf(buf, sizeof buf, "[%lu,%lu] {", r0 + p.first, r0 + p.second); line += buf;
                for (u64 q = 0; q < docs.size(); ++q) { snprintf(buf, sizeof buf, q ? ",%lu" : "%lu", docs[q]); line += buf; }
                line += "} "; ++st.nfeat; st.totlen += p.second - p.first + 1; st.nlisted += docs.size();
            }
            fr.clear();
        }
        // keep span-maximal MEMs only (a MEM found earlier can have its span inside a later one's)
        for (u64 p = 0; p < fr.size(); ++p) {
            bool in = false;
            for (u64 q = 0; q < fr.size() && !in; ++q)
                if (q != p && fr[q].l <= fr[p].l && fr[p].r <= fr[q].r && (fr[q].l < fr[p].l || fr[p].r < fr[q].r || q < p)) in = true;
            if (in) continue;
            snprintf(buf, sizeof buf, "[%lu,%lu] ", r0 + fr[p].l, r0 + fr[p].r); line += buf; line += fr[p].d;
            ++st.nfeat; st.totlen += fr[p].r - fr[p].l + 1;
            if (list) st.nlisted += std::count(fr[p].d.begin(), fr[p].d.end(), ',') + 1;
            if (check) found.push_back({fr[p].l, fr[p].r});
        }
        if (check) {                                 // brute force: every digested MEM, with its span
            std::vector<std::pair<u64, u64>> all;
            u64 prev_a = 0; std::vector<u64> A(nd);
            for (u64 bb = 0; bb < nd; ++bb) {        // leftmost start of a match ending at bb
                St c = E_.full(); u64 j = bb + 1;
                while (j > 0 && E_.step(c, (unsigned char)dg[j - 1])) --j;
                A[bb] = j;                           // == bb + 1 if dg[bb] does not occur
            }
            for (u64 bb = 0; bb < nd; ++bb) {
                if (A[bb] > bb) continue;
                if (bb + 1 < nd && A[bb + 1] <= A[bb]) continue;      // not right-maximal
                u64 aa = A[bb];
                u64 left = std::lower_bound(lo.begin(), lo.end(), aa) - lo.begin();
                u64 yR = std::upper_bound(hi.begin(), hi.end(), bb) - hi.begin();
                u64 right = yR >= nwin ? n - 1 : yR + w - 2;
                if (right >= left && right - left + 1 >= L) all.push_back({left, right});
            }
            (void)prev_a;
            {   // keep only span-maximal ones (spans not contained in another MEM's span)
                std::vector<std::pair<u64, u64>> mx;
                for (auto &p : all) { bool in = false; for (auto &q : all) if (q != p && q.first <= p.first && p.second <= q.second) { in = true; break; } if (!in) mx.push_back(p); }
                all.swap(mx);
            }
            if (all != found) {
                static thread_local u64 bad = 0;
                if (++bad <= 5) { fprintf(stderr, "BML CHECK MISMATCH read %lu: bml", st.nreads); for (auto &p : found) fprintf(stderr, " [%lu,%lu]", p.first, p.second);
                    fprintf(stderr, " brute"); for (auto &p : all) fprintf(stderr, " [%lu,%lu]", p.first, p.second); fprintf(stderr, "\n"); }
            }
            static thread_local u64 runs = 0; ++runs;
        }
    }
    if (st.nfeat == nfeat0) ++st.nuncl;
}

// ---- Boyer-Moore-Li on the plain DNA (no digest), threshold L in bases ---------------------------
// candidate x: backward-search R[x..x+L) from its right end; on failure at position j (R[j..x+L) does not
// occur) every start in [x, j] fails too, so jump to j + 1; on success extend left and right (reverse
// complement) to the MEM R[a..e) and continue from the first start whose window ends beyond it, e - L + 1.
// This finds every MEM of length at least L.
template <class E>
void classify_bml_dna(E &E_, rz::index &Z, const string &read, u64 L, bool list, string &line, Stats &st) {
    typedef typename E::St St;
    static thread_local std::vector<St> ivs, iv2, xiv;
    static thread_local std::vector<u64> pos, docs;
    char buf[96];
    static const bool check = getenv("RZ_BML_CHECK") != nullptr;
    u64 nfeat0 = st.nfeat, i = 0, n0 = read.size();
    std::vector<std::pair<u64, u64>> segs;           // (start, length): maximal A/C/G/T runs, or pseudo-MEMs (RZ_KEBAB)
    while (i < n0) {
        while (i < n0 && rz::dg_b2(read[i]) < 0) ++i;
        u64 r0 = i;
        while (i < n0 && rz::dg_b2(read[i]) >= 0) ++i;
        u64 n = i - r0;
        KB_total += n;
        if (KBON) {     // real KeBaB: maximal substrings all of whose k-mers pass the filter (timed as search)
            double tk = rz::now(); u64 K = KBF.k;
            if (n < K) { st.tsearch += rz::now() - tk; continue; }
            static thread_local std::vector<uint64_t> codes; codes.assign(n - K + 1, 0);
            KBF.kmers((const unsigned char *)read.data() + r0, n, [&](uint64_t p, uint64_t c) { codes[p] = c; });
            for (u64 p = 0; p < codes.size() && p < 8; ++p) KBF.prefetch(codes[p]);
            u64 s0 = 0; bool open = false;
            for (u64 p = 0; p < codes.size(); ++p) {
                if (p + 8 < codes.size()) KBF.prefetch(codes[p + 8]);
                bool ok = KBF.has(codes[p]);
                if (ok && !open) { s0 = p; open = true; }
                if (!ok && open) { segs.push_back({r0 + s0, p + K - 1 - s0}); open = false; }
            }
            if (open) segs.push_back({r0 + s0, n - s0});
            st.tsearch += rz::now() - tk; continue;
        }
        if (!KEBAB_K || n < KEBAB_K) { segs.push_back({r0, n}); continue; }
        // ideal KeBaB (exact k-mer membership, i.e. a filter without false positives; not timed, not counted
        // as steps): split the run into maximal substrings all of whose k-mers occur in the database
        u64 K = KEBAB_K, s0 = 0; bool open = false;
        for (u64 p = 0; p + K <= n; ++p) {
            St c = E_.full(); bool ok = true;
            for (u64 j = p + K; j > p; --j) if (!E_.step(c, (unsigned char)read[r0 + j - 1])) { ok = false; break; }
            if (ok && !open) { s0 = p; open = true; }
            if (!ok && open) { segs.push_back({r0 + s0, p + K - 1 - s0}); open = false; }
        }
        if (open) segs.push_back({r0 + s0, n - s0});
    }
    u64 kb_end = 0;                                  // pseudo-MEMs overlap: count the bases of their union
    for (auto sg : segs) {
        u64 r0 = sg.first, n = sg.second;
        if (n < L) continue;
        KB_kept += r0 + n > kb_end ? r0 + n - std::max(r0, kb_end) : 0; kb_end = std::max(kb_end, r0 + n);
        const char *R = read.data() + r0;
        std::vector<std::pair<u64, u64>> found;
        static thread_local std::vector<char> tsel;               // RZ_TRIM: minimizer positions (all tied minima) of this run
        if (TRIM_K && n >= TRIM_K) {
            u64 k = TRIM_K, W = TRIM_W - TRIM_K + 1, nk = n - k + 1;
            std::vector<u64> h(nk); uint32_t mask = (uint32_t)((1ull << (2 * k)) - 1), fw = 0;
            for (u64 t = 0; t < n; ++t) {
                fw = ((fw << 2) | (uint32_t)rz::dg_b2(R[t])) & mask;
                if (t + 1 >= k) { uint32_t rc = rz::dg_rc_code(fw, (int)k); h[t + 1 - k] = rz::dg_mix(fw < rc ? fw : rc); }
            }
            tsel.assign(nk, 0);
            for (u64 y = 0; y + W <= nk; ++y) {
                u64 mn = ~0ull; for (u64 q = y; q < y + W; ++q) mn = std::min(mn, h[q]);
                for (u64 q = y; q < y + W; ++q) if (h[q] == mn) tsel[q] = 1;
            }
        }
        u64 x = 0;
        while (x + L <= n) {
            double t1 = rz::now();
            St c = E_.full(); u64 j = x + L;
            if (E::isrz || HYB) { ivs.clear(); ivs.push_back(c); }
            auto &G = gtab<St>(); bool useG = G.t && !HYB && (!E::isrz || TAGX);   // the rz-index's grids need every suffix's interval
            if (useG && L >= (u64)G.t) {
                St c2 = c; if (G.get((const unsigned char *)R + x + L - G.t, c2)) { c = c2; j = x + L - G.t; ++st.steps; }
            }
            bool fail = false;
            while (j > x) { ++st.steps; if (!E_.step(c, (unsigned char)R[j - 1])) { fail = true; break; } --j; if (E::isrz || HYB) ivs.push_back(c); }
            if (fail) { st.tsearch += rz::now() - t1; x = j; continue; }      // R[j-1..x+L) does not occur
            u64 b = x + L;                                                       // ivs: suffixes of R[..b)
            while (j > 0) { ++st.steps; if (!E_.step(c, (unsigned char)R[j - 1])) break; --j; if (E::isrz || HYB) ivs.push_back(c); }
            u64 a = j, e = a;
            St xr = E_.full();
            if (E::isrz || HYB) { xiv.clear(); xiv.push_back(xr); }
            if (useG && e + G.t <= n) {                                      // rc(R[a..a+t)) by one lookup
                unsigned char rcb[32]; for (int q = 0; q < G.t; ++q) rcb[q] = rz::comp((unsigned char)R[e + G.t - 1 - q]);
                St x2 = xr; if (G.get(rcb, x2)) { xr = x2; e += G.t; ++st.steps; }
            }
            while (e < n) { ++st.steps; if (!E_.step(xr, rz::comp((unsigned char)R[e]))) break; ++e; if (E::isrz || HYB) xiv.push_back(xr); }
            u64 m = e - a;                                                       // MEM R[a..e)
            double t2;
            if constexpr (E::isrz) {
              if (TAGX) {
                t2 = rz::now(); st.tsearch += t2 - t1;
                u64 oa = a, oe = e;                  // reported interval (read-run coordinates, half-open)
                St li = xr; bool skip = false;
                if (TRIM_K) {                        // interior: full phrases whose boundaries every occurrence of the MEM shares
                    u64 k = TRIM_K, W = TRIM_W - TRIM_K + 1;
                    long zlo = (long)(a + W - 1), zhi = (long)e - (long)k - (long)W + 1;   // boundaries determined inside R[a..e)
                    long b1 = -1, b2 = -1;
                    for (long q = std::max(zlo, 0L); q <= zhi && q < (long)tsel.size(); ++q) if (tsel[q]) { if (b1 < 0) b1 = q; b2 = q; }
                    if (b1 < 0 || b2 == b1) skip = true;
                    else {
                        oa = (u64)b1; oe = (u64)b2 + k;
                        St c = E_.full(); bool ok = true;
                        for (u64 q = oe; q > oa; --q) if (!E_.step(c, (unsigned char)R[q - 1])) { ok = false; break; }
                        if (!ok) skip = true; else li = c;
                        ++TR_mems; TR_trim += (e - a) - (oe - oa);
                    }
                    if (skip) ++TR_skip;
                    t2 = rz::now();                  // the emulation's extra search is not timed
                }
                if (TRIM_FIXON) {                    // RZ_TRIMFIX=t: drop t bases at each end
                    if (e - a <= 2 * TRIM_FIX) skip = true;
                    else { oa = a + TRIM_FIX; oe = e - TRIM_FIX; St c = E_.full(); for (u64 q = oe; q > oa; --q) E_.step(c, (unsigned char)R[q - 1]); li = c; }
                    t2 = rz::now();
                }
                if (skip) { x = e - L + 1; continue; }
                snprintf(buf, sizeof buf, "[%lu,%lu] {", r0 + oa, r0 + oe - 1); line += buf;
                tag_answer(E::sp(li), E::ep(li), line, st.nlisted);
                line += "} ";
              } else {
                const std::vector<St> *suf = &ivs;
                if (e != b) {
                    St y = E_.full(); iv2.clear(); iv2.push_back(y);
                    for (u64 q = e; q > a; --q) { ++st.steps; E_.step(y, (unsigned char)R[q - 1]); iv2.push_back(y); }
                    suf = &iv2;
                }
                t2 = rz::now(); st.tsearch += t2 - t1;
                Z.fs.resize(m + 1); Z.fe.resize(m + 1); Z.xs.resize(m + 1); Z.xe.resize(m + 1);
                for (u64 q = 0; q <= m; ++q) { Z.fs[q] = (*suf)[m - q].sp; Z.fe[q] = (*suf)[m - q].ep; Z.xs[q] = xiv[q].sp; Z.xe[q] = xiv[q].ep; }
                auto r = Z.query_with_both(m);
                if (r.first == rz::NONE) { ++st.nfail; st.tquery += rz::now() - t2; x = e - L + 1; continue; }
                snprintf(buf, sizeof buf, "[%lu,%lu] {%lu,%lu} ", r0 + a, r0 + e - 1, Z.genome_of(r.first), Z.genome_of(r.second)); line += buf;
              }
            } else {
                t2 = rz::now(); st.tsearch += t2 - t1;
                if (!list) {
                    auto r = E_.ends(xr);
                    snprintf(buf, sizeof buf, "[%lu,%lu] {%lu,%lu} ", r0 + a, r0 + e - 1, Z.genome_of(r.first), Z.genome_of(r.second)); line += buf;
                } else if (std::pair<u64, u64> rr; HYB && hyb_lca<E>(E_, Z, xr, ivs, iv2, xiv, m, b - 1, e - 1, a, [&](u64 q) { return (unsigned char)R[q]; }, st, rr)) {
                    snprintf(buf, sizeof buf, "[%lu,%lu] <%lu,%lu> ", r0 + a, r0 + e - 1, Z.genome_of(rr.first), Z.genome_of(rr.second)); line += buf;
                } else {
                    E_.list(xr, pos);
                    snprintf(buf, sizeof buf, "[%lu,%lu] {", r0 + a, r0 + e - 1); line += buf;
                    append_docs(pos, Z, docs, line, st.nlisted);
                    line += "} ";
                }
            }
            st.tquery += rz::now() - t2;
            ++st.nfeat; st.totlen += m;
            if (check) found.push_back({a, e - 1});
            x = e - L + 1;
        }
        if (check) {                                 // brute force: every MEM of length >= L
            std::vector<std::pair<u64, u64>> all; std::vector<u64> A(n);
            for (u64 bb = 0; bb < n; ++bb) { St c = E_.full(); u64 j = bb + 1; while (j > 0 && E_.step(c, (unsigned char)R[j - 1])) --j; A[bb] = j; }
            for (u64 bb = 0; bb < n; ++bb) {
                if (A[bb] > bb) continue;
                if (bb + 1 < n && A[bb + 1] <= A[bb]) continue;
                if (bb - A[bb] + 1 >= L) all.push_back({A[bb], bb});
            }
            if (all != found) { static u64 bad = 0; if (++bad <= 5) fprintf(stderr, "BML CHECK MISMATCH read %lu (%zu vs %zu MEMs)\n", st.nreads, found.size(), all.size()); }
        }
    }
    if (st.nfeat == nfeat0) ++st.nuncl;
}

static thread_local u64 PX_text = 0, PX_phr = 0, PX_map = 0;   // search steps (forward-backward)

// -W: forward-backward (Li 2012) instead of BML: find every MEM, left to right (forward extension from the
// start of a MEM, then backward search from one past its end for the start of the next), but list only the
// MEMs of at least L bases.  Uses the lookup table (-F) when given.
static bool FB = false;
template <class E>
u64 ext_left_plain(E &E_, const unsigned char *R, u64 y, typename E::St &I) {
    I = E_.full(); u64 j = y;
    auto &G = gtab<typename E::St>();
    if (G.t && y >= (u64)G.t) { typename E::St c2 = I; if (G.get(R + y - G.t, c2)) { I = c2; j = y - G.t; } }
    while (j > 0 && E_.step(I, R[j - 1])) { --j; ++PX_text; }
    return j;
}
template <class E>
void classify_fb(E &E_, rz::index &Z, const string &read, u64 L, string &line, Stats &st) {
    typedef typename E::St St;
    static thread_local string rc;
    char buf[96]; u64 nfeat0 = st.nfeat, i = 0, n0 = read.size();
    while (i < n0) {
        while (i < n0 && rz::dg_b2(read[i]) < 0) ++i;
        u64 r0 = i; while (i < n0 && rz::dg_b2(read[i]) >= 0) ++i;
        u64 m = i - r0; if (m < L) continue;
        double t1 = rz::now();
        const unsigned char *R = (const unsigned char *)read.data() + r0;
        rc.resize(m); for (u64 t = 0; t < m; ++t) rc[t] = rz::comp((unsigned char)R[m - 1 - t]);
        const unsigned char *RC = (const unsigned char *)rc.data();
        auto ext = [&](const unsigned char *S, u64 y, St &I) { return ext_left_plain(E_, S, y, I); };
        u64 s0 = 0;
        while (s0 < m) {
            St I, J;
            u64 e = m - ext(RC, m - s0, J);                 // R[s0..e) is the MEM starting at s0
            if (e <= s0) { ++s0; continue; }
            if (e - s0 >= L) {
                double t2 = rz::now(); st.tsearch += t2 - t1;
                snprintf(buf, sizeof buf, "[%lu,%lu] {", r0 + s0, r0 + e - 1); line += buf;
                tag_answer(J.sp, J.ep, line, st.nlisted);
                line += "} ";
                t1 = rz::now(); st.tquery += t1 - t2;
                ++st.nfeat; st.totlen += e - s0;
            }
            if (e >= m) break;
            u64 s2 = ext(R, e + 1, I);                    // leftmost start of a match ending at R[e]
            s0 = s2 > e ? e + 1 : s2;
        }
        st.tsearch += rz::now() - t1;
    }
    st.steps = PX_text + PX_phr;
    if (st.nfeat == nfeat0) ++st.nuncl;
}

static bool HPC = false;                     // -U: homopolymer-compressed index and reads
static void hpc_read(const string &s, string &h, std::vector<uint32_t> &ofs) {   // ofs[i] = read position of h[i]
    h.clear(); ofs.clear();
    for (size_t i = 0; i < s.size(); ++i) if (i == 0 || s[i] != s[i - 1]) { h.push_back(s[i]); ofs.push_back((uint32_t)i); }
    ofs.push_back((uint32_t)s.size());
}
static void hpc_remap(string &line, const std::vector<uint32_t> &ofs) {          // [a,b] in h -> its span in the read
    string o; o.reserve(line.size() + 16); size_t i = 0;
    while (i < line.size()) {
        if (line[i] == '[') {
            size_t c = line.find(',', i), e = line.find(']', i);
            if (c != string::npos && e != string::npos && c < e) {
                u64 a = std::stoull(line.substr(i + 1, c - i - 1)), b = std::stoull(line.substr(c + 1, e - c - 1));
                o += "[" + std::to_string(ofs[a]) + "," + std::to_string(ofs[b + 1] - 1) + "]"; i = e + 1; continue;
            }
        }
        o += line[i++];
    }
    line.swap(o);
}
static int NJ = 1;                           // -j: threads
static std::mutex CNT_m;
static u64 T_TR_mems = 0, T_TR_trim = 0, T_TR_skip = 0, T_KB_total = 0, T_KB_kept = 0, T_NHYB = 0, T_PX_text = 0, T_PX_phr = 0, T_PX_map = 0, T_TAG_lf = 0;
static void flush_counters() {                   // add this thread's counters to the totals
    std::lock_guard<std::mutex> g(CNT_m);
    T_TR_mems += TR_mems; T_TR_trim += TR_trim; T_TR_skip += TR_skip; T_KB_total += KB_total; T_KB_kept += KB_kept; T_NHYB += NHYB;
    T_PX_text += PX_text; T_PX_phr += PX_phr; T_PX_map += PX_map; if (TAGX) T_TAG_lf += TAGX->lf_steps;
    TR_mems = TR_trim = TR_skip = KB_total = KB_kept = NHYB = PX_text = PX_phr = PX_map = 0; if (TAGX) TAGX->lf_steps = 0;
}
static void add_stats(Stats &a, const Stats &b) {
    a.nreads += b.nreads; a.nfeat += b.nfeat; a.totlen += b.totlen; a.nfail += b.nfail; a.steps += b.steps; a.nuncl += b.nuncl;
    a.nlisted += b.nlisted; a.nocc += b.nocc; a.nver += b.nver; a.tsearch += b.tsearch; a.tquery += b.tquery; a.tver += b.tver;
}

int main(int argc, char **argv) {
    string kbffile, auxfile, csafile, mapfile, srfile, rixfile, vfyfile, nseqfile, tagfile; u64 minlen = 1, bmlL = 0; int opt; bool mem = false, list = false;
    const char *usage = "usage: rz-classify [-j threads] [-M | -L bases [-l]] [-x aux] [-m minlen] [-C csa] [-B map] [-S sri -R rix] index.rz reads.fq out.listings\n";
    while ((opt = getopt(argc, argv, "MWUlx:m:C:B:S:R:L:V:H:A:N:T:F:j:K:")) != -1) {
        if (opt == 'U') { HPC = true; continue; }
        if (opt == 'x') auxfile = optarg; else if (opt == 'm') minlen = std::stoull(optarg); else if (opt == 'M') mem = true; else if (opt == 'W') FB = true; else if (opt == 'l') list = true; else if (opt == 'L') bmlL = std::stoull(optarg); else if (opt == 'V') vfyfile = optarg; else if (opt == 'H') HYB_T = std::stoull(optarg); else if (opt == 'A') HYB_A = std::stod(optarg); else if (opt == 'N') nseqfile = optarg; else if (opt == 'T') tagfile = optarg; else if (opt == 'F') FT_t = std::stoi(optarg); else if (opt == 'j') NJ = std::max(1, std::stoi(optarg)); else if (opt == 'C') csafile = optarg; else if (opt == 'B') mapfile = optarg; else if (opt == 'S') srfile = optarg; else if (opt == 'R') rixfile = optarg; else if (opt == 'K') kbffile = optarg;
        else { fputs(usage, stderr); return 1; }
    }
    if (argc - optind != 3) { fputs(usage, stderr); return 1; }
    if (FT_t && (!bmlL || FT_t > (rz::bytemap().on ? 3 : 13))) { fprintf(stderr, "-F needs -L, and t <= 13 (t <= 3 with -B)\n"); return 1; }
    if (HPC && (!bmlL || FB || !mapfile.empty())) { fprintf(stderr, "-U needs -L, without -W or -B\n"); return 1; }
    if (FB && (!bmlL || tagfile.empty() || rz::bytemap().on)) { fprintf(stderr, "-W needs -L and -T, without -B\n"); return 1; }
    if (const char *kb = getenv("RZ_KEBAB")) KEBAB_K = std::stoull(kb);
    if (!kbffile.empty()) {
        if (!bmlL || KEBAB_K) { fprintf(stderr, "-K needs -L, and not RZ_KEBAB\n"); return 1; }
        if (!KBF.load(kbffile)) { fprintf(stderr, "cannot load %s\n", kbffile.c_str()); return 1; }
        if (KBF.k > bmlL) { fprintf(stderr, "-K: the filter's k = %u exceeds L\n", KBF.k); return 1; }
        KBON = true; fprintf(stderr, "KeBaB filter %s: k=%u, h=%u, %.3f GB\n", kbffile.c_str(), KBF.k, KBF.h, KBF.bytes() / 1e9);
    }
    COUNTS = getenv("RZ_COUNTS") != nullptr;
    if (const char *lf = getenv("RZ_IVLOG")) { IVLOG = fopen(lf, "wb"); if (!IVLOG) { fprintf(stderr, "cannot write %s\n", lf); return 1; } atexit([] { fclose(IVLOG); }); }
    if (const char *tf = getenv("RZ_TRIMFIX")) { TRIM_FIX = std::stoull(tf); TRIM_FIXON = true; }
    if (const char *tr = getenv("RZ_TRIM")) { TRIM_K = std::stoull(tr); const char *c = strchr(tr, ','); TRIM_W = c ? std::stoull(c + 1) : 11; }
    // the rz-index's grids are needed only for its own LCA queries (no -T or -S, or -H/-A, or -x)
    const bool grids = !bmlL || (tagfile.empty() && srfile.empty()) || HYB_T || HYB_A > 0 || !auxfile.empty();
    rz::index Z; { std::ifstream in(argv[optind], std::ios::binary); Z.load(in, grids); }
    if (grids && !Z.has_grids) { fprintf(stderr, "%s has no grids (built with rz-build -G): use -L with -T or -S\n", argv[optind]); return 1; }
    rz::aux_index A;
    if (!auxfile.empty()) { std::ifstream in(auxfile, std::ios::binary); A.load(in); Z.aux = &A; }
    rz::rlcsa_bwt X;
    if (!csafile.empty()) { X.load(csafile); Z.csa = &X; }
    if (!mapfile.empty() && !rz::bytemap().load(mapfile)) { fprintf(stderr, "cannot load map %s\n", mapfile.c_str()); return 1; }
    if (!srfile.empty() && rixfile.empty()) { fprintf(stderr, "-S needs -R index.rix\n"); return 1; }
    if (bmlL && rz::bytemap().on && bmlL < (u64)rz::bytemap().w) { fprintf(stderr, "-L on a digest needs L >= w\n"); return 1; }
    if (list && (!bmlL || (srfile.empty() && tagfile.empty()))) { fprintf(stderr, "-l needs -L and -S or -T\n"); return 1; }
    static rz::tag_index TX;
    if (!tagfile.empty()) {
        if (!bmlL || !srfile.empty()) { fprintf(stderr, "-T needs -L and no -S\n"); return 1; }
        if (!TX.load(tagfile, !(list && COUNTS))) {   // counting needs no RMQ
            fprintf(stderr, "cannot load %s\n", tagfile.c_str()); return 1; }
        TX.bwt = &Z.bwt; TAGX = &TX; TAGLIST = list;
    }
    if (TAGX) list = false;                          // the rz engine does the search; listing comes from the tags
    HYB = HYB_T || HYB_A > 0;
    if (HYB && !list) { fprintf(stderr, "-H/-A need -l\n"); return 1; }
    if (HYB_A > 0) {
        if (nseqfile.empty()) { fprintf(stderr, "-A needs -N sequences-per-document file\n"); return 1; }
        std::ifstream nf(nseqfile); u64 c; NSEQ.assign(1, 0);
        while (nf >> c) NSEQ.push_back(NSEQ.back() + c);
        if (NSEQ.size() != Z.nstr / 2 + 1) { fprintf(stderr, "-N: %zu documents, index has %lu\n", NSEQ.size() - 1, Z.nstr / 2); return 1; }
    }
    static rz::vfy_index VV;
    if (!vfyfile.empty()) {
        if (!list || !rz::bytemap().on) { fprintf(stderr, "-V needs -l and -B\n"); return 1; }
        if (!VV.load(vfyfile)) { fprintf(stderr, "cannot load %s\n", vfyfile.c_str()); return 1; }
        VFY = &VV;
    }
    if (NJ > 1 && (!bmlL || (srfile.empty() && tagfile.empty()) || VFY || HYB)) { fprintf(stderr, "-j needs -L with -T or -S, without -V, -H or -A\n"); return 1; }
    std::ifstream fq(argv[optind + 1]);
    FILE *out = fopen(argv[optind + 2], "w");
    if (!fq || !out) { fprintf(stderr, "cannot open reads or output\n"); return 1; }

    Stats st; double WALL = 0;
    auto loop = [&](auto &E_) {
        if (FT_t) {
            auto &G = gtab<typename std::decay_t<decltype(E_)>::St>();
            double t0 = rz::now(); G.build(E_, FT_t, rz::bytemap().on ? 256 : 4); FT_bytes = G.bytes();
            fprintf(stderr, "lookup table: t=%d, %.1f MB, built in %.1f s\n", FT_t, G.bytes() / 1e6, rz::now() - t0);
        }
        // reads are taken in chunks of CH by NJ workers and written back in input order
        const size_t CH = 256; std::mutex in_m, out_m; std::condition_variable out_cv; u64 rseq = 0, wseq = 0; bool eof = false;
        std::vector<Stats> tst(NJ);
        auto work = [&](int tid) {
            Stats &ts = tst[tid]; std::vector<string> nm, sq, ln; string h, s, plus, q;
            for (;;) {
                u64 my;
                {
                    std::lock_guard<std::mutex> g(in_m); nm.clear(); sq.clear();
                    while (!eof && nm.size() < CH) {
                        if (!(std::getline(fq, h) && std::getline(fq, s) && std::getline(fq, plus) && std::getline(fq, q))) { eof = true; break; }
                        for (auto &c : s) c = (char)toupper((unsigned char)c);
                        size_t sp_ = h.find_first_of(" \t");
                        nm.push_back(h.substr(1, sp_ == string::npos ? string::npos : sp_ - 1)); sq.push_back(s);
                    }
                    my = rseq++;
                }
                ln.resize(nm.size());
                for (size_t i = 0; i < nm.size(); ++i) {
                    ++ts.nreads; ln[i].clear();
                    if (bmlL) {
                        bool l2 = std::decay_t<decltype(E_)>::isrz ? false : list;
                        if (rz::bytemap().on) classify_bml(E_, Z, sq[i], bmlL, l2, ln[i], ts);
                        else if (FB) { if constexpr (std::decay_t<decltype(E_)>::isrz) classify_fb(E_, Z, sq[i], bmlL, ln[i], ts); }
                        else if (HPC) {
                            static thread_local string hs; static thread_local std::vector<uint32_t> ofs;
                            double t0 = rz::now(); hpc_read(sq[i], hs, ofs); ts.tsearch += rz::now() - t0;
                            classify_bml_dna(E_, Z, hs, bmlL, l2, ln[i], ts);
                            hpc_remap(ln[i], ofs);
                        }
                        else classify_bml_dna(E_, Z, sq[i], bmlL, l2, ln[i], ts);
                    } else {
                        if (rz::bytemap().on) sq[i] = rz::digest_mapped(sq[i]);
                        classify_read(E_, Z, sq[i], mem, minlen, ln[i], ts);
                    }
                }
                {
                    std::unique_lock<std::mutex> lk(out_m); out_cv.wait(lk, [&] { return wseq == my; });
                    for (size_t i = 0; i < nm.size(); ++i) fprintf(out, ">%s\n%s\n", nm[i].c_str(), ln[i].c_str());
                    ++wseq;
                }
                out_cv.notify_all();
                if (nm.empty()) break;
            }
            flush_counters();
        };
        double w0 = rz::now();
        if (NJ == 1) work(0);
        else { std::vector<std::thread> th; for (int t = 0; t < NJ; ++t) th.emplace_back(work, t); for (auto &x : th) x.join(); }
        WALL = rz::now() - w0;
        for (auto &ts : tst) add_stats(st, ts);
    };
    if (srfile.empty()) {
        RzEng E_{Z};
        loop(E_);
    }
    else {
        rz::rlbwt SB;
        { std::ifstream in(rixfile, std::ios::binary); u64 nn; in.read((char *)&nn, 8); SB.load(in); }
        if (csafile.empty()) {
            rz::rl_nav NV(SB); rz::srindex_rl SR; { std::ifstream in(srfile, std::ios::binary); SR.load(in); } SR.mv = &NV;
            SrEng<rz::rl_nav> E_{SR}; loop(E_);
        } else {
            rz::csa_nav NV(SB, X); rz::srindex_csa SR; { std::ifstream in(srfile, std::ios::binary); SR.load(in); } SR.mv = &NV;
            SrEng<rz::csa_nav> E_{SR}; loop(E_);
        }
    }
    fclose(out);
    if (TAGX) fprintf(stderr, "tags s=%lu %s, LF steps %lu | ", TAGX->s, TAGLIST ? "list" : "lca", T_TAG_lf);
    if (bmlL) fprintf(stderr, "L=%lu %s unclassified=%lu docs/MEM=%.1f | ", bmlL, VFY ? "verified-list" : (list || TAGLIST) ? "list" : "lca", st.nuncl, st.nfeat ? (double)st.nlisted / st.nfeat : 0.0);
    if (HYB) fprintf(stderr, "hybrid T=%lu A=%.2f: %lu MEMs by LCA | ", HYB_T, HYB_A, T_NHYB);
    if (VFY) fprintf(stderr, "occurrences checked %lu, verified %lu (%.1f%%), verification %.2f s | ", st.nocc, st.nver, st.nocc ? 100.0 * st.nver / st.nocc : 0.0, st.tver);
    fprintf(stderr, "reads=%lu %s=%lu avg_len=%.1f failed=%lu steps/read=%.1f | search %.2f s, %s %.2f s | %.1f us/read (search incl.), %.2f us/feature query\n",
            st.nreads, bmlL ? "MEMs" : mem ? "MEMs" : "phrases", st.nfeat, st.nfeat ? (double)st.totlen / st.nfeat : 0.0, st.nfail,
            st.nreads ? (double)st.steps / st.nreads : 0.0, st.tsearch, srfile.empty() ? "rz ends" : "sr ends", st.tquery,
            1e6 * (st.tsearch + st.tquery) / (st.nreads ? st.nreads : 1), 1e6 * st.tquery / (st.nfeat ? st.nfeat : 1));
    fprintf(stderr, "threads=%d wall %.2f s, %.0f reads/s\n", NJ, WALL, WALL > 0 ? st.nreads / WALL : 0.0);
    if (FT_t) fprintf(stderr, "lookups: table %.1f MB\n", FT_bytes / 1e6);
    if (TRIM_K) fprintf(stderr, "trim k=%lu w=%lu: %lu MEMs trimmed by %.1f bases on average, %lu MEMs dropped\n", TRIM_K, TRIM_W, T_TR_mems, T_TR_mems ? (double)T_TR_trim / T_TR_mems : 0.0, T_TR_skip);
    if (KEBAB_K || KBON) fprintf(stderr, "kebab k=%lu: %.1f%% of read bases in pseudo-MEMs of length >= L\n", (u64)(KBON ? KBF.k : KEBAB_K), T_KB_total ? 100.0 * T_KB_kept / T_KB_total : 0.0);
    return st.nfail ? 2 : 0;
}
