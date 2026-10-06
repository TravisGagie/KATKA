// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Two-level indexing with closed syncmers (after Hong et al., AMB 2024): the text S is cut into phrases
// that start at canonical closed syncmers (a k-mer whose smallest s-mer, by the hash of its canonical form,
// is at its start or its end; a k-mer is one iff its reverse complement is).  A phrase runs from one
// syncmer start to the end of the next syncmer's k-mer (consecutive phrases overlap by k), and the last
// phrase of an A/C/G/T run runs to the separator X after it, inclusive.  Phrases get ids 1..d in
// lexicographic order (0 = no phrase).  We store
//   B   : n-bit bitvector marking the BWT rows of S whose suffixes start with a closed syncmer;
//   the parse's BWT, in the same row order restricted to B, as a RunRank (RLCSA) with counts C;
//   the dictionary as (64-bit hash -> id).
// Since closed syncmers depend only on their k-mers, every occurrence of a string that starts with one
// starts a phrase, so the BWT interval of such a string is all marked and maps to a parse interval by rank
// on B, and back by select.
#pragma once
#include <cstdint>
#include <vector>
#include <deque>
#include <algorithm>
#include <cstdio>
#include <string>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include "runrank.hpp"
#include "digest.hpp"

namespace rz {

inline uint64_t ph_hash(const unsigned char *s, uint64_t len) {
    uint64_t h = 1469598103934665603ULL ^ (len * 0x9e3779b97f4a7c15ULL);
    for (uint64_t i = 0; i < len; ++i) { h ^= s[i]; h *= 1099511628211ULL; }
    return dg_mix(h);
}

// sync[p] = 1 iff R[p..p+k) is a canonical closed syncmer (R over A/C/G/T only, length n)
inline void closed_syncmers(const unsigned char *R, uint64_t n, int k, int s, std::vector<char> &sync) {
    sync.assign(n, 0);
    if (n < (uint64_t)k) return;
    uint64_t ns = n - s + 1;
    static thread_local std::vector<uint64_t> h, dq;
    if (h.size() < ns) { h.resize(ns); dq.resize(ns); }
    uint32_t mask = (uint32_t)((1ull << (2 * s)) - 1), fw = 0;
    static thread_local std::vector<uint64_t> tab; static thread_local int tabs = -1;   // canonical hash of every s-mer code
    if (s <= 10 && tabs != s) {
        tab.resize(1ull << (2 * s));
        for (uint32_t c = 0; c <= mask; ++c) { uint32_t rc = dg_rc_code(c, s); tab[c] = dg_mix(c < rc ? c : rc); }
        tabs = s;
    }
    static const int8_t b2[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
    for (uint64_t t = 0; t < n; ++t) {
        fw = ((fw << 2) | (uint32_t)b2[R[t]]) & mask;
        if (t + 1 >= (uint64_t)s) {
            if (s <= 10) h[t + 1 - s] = tab[fw];
            else { uint32_t rc = dg_rc_code(fw, s); h[t + 1 - s] = dg_mix(fw < rc ? fw : rc); }
        }
    }
    uint64_t W = k - s + 1, hd = 0, tl = 0;           // s-mers per k-mer; monotone queue dq[hd..tl)
    for (uint64_t q = 0; q < ns; ++q) {
        while (tl > hd && h[dq[tl - 1]] > h[q]) --tl;
        dq[tl++] = q;
        if (q + 1 >= W) {
            uint64_t p = q + 1 - W;                   // k-mer at p covers s-mers p..q
            while (dq[hd] < p) ++hd;
            uint64_t mn = h[dq[hd]];
            if (h[p] == mn || h[q] == mn) sync[p] = 1;
        }
    }
}

struct parse_index {
    uint64_t n = 0, np = 0, d = 0; int k = 0, s = 0;
    sdsl::bit_vector B; sdsl::rank_support_v5<> Br; sdsl::select_support_mcl<> Bs;   // phrase starts are dense
    std::vector<uint64_t> C;                      // d+2: C[c] = phrase occurrences with id < c
    RunRank rr;                                   // parse BWT (ids of preceding phrases; 0 at run starts)
    std::vector<uint64_t> dh; std::vector<uint32_t> did;   // dictionary, sorted by hash

    std::vector<uint64_t> th; std::vector<uint32_t> tid; uint64_t tmask = 0;   // open addressing, built at load
    void build_table() {
        uint64_t sz = 1; while (sz < 2 * dh.size() + 2) sz <<= 1;
        th.assign(sz, 0); tid.assign(sz, 0); tmask = sz - 1;
        for (uint64_t j = 0; j < dh.size(); ++j) {
            uint64_t q = dh[j] & tmask; while (tid[q]) q = (q + 1) & tmask;
            th[q] = dh[j]; tid[q] = did[j];
        }
    }
    uint32_t id(uint64_t h) const {
        if (!tmask) { auto it = std::lower_bound(dh.begin(), dh.end(), h); return (it != dh.end() && *it == h) ? did[it - dh.begin()] : 0; }
        for (uint64_t q = h & tmask; tid[q]; q = (q + 1) & tmask) if (th[q] == h) return tid[q];
        return 0;
    }
    // backward step by phrase c on parse interval [l, r)
    bool step(uint64_t &l, uint64_t &r, uint32_t c) const {
        if (c == 0 || c > d) return false;
        uint64_t a = C[c] + rr.rank(c, l), b = C[c] + rr.rank(c, r);
        if (a >= b) return false;
        l = a; r = b; return true;
    }
    void to_parse(uint64_t sp, uint64_t ep, uint64_t &l, uint64_t &r) const { l = Br(sp); r = Br(ep); }
    void to_text(uint64_t l, uint64_t r, uint64_t &sp, uint64_t &ep) const { sp = Bs(l + 1); ep = Bs(r) + 1; }

    uint64_t bytes() const {
        return sdsl::size_in_bytes(B) + C.size() * 8 + rr.bytes() + dh.size() * 8 + did.size() * 4;
    }
    void save(const std::string &f) const {
        FILE *o = fopen(f.c_str(), "wb"); fwrite("RZPIX001", 1, 8, o);
        fwrite(&n, 8, 1, o); fwrite(&np, 8, 1, o); fwrite(&d, 8, 1, o); fwrite(&k, 4, 1, o); fwrite(&s, 4, 1, o);
        RunRank::wv(o, C); rr.save(o); RunRank::wv(o, dh); RunRank::wv(o, did); fclose(o);
        std::ofstream b(f + ".B", std::ios::binary); B.serialize(b);
    }
    bool load(const std::string &f) {
        FILE *i = fopen(f.c_str(), "rb"); char mg[8];
        if (!i || fread(mg, 1, 8, i) != 8 || std::string(mg, 8) != "RZPIX001") { if (i) fclose(i); return false; }
        if (fread(&n, 8, 1, i) != 1 || fread(&np, 8, 1, i) != 1 || fread(&d, 8, 1, i) != 1 || fread(&k, 4, 1, i) != 1 || fread(&s, 4, 1, i) != 1) return false;
        RunRank::rv(i, C); rr.load(i); RunRank::rv(i, dh); RunRank::rv(i, did); fclose(i);
        std::ifstream b(f + ".B", std::ios::binary); if (!b) return false; B.load(b);
        sdsl::util::init_support(Br, &B); sdsl::util::init_support(Bs, &B);
        build_table();
        return true;
    }
};

}  // namespace rz
