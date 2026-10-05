// Movi-style move structure over a (multi-string) BWT, for backward search.
//
// One move row per maximal BWT run (so every row end is a run end):
//   start : first BWT row of the run (40 bits)
//   id    : move row containing LF(start)
//   dofs  : offset of LF(start) inside move row id
//   c     : run character (terminators 0x00-0x02 are mapped to 0x01, as in rlbwt)
// 14 bytes per run.  LF(k, o) = (id[k], dofs[k] + o), then "fast forward" over following rows.
// Backward search keeps interval ends as (move row, offset) pairs and scans run heads to
// reposition them (Movi's scanning mode), so no rank structure is needed.
#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

namespace rz {

struct move_bwt {
    typedef uint64_t u64;
    struct __attribute__((packed)) row {
        uint32_t s0; uint8_t s1; uint32_t id; uint32_t dofs; uint8_t c;
        inline u64 start() const { return s0 | ((u64)s1 << 32); }
    };
    u64 n = 0, R = 0;
    u64 C[257];
    std::vector<row> rows;               // R + 1 (sentinel with start = n)
    u64 ff_steps = 0;                    // statistics

    inline u64 start(u64 k) const { return rows[k].start(); }
    inline u64 len(u64 k) const { return rows[k + 1].start() - rows[k].start(); }
    inline unsigned char chr(u64 k) const { return rows[k].c; }
    inline u64 at(u64 k, u64 o) const { return rows[k].start() + o; }

    template <class F> static void scan(const std::string &file, F f) {
        FILE *fp = fopen(file.c_str(), "rb");
        if (!fp) { std::cerr << "cannot open " << file << "\n"; exit(1); }
        std::vector<unsigned char> buf(1 << 24);
        size_t got; u64 pos = 0;
        while ((got = fread(buf.data(), 1, buf.size(), fp)) > 0)
            for (size_t i = 0; i < got; ++i) { unsigned char c = buf[i] > 2 ? buf[i] : 1; f(pos++, c); }
        fclose(fp);
    }

    void build(const std::string &file) {
        u64 cnt[256] = {0};
        int last = -1; R = 0; n = 0;
        scan(file, [&](u64 pos, unsigned char c) { if ((int)c != last) { ++R; last = c; } cnt[c]++; n = pos + 1; });
        if (n >= (1ULL << 40) || R >= (1ULL << 32)) { std::cerr << "move_bwt: BWT too large\n"; exit(1); }
        C[0] = 0;
        for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + cnt[c];
        rows.assign(R + 1, row{0, 0, 0, 0, 0});
        u64 k = 0; last = -1;
        auto set_start = [&](u64 k, u64 s) { rows[k].s0 = (uint32_t)s; rows[k].s1 = (uint8_t)(s >> 32); };
        scan(file, [&](u64 pos, unsigned char c) {
            if ((int)c != last) { set_start(k, pos); rows[k].c = c; ++k; last = c; }
        });
        set_start(R, n);
        // LF of each run head: F-positions of runs of the same character increase with k,
        // so one pointer per character sweeps the move rows.
        u64 before[256] = {0};
        std::vector<u64> ptr(256, 0);
        auto run_containing = [&](u64 x) {
            u64 lo = 0, hi = R;                          // largest k with start(k) <= x
            while (hi - lo > 1) { u64 mid = (lo + hi) / 2; if (start(mid) <= x) lo = mid; else hi = mid; }
            return lo;
        };
        for (int c = 0; c < 256; ++c) if (cnt[c]) ptr[c] = run_containing(C[c]);
        for (u64 k = 0; k < R; ++k) {
            unsigned char c = rows[k].c;
            u64 f = C[c] + before[c];
            u64 &p = ptr[c];
            while (start(p + 1) <= f) ++p;
            u64 o = f - start(p);
            if (o >= (1ULL << 32)) { std::cerr << "move_bwt: run too long\n"; exit(1); }
            rows[k].id = (uint32_t)p; rows[k].dofs = (uint32_t)o;
            before[c] += len(k);
        }
    }

    // LF of row (k, o), in place
    inline void LF(u64 &k, u64 &o) {
        const row &r = rows[k];
        k = r.id; o = r.dofs + o;
        u64 l;
        while (o >= (l = len(k))) { o -= l; ++k; ++ff_steps; }
    }

    // Interval [(ks,os), (ke,oe)] (inclusive).  Returns false (and leaves the interval
    // unspecified) if cX does not occur.
    struct interval { u64 ks, os, ke, oe; };
    inline interval full() const { return {0, 0, R - 1, len(R - 1) - 1}; }
    inline bool step(interval &I, unsigned char c) {
        while (I.ks <= I.ke && rows[I.ks].c != c) { ++I.ks; I.os = 0; }
        if (I.ks > I.ke || (I.ks == I.ke && I.os > I.oe)) return false;
        if (rows[I.ke].c != c) { do --I.ke; while (rows[I.ke].c != c); I.oe = len(I.ke) - 1; }
        if (I.ks > I.ke || (I.ks == I.ke && I.os > I.oe)) return false;
        LF(I.ks, I.os); LF(I.ke, I.oe);
        return true;
    }
    // Same, also tracking the toehold: SA[last row] = SA[row (tk, to)] - td.
    inline bool step_toehold(interval &I, unsigned char c, u64 &tk, u64 &to, u64 &td) {
        while (I.ks <= I.ke && rows[I.ks].c != c) { ++I.ks; I.os = 0; }
        if (I.ks > I.ke || (I.ks == I.ke && I.os > I.oe)) return false;
        if (rows[I.ke].c != c) {
            do --I.ke; while (rows[I.ke].c != c);
            I.oe = len(I.ke) - 1;
            if (I.ks > I.ke || (I.ks == I.ke && I.os > I.oe)) return false;
            tk = I.ke; to = I.oe; td = 1;              // end of a BWT run with character c
        } else ++td;
        LF(I.ks, I.os); LF(I.ke, I.oe);
        return true;
    }

    inline u64 rows_sp(const interval &I) const { return at(I.ks, I.os); }
    inline u64 rows_ep(const interval &I) const { return at(I.ke, I.oe) + 1; }
    inline void first(const interval &I, u64 &k, u64 &o) const { k = I.ks; o = I.os; }
    inline void init_toehold(u64 &tk, u64 &to) const { tk = R - 1; to = len(R - 1) - 1; }
    u64 bytes() const { return 16 + 257 * 8 + rows.size() * sizeof(row); }
    u64 serialize(std::ostream &out) const {
        out.write((char *)&n, 8); out.write((char *)&R, 8); out.write((char *)C, 257 * 8);
        out.write((char *)rows.data(), rows.size() * sizeof(row));
        return bytes();
    }
    void load(std::istream &in) {
        in.read((char *)&n, 8); in.read((char *)&R, 8); in.read((char *)C, 257 * 8);
        rows.resize(R + 1);
        in.read((char *)rows.data(), rows.size() * sizeof(row));
    }
};

}  // namespace rz
