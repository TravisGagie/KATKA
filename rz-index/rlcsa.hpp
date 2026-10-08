// RLCSA backend for the rz-index's backward searches: the rank structures of Hybrid-FM-CSA
// (runrank.hpp = explicit per-character run lists, an uncompressed RLCSA; efrank.hpp = Elias-Fano run
// lists, a compressed RLCSA after Brown, Gagie, Manzini, Navarro and Sciortino), built directly from
// the BWT file's runs, without the text.  A backward step is a binary search (or Elias-Fano
// predecessor query) over the runs of one character, so its cost depends on that character's number
// of runs, not on the alphabet size.  It gives exactly the same intervals as rlbwt::extend.
// The explicit version stores positions in 32 bits if n < 2^32 and in 64 bits otherwise (`wide`); compiled
// with -DRZ_NO64 (make NO64=1), it supports only n < 2^32, as before, and never uses the 64-bit code.
#pragma once
#include "runrank.hpp"
#include "efrank.hpp"
#include <string>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace rz {

struct rlcsa_bwt {
    uint64_t n = 0, R = 0, C[257];
    bool ef = false, wide = false;
    RunRank rr;
    RunRank64 rr64;
    EFRank er;

    // one pass to count, one to fill; `runs` is rlbwt::runs (terminators already mapped to 0x01)
    template <class RUNS> void build(const std::string &bwtfile, bool use_ef, RUNS runs) {
        ef = use_ef;
        uint64_t cnt[256] = {0}, rc[256] = {0};
        runs(bwtfile, [&](uint64_t k, uint64_t s, uint64_t len, unsigned char c) { cnt[c] += len; rc[c]++; R = k + 1; n = s + len; });
        C[0] = 0; for (int c = 0; c < 256; ++c) C[c + 1] = C[c] + cnt[c];
        if (!ef) {
#ifdef RZ_NO64
            if (n >> 32) throw std::runtime_error("explicit RLCSA needs n < 2^32 when compiled with RZ_NO64; use the Elias-Fano version");
            wide = false;
#else
            wide = (n >> 32) || getenv("RZ_CSA64");     // RZ_CSA64=1: force 64-bit positions (testing)
#endif
            if (wide) build_explicit(rr64, bwtfile, rc, runs); else build_explicit(rr, bwtfile, rc, runs);
        } else {
            er.n = n; er.sigma = 256; er.nruns = R;
            std::vector<uint64_t> slot(257, 0), seen(256, 0);
            for (int c = 0; c < 256; ++c) slot[c + 1] = slot[c] + rc[c];
            std::vector<uint64_t> ks(R), fp(R);
            runs(bwtfile, [&](uint64_t, uint64_t s, uint64_t len, unsigned char c) {
                uint64_t k = slot[c]++; ks[k] = (uint64_t)c * n + s; fp[k] = C[c] + seen[c]; seen[c] += len; });
            er.KS.build(ks, 256 * n);
            { std::vector<uint64_t>().swap(ks); }
            er.FP.build(fp, n);
        }
    }
    template <class RR, class RUNS> void build_explicit(RR &rr, const std::string &bwtfile, const uint64_t *rc, RUNS &runs) {
        typedef decltype(rr.runs[0].start) W;
        rr.n = n; rr.sigma = 256; rr.nruns = R;
        rr.off.assign(257, 0);
        for (int c = 0; c < 256; ++c) rr.off[c + 1] = rr.off[c] + rc[c] + 1;
        rr.runs.assign(rr.off[256], typename RR::Run{0, 0});
        std::vector<uint64_t> pos(rr.off.begin(), rr.off.end() - 1), seen(256, 0);
        runs(bwtfile, [&](uint64_t, uint64_t s, uint64_t len, unsigned char c) {
            rr.runs[pos[c]++] = typename RR::Run{(W)s, (W)seen[c]}; seen[c] += len; });
        for (int c = 0; c < 256; ++c) rr.runs[pos[c]] = typename RR::Run{(W)n, (W)seen[c]};
        // jump tables, as in RunRank::build
        rr.joff.assign(257, 0); rr.shift.assign(256, 0); rr.jumps.clear();
        for (int c = 0; c < 256; ++c) {
            rr.joff[c] = rr.jumps.size();
            uint64_t Rc = rr.off[c + 1] - rr.off[c] - 1;
            if (Rc < RR::JUMP_MIN_RUNS) continue;
            uint64_t target = Rc / 8; uint8_t sh = 0;
            while ((n >> sh) > target) ++sh;
            rr.shift[c] = sh;
            uint64_t nb = (n >> sh) + 2, k = 0;
            const typename RR::Run *r = &rr.runs[rr.off[c]];
            for (uint64_t b = 0; b < nb; ++b) {
                uint64_t lim = b << sh;
                while (k < Rc && r[k].start < lim) ++k;
                rr.jumps.push_back((uint32_t)k);
            }
        }
        rr.joff[256] = rr.jumps.size();
    }
    template <class RR> static inline bool last_run_explicit(const RR &rr, unsigned char c, uint64_t i, uint64_t &start, uint64_t &len, uint64_t &cum) {
        const typename RR::Run *r = rr.runs.data() + rr.off[c];
        uint64_t R = rr.off[c + 1] - rr.off[c] - 1;
        if (R == 0) return false;
        uint64_t lo = 0, hi = R;
        if (rr.joff[c + 1] != rr.joff[c]) { const uint32_t *J = rr.jumps.data() + rr.joff[c]; uint64_t bk = i >> rr.shift[c]; lo = J[bk]; hi = J[bk + 1]; }
        while (lo < hi) { uint64_t mid = (lo + hi) >> 1; if (r[mid].start < i) lo = mid + 1; else hi = mid; }
        if (lo == 0) return false;
        start = r[lo - 1].start; cum = r[lo - 1].cum; len = r[lo].cum - cum;
        return true;
    }
    // the last run of c starting before position i (i >= 1): its start, length, and the number of c's
    // before it; false if there is none
    inline bool last_run(unsigned char c, uint64_t i, uint64_t &start, uint64_t &len, uint64_t &cum) const {
        if (i == 0) return false;
        if (!ef) return wide ? last_run_explicit(rr64, c, i, start, len, cum) : last_run_explicit(rr, c, i, start, len, cum);
        uint64_t base = (uint64_t)c * n, b0 = er.KS.count_less(base), a = er.KS.count_less(base + i);
        if (a == b0) return false;
        uint64_t j = a - 1, f, fn;
        start = er.KS.access(j) - base;
        if (j + 1 < er.nruns) er.FP.access2(j, f, fn); else { f = er.FP.access(j); fn = n; }
        cum = f - er.FP.access(b0); len = fn - f;
        return true;
    }
    inline uint64_t rank(unsigned char c, uint64_t i) const { return ef ? er.rank(c, i) : wide ? rr64.rank(c, i) : rr.rank(c, i); }
    // backward step with toehold information for the r-index and sr-index: [sp, ep) -> interval of cX;
    // covers = (BWT[ep-1] == c); otherwise y = position of the last c in [sp, ep), which ends a run
    inline bool step_th(uint64_t &sp, uint64_t &ep, unsigned char c, bool &covers, uint64_t &y) const {
        uint64_t st, len, cum;
        if (!last_run(c, ep, st, len, cum)) return false;
        uint64_t b = cum + std::min(len, ep - st);
        uint64_t a = rank(c, sp);
        if (a == b) return false;
        covers = ep - st <= len;
        y = st + len - 1;
        sp = C[c] + a; ep = C[c] + b;
        return true;
    }
    // [sp, ep) -> interval of cX
    inline void extend(uint64_t &sp, uint64_t &ep, unsigned char c) const {
        uint64_t a, b;
        if (ef) er.rank2(c, sp, ep, a, b); else if (wide) { a = rr64.rank(c, sp); b = rr64.rank(c, ep); } else { a = rr.rank(c, sp); b = rr.rank(c, ep); }
        sp = C[c] + a; ep = C[c] + b;
    }
    uint64_t bytes() const { return 8 * 259 + (ef ? er.bytes() : wide ? rr64.bytes() : rr.bytes()); }
    void save(const std::string &file) const {
        FILE *f = fopen(file.c_str(), "wb");
        fwrite("RZCSA001", 1, 8, f); fwrite(&n, 8, 1, f); fwrite(&R, 8, 1, f); fwrite(C, 8, 257, f);
        uint8_t e = ef ? 1 : wide ? 2 : 0; fwrite(&e, 1, 1, f);   // 2: explicit with 64-bit positions
        if (ef) er.save(f); else if (wide) rr64.save(f); else rr.save(f);
        fclose(f);
    }
    void load(const std::string &file) {
        FILE *f = fopen(file.c_str(), "rb"); char mg[8];
        if (!f || fread(mg, 1, 8, f) != 8 || memcmp(mg, "RZCSA001", 8)) throw std::runtime_error("bad RLCSA file");
        uint8_t e = 0;
        if (fread(&n, 8, 1, f) != 1 || fread(&R, 8, 1, f) != 1 || fread(C, 8, 257, f) != 257 || fread(&e, 1, 1, f) != 1)
            throw std::runtime_error("bad RLCSA file");
        ef = e == 1; wide = e == 2;
        if (ef) er.load(f); else if (wide) rr64.load(f); else rr.load(f);
        fclose(f);
    }
};

}  // namespace rz
