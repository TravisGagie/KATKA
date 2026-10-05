// Verification structures for digested indexes: random access to the undigested reference (RLZ) and
// the DNA position of every digest symbol, so that occurrences of digested MEMs found in the digest
// can be checked (and extended) against the DNA.
//
// The digested reference has one document per genus file: records joined by separators, the digest
// computed on each maximal A/C/G/T run of length >= k, run digests joined by X (rz::digest_codes).
// For each document we keep its runs in order; for each run:
//   tstart : start of the run in T, the concatenation of all runs of all documents
//   len    : its length in bases
//   dstart : start of its digest within the document's digest
//   dlen   : the digest's length (symbols, after collapsing equal neighbours)
// Symbol positions: for each symbol (in run order) the first selected k-mer of its group, as an offset
// from the run start; stored as 4-bit gaps from the previous symbol's position (15 = escape: the gap is
// in `esc`), with an absolute sample (offset, escape rank) every 32 symbols.
// RLZ: T is parsed greedily against R (the longest run of each document); phrase f covers
// T[pstart[f] .. pstart[f+1]) and copies R[psrc[f] ..].  Phrases never cross runs.
#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

namespace rz {

struct vfy_index {
    std::string R;                                   // reference (A/C/G/T)
    std::vector<uint32_t> pstart, psrc;              // phrases (pstart has a final sentinel = |T|)
    std::vector<uint32_t> dfirst;                    // per document: first run (with a final sentinel)
    std::vector<uint32_t> dlg;                       // per document: digest length
    std::vector<uint32_t> tstart, rlen, dstart, dlen;// per run
    std::vector<uint64_t> sfirst;                    // per run: index of its first symbol (sentinel at end)
    std::vector<uint8_t> gaps;                       // 4-bit gaps, two per byte
    std::vector<uint16_t> esc;                       // escaped gaps
    std::vector<uint16_t> samp_pos;                  // per 32 symbols: absolute offset in the run
    std::vector<uint32_t> samp_esc;                  // per 32 symbols: escapes before it
    // the 32-symbol blocks are aligned to runs: each run starts a new block (sblock = first block of the run)
    std::vector<uint64_t> sblock;

    inline unsigned gap(uint64_t i) const { return (gaps[i >> 1] >> ((i & 1) * 4)) & 15; }
    // offset (in the run) of the first k-mer of symbol q of run j
    uint32_t sym_pos(uint32_t j, uint64_t q) const {
        uint64_t b = sblock[j] + q / 32, i0 = sfirst[j] + (q / 32) * 32, i1 = sfirst[j] + q;
        uint32_t p = samp_pos[b], e = samp_esc[b];
        for (uint64_t i = i0 + 1; i <= i1; ++i) { unsigned g = gap(i); if (g == 15) p += esc[e++]; else p += g; }
        if (i1 > i0) { /* the sample already covers i0 */ }
        return p;
    }
    // copy run j of T into out
    void run(uint32_t j, std::string &out) const {
        uint32_t a = tstart[j], n = rlen[j];
        out.resize(n);
        uint64_t f = std::upper_bound(pstart.begin(), pstart.end(), a) - pstart.begin() - 1;
        uint32_t i = 0;
        while (i < n) {
            uint32_t ps = pstart[f], pe = pstart[f + 1], off = a + i - ps, take = std::min(pe - ps - off, n - i);
            std::copy(R.begin() + psrc[f] + off, R.begin() + psrc[f] + off + take, out.begin() + i);
            i += take; ++f;
        }
    }
    template <class V> static void wv(FILE *f, const V &v) { uint64_t n = v.size(); fwrite(&n, 8, 1, f); fwrite(v.data(), sizeof(v[0]), n, f); }
    template <class V> static bool rv(FILE *f, V &v) { uint64_t n; if (fread(&n, 8, 1, f) != 1) return false; v.resize(n); return fread(v.data(), sizeof(v[0]), n, f) == n; }
    uint64_t bytes() const {
        return R.size() + 4 * (pstart.size() + psrc.size() + dfirst.size() + dlg.size() + tstart.size() + rlen.size() + dstart.size() + dlen.size())
               + 8 * (sfirst.size() + sblock.size()) + gaps.size() + 2 * esc.size() + 2 * samp_pos.size() + 4 * samp_esc.size();
    }
    void save(const std::string &fn) const {
        FILE *f = fopen(fn.c_str(), "wb"); fwrite("RZVFY001", 1, 8, f);
        uint64_t n = R.size(); fwrite(&n, 8, 1, f); fwrite(R.data(), 1, n, f);
        wv(f, pstart); wv(f, psrc); wv(f, dfirst); wv(f, dlg); wv(f, tstart); wv(f, rlen); wv(f, dstart); wv(f, dlen);
        wv(f, sfirst); wv(f, gaps); wv(f, esc); wv(f, samp_pos); wv(f, samp_esc); wv(f, sblock); fclose(f);
    }
    bool load(const std::string &fn) {
        FILE *f = fopen(fn.c_str(), "rb"); char mg[8];
        if (!f || fread(mg, 1, 8, f) != 8 || std::string(mg, 8) != "RZVFY001") { if (f) fclose(f); return false; }
        uint64_t n; bool ok = fread(&n, 8, 1, f) == 1; R.resize(n); ok = ok && fread(&R[0], 1, n, f) == n;
        ok = ok && rv(f, pstart) && rv(f, psrc) && rv(f, dfirst) && rv(f, dlg) && rv(f, tstart) && rv(f, rlen) && rv(f, dstart) && rv(f, dlen)
                && rv(f, sfirst) && rv(f, gaps) && rv(f, esc) && rv(f, samp_pos) && rv(f, samp_esc) && rv(f, sblock);
        fclose(f); return ok;
    }
};

}  // namespace rz
