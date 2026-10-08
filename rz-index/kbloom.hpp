// Blocked Bloom filter of the canonical k-mers (k <= 32) of a text over A/C/G/T, for KeBaB-style filtering
// (Brown et al.): a read position starts a k-mer that may occur in the text iff the filter says so, and a
// maximal run of positive k-mers covering at least L bases (a pseudo-MEM) contains every MEM of length at
// least L that overlaps it.  All h bits of a k-mer lie in one 512-bit block (one cache line), so a query
// costs one memory access whatever h is; this allows the optimal number of hash functions, and thus about
// half the space of a one-hash filter at the same false-positive rate (blocking costs a little).
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

namespace rz {

struct kbloom {
    uint32_t k = 0, h = 0;
    uint64_t nblocks = 0;
    std::vector<uint64_t> bits;                     // nblocks * 8 words

    static inline uint64_t mix(uint64_t x) {        // splitmix64 finalizer
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL; x ^= x >> 27; x *= 0x94d049bb133111ebULL; x ^= x >> 31; return x;
    }
    static inline int b2(unsigned char c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; } }
    inline uint64_t block_of(uint64_t hv) const { return (uint64_t)(((unsigned __int128)hv * nblocks) >> 64); }
    inline void add(uint64_t code) {
        uint64_t hv = mix(code), g = mix(hv ^ 0x9e3779b97f4a7c15ULL); uint64_t *B = &bits[block_of(hv) * 8];
        for (uint32_t i = 0; i < h; ++i) { uint32_t b = (g >> (9 * (i % 7))) & 511; if (i % 7 == 6) g = mix(g); B[b >> 6] |= 1ULL << (b & 63); }
    }
    inline bool has(uint64_t code) const {
        uint64_t hv = mix(code), g = mix(hv ^ 0x9e3779b97f4a7c15ULL); const uint64_t *B = &bits[block_of(hv) * 8];
        for (uint32_t i = 0; i < h; ++i) { uint32_t b = (g >> (9 * (i % 7))) & 511; if (i % 7 == 6) g = mix(g); if (!(B[b >> 6] >> (b & 63) & 1)) return false; }
        return true;
    }
    inline void prefetch(uint64_t code) const { __builtin_prefetch(&bits[block_of(mix(code)) * 8]); }
    // canonical codes of the k-mers of s[0..n): calls f(i, code) for each k-mer s[i..i+k) of A/C/G/T only
    template <class F> void kmers(const unsigned char *s, uint64_t n, F f) const {
        uint64_t fw = 0, rc = 0, len = 0, mask = k == 32 ? ~0ULL : (1ULL << (2 * k)) - 1;
        for (uint64_t i = 0; i < n; ++i) {
            int c = b2(s[i]);
            if (c < 0) { len = 0; fw = rc = 0; continue; }
            fw = ((fw << 2) | (uint64_t)c) & mask; rc = (rc >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1))); ++len;
            if (len >= k) f(i + 1 - k, fw < rc ? fw : rc);
        }
    }
    void init(uint64_t nkeys, double bits_per_key, uint32_t hh) {
        h = hh; nblocks = std::max<uint64_t>(1, (uint64_t)std::ceil(nkeys * bits_per_key / 512.0)); bits.assign(nblocks * 8, 0);
    }
    uint64_t bytes() const { return bits.size() * 8; }
    void save(const std::string &f) const {
        FILE *o = fopen(f.c_str(), "wb"); if (!o) throw std::runtime_error("cannot write " + f);
        fwrite("RZKBLM01", 1, 8, o); fwrite(&k, 4, 1, o); fwrite(&h, 4, 1, o); fwrite(&nblocks, 8, 1, o);
        fwrite(bits.data(), 8, bits.size(), o); fclose(o);
    }
    bool load(const std::string &f) {
        FILE *i = fopen(f.c_str(), "rb"); char mg[8]; if (!i) return false;
        if (fread(mg, 1, 8, i) != 8 || memcmp(mg, "RZKBLM01", 8) || fread(&k, 4, 1, i) != 1 || fread(&h, 4, 1, i) != 1 || fread(&nblocks, 8, 1, i) != 1) { fclose(i); return false; }
        bits.resize(nblocks * 8); bool ok = fread(bits.data(), 8, bits.size(), i) == bits.size(); fclose(i); return ok;
    }
};

}  // namespace rz
