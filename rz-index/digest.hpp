// Strand-symmetric minimizer digests (k <= 3), and the reverse-complement map on the digest alphabet.
//
// A k-mer is hashed by its canonical form (the smaller of its 2-bit code and that of its reverse
// complement), and in every window of w bases (w - k + 1 k-mers) ALL positions attaining the minimum
// hash are selected, so ties are broken symmetrically.  The selected k-mers, in text order and in the
// orientation of the strand, become one symbol each; runs of equal consecutive symbols are collapsed
// (as in SPUMONI 2).  Maximal runs of A/C/G/T are digested separately and joined with X.
//
// Hence digest(rc(G)) = rcsym(reverse(digest(G))), where rcsym maps the symbol of a k-mer to the
// symbol of its reverse complement: the digested collection stays closed under reverse complement.
//
// Symbols are bytes below 128 (pfp-merge reads bytes as signed chars), never 0..2 (reserved by
// Big-BWT / pfp-merge) and never A, C, G, T or X:  k = 3 -> 32..63 and 96..127,  k = 2 -> 16..31,
// k = 1 -> 12..15.  A, C, G, T and X keep their usual complements.
#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <cstdio>

namespace rz {

inline unsigned char dg_sym(int k, uint32_t f) {
    if (k == 3) return (unsigned char)(f < 32 ? 32 + f : 96 + (f - 32));
    return (unsigned char)(k == 2 ? 16 + f : 12 + f);
}
// (k, code) of a digest symbol; k = 0 if c is not one
inline void dg_unsym(unsigned char c, int &k, uint32_t &f) {
    if (c >= 32 && c < 64) { k = 3; f = c - 32; }
    else if (c >= 96 && c < 128) { k = 3; f = c - 96 + 32; }
    else if (c >= 16 && c < 32) { k = 2; f = c - 16; }
    else if (c >= 12 && c < 16) { k = 1; f = c - 12; }
    else k = 0;
}
inline uint32_t dg_rc_code(uint32_t f, int k) {       // 2-bit code of the reverse complement
    uint32_t r = 0;
    for (int i = 0; i < k; ++i) { r = (r << 2) | (3 - (f & 3)); f >>= 2; }
    return r;
}
// Byte-mapped digests (any k <= 4): each minimizer code that occurs gets its own byte (3..255, skipping X),
// the separator is X.  The map is written by rz-prep -B and loaded by the query tools; while loaded, it
// defines rc_symbol.
struct dg_bytemap {
    bool on = false;
    int k = 0, w = 0;
    unsigned char code2byte[257];      // codes 0..4^k-1 and the separator 4^k; 0 = code never seen
    unsigned char comp[256];           // byte -> byte of the reverse-complement k-mer
    void save(const std::string &f) const {
        FILE *o = fopen(f.c_str(), "wb"); fwrite("RZDGMAP1", 1, 8, o); fwrite(&k, 4, 1, o); fwrite(&w, 4, 1, o);
        fwrite(code2byte, 1, 257, o); fwrite(comp, 1, 256, o); fclose(o);
    }
    bool load(const std::string &f) {
        FILE *i = fopen(f.c_str(), "rb"); char mg[8];
        if (!i || fread(mg, 1, 8, i) != 8 || std::string(mg, 8) != "RZDGMAP1" || fread(&k, 4, 1, i) != 1 || fread(&w, 4, 1, i) != 1
            || fread(code2byte, 1, 257, i) != 257 || fread(comp, 1, 256, i) != 256) { if (i) fclose(i); return false; }
        fclose(i); on = true; return true;
    }
};
inline dg_bytemap &bytemap() { static dg_bytemap m; return m; }

// reverse complement of one symbol of either alphabet (DNA or digest); anything else maps to itself
inline unsigned char rc_symbol(unsigned char c) {
    if (bytemap().on) return bytemap().comp[c];
    switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: break; }
    int k; uint32_t f;
    dg_unsym(c, k, f);
    return k ? dg_sym(k, dg_rc_code(f, k)) : c;
}
inline uint64_t dg_mix(uint64_t x) {                   // splitmix64 finalizer
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
inline int dg_b2(char c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; } }

// digest of a string over A/C/G/T only
inline void digest_run(const char *s, uint64_t n, int k, int w, std::string &out) {
    if (n < (uint64_t)k) return;
    uint64_t nk = n - k + 1, W = (uint64_t)(w - k + 1);
    if (W > nk) W = nk;
    uint32_t mask = (1u << (2 * k)) - 1;
    std::vector<uint64_t> h(nk);
    std::vector<uint32_t> f(nk);
    uint32_t fw = 0;
    for (uint64_t i = 0; i < n; ++i) {
        fw = ((fw << 2) | (uint32_t)dg_b2(s[i])) & mask;
        if (i + 1 >= (uint64_t)k) {
            uint64_t p = i + 1 - k;
            uint32_t rc = dg_rc_code(fw, k);
            f[p] = fw;
            h[p] = dg_mix(fw < rc ? fw : rc);
        }
    }
    std::vector<char> sel(nk, 0);
    std::deque<uint64_t> dq;                           // positions with non-decreasing hashes
    for (uint64_t p = 0; p < nk; ++p) {
        while (!dq.empty() && h[dq.back()] > h[p]) dq.pop_back();
        dq.push_back(p);
        if (p + 1 >= W) {
            uint64_t s0 = p + 1 - W;
            while (dq.front() < s0) dq.pop_front();
            uint64_t mn = h[dq.front()];
            for (uint64_t q : dq) { if (h[q] != mn) break; sel[q] = 1; }   // all tied minima
        }
    }
    for (uint64_t p = 0; p < nk; ++p) if (sel[p]) {
        char c = (char)dg_sym(k, f[p]);
        if (out.empty() || out.back() != c) out.push_back(c);
    }
}
// digest of a genome in which non-ACGT characters (including contig separators) are X
inline std::string digest(const std::string &g, int k, int w) {
    std::string out, part;
    uint64_t i = 0, n = g.size();
    while (i < n) {
        while (i < n && dg_b2(g[i]) < 0) ++i;
        uint64_t j = i;
        while (j < n && dg_b2(g[j]) >= 0) ++j;
        if (j > i) {
            part.clear();
            digest_run(g.data() + i, j - i, k, w, part);
            if (!part.empty()) { if (!out.empty()) out.push_back('X'); out += part; }
        }
        i = j;
    }
    return out;
}
// The same digest as integer codes (the k-mer's 2-bit code, in the orientation of the strand), for any
// k <= 12, e.g. SPUMONI 2's k = 4.  Maximal A/C/G/T runs are digested separately and joined by `sep`.
inline void digest_codes(const std::string &g, int k, int w, uint32_t sep, std::vector<uint32_t> &out) {
    out.clear();
    uint64_t i = 0, n = g.size();
    while (i < n) {
        while (i < n && dg_b2(g[i]) < 0) ++i;
        uint64_t j = i;
        while (j < n && dg_b2(g[j]) >= 0) ++j;
        uint64_t len = j - i;
        if (len >= (uint64_t)k) {
            uint64_t nk = len - k + 1, W = (uint64_t)(w - k + 1);
            if (W > nk) W = nk;
            uint32_t mask = (uint32_t)((1ull << (2 * k)) - 1), fw = 0;
            std::vector<uint64_t> h(nk); std::vector<uint32_t> f(nk);
            for (uint64_t t = 0; t < len; ++t) {
                fw = ((fw << 2) | (uint32_t)dg_b2(g[i + t])) & mask;
                if (t + 1 >= (uint64_t)k) {
                    uint64_t p = t + 1 - k; uint32_t rc = dg_rc_code(fw, k);
                    f[p] = fw; h[p] = dg_mix(fw < rc ? fw : rc);
                }
            }
            std::vector<char> sel(nk, 0); std::deque<uint64_t> dq;
            for (uint64_t p = 0; p < nk; ++p) {
                while (!dq.empty() && h[dq.back()] > h[p]) dq.pop_back();
                dq.push_back(p);
                if (p + 1 >= W) {
                    uint64_t s0 = p + 1 - W;
                    while (dq.front() < s0) dq.pop_front();
                    uint64_t mn = h[dq.front()];
                    for (uint64_t q : dq) { if (h[q] != mn) break; sel[q] = 1; }
                }
            }
            bool first = true;
            for (uint64_t p = 0; p < nk; ++p) if (sel[p]) {
                if (first) { if (!out.empty()) out.push_back(sep); first = false; }
                if (out.empty() || out.back() != f[p]) out.push_back(f[p]);
            }
        }
        i = j;
    }
}
inline std::string digest_mapped(const std::string &g) {
    const dg_bytemap &M = bytemap();
    std::vector<uint32_t> c; uint32_t sep = 1u << (2 * M.k);
    digest_codes(g, M.k, M.w, sep, c);
    std::string out(c.size(), 0);
    for (size_t i = 0; i < c.size(); ++i) out[i] = (char)M.code2byte[c[i]];
    return out;
}
inline std::string rc_string(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (auto &c : r) c = (char)rc_symbol((unsigned char)c);
    return r;
}

}  // namespace rz
