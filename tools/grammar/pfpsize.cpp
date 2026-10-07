// pfpsize: sizes of dictionary + parse decompositions of the run-length encoded tag array, one or two levels.
// usage: pfpsize prefix.tagruns mode scheme a b [levels]
//   mode:   both (len_1, g_1, len_2, g_2, ...; lengths and genera with disjoint token values) | genus | len
//   scheme: pfp  a=w b=p   phrase boundaries where the Karp-Rabin hash of the last w tokens is 0 mod p
//                          (prefix-free parsing: consecutive phrases overlap by the w-token trigger)
//           sync a=k b=s   phrase boundaries at closed syncmers: k-token windows whose smallest s-token hash is at
//                          the start or the end of the window (phrases overlap by k tokens, as with PFP)
//   levels: 1 (default) or 2 (the parse, a sequence of phrase IDs, is parsed again with the same scheme)
// Prints, per level: tokens, phrases |P|, distinct phrases |D|, total tokens in D (with and without the overlaps),
// and a space estimate: D stored plainly at the given bits per token plus P at ceil(log2 |D|) bits per phrase.
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <deque>
#include <unordered_map>
#include <algorithm>
#include <cmath>

static void die(const char *m) { fprintf(stderr, "pfpsize: %s\n", m); exit(1); }

struct Runs {
    std::string f; uint64_t n = 0, r = 0;
    explicit Runs(const std::string &file) : f(file) {
        FILE *h = fopen(f.c_str(), "rb"); if (!h) die("cannot open the .tagruns file");
        if (fread(&n, 8, 1, h) != 1 || fread(&r, 8, 1, h) != 1) die("bad header"); fclose(h);
    }
    template <class F> void each(F fn) const {
        FILE *a = fopen(f.c_str(), "rb"), *b = fopen(f.c_str(), "rb");
        fseek(a, 16, SEEK_SET); fseek(b, 16 + 4 * r, SEEK_SET);
        const uint64_t B = 1 << 20; std::vector<uint32_t> s(B + 1); std::vector<uint16_t> g(B);
        uint32_t next; if (fread(&next, 4, 1, a) != 1) die("read");
        for (uint64_t j0 = 0; j0 < r; j0 += B) {
            uint64_t m = std::min(B, r - j0);
            s[0] = next;
            if (j0 + m < r) { if (fread(&s[1], 4, m, a) != m) die("read"); next = s[m]; }
            else { if (m > 1 && fread(&s[1], 4, m - 1, a) != m - 1) die("read"); s[m] = (uint32_t)n; }
            if (fread(g.data(), 2, m, b) != m) die("read");
            for (uint64_t i = 0; i < m; ++i) fn((uint64_t)s[i + 1] - s[i], g[i]);
        }
        fclose(a); fclose(b);
    }
};

static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }

struct Hash128 { uint64_t a, b; bool operator==(const Hash128 &o) const { return a == o.a && b == o.b; } };
struct H128h { size_t operator()(const Hash128 &h) const { return h.a ^ (h.b * 0x9e3779b97f4a7c15ULL); } };

// Parses a token stream into phrases; feed tokens with push(), then finish().
struct Parser {
    std::string scheme; uint64_t A, Bp;
    // pfp state
    std::deque<uint64_t> win; uint64_t kr = 0, powA = 1; const uint64_t MOD = (1ULL << 61) - 1, BASE = 1000003;
    // sync state: last k token hashes, and a deque for sliding-window minima of s-mer hashes
    std::deque<uint64_t> toks;
    std::vector<uint64_t> cur;           // tokens of the current phrase
    std::unordered_map<Hash128, uint32_t, H128h> dict; std::vector<uint32_t> dlen;
    std::vector<uint32_t> parse;         // phrase IDs
    uint64_t ntok = 0, dict_tokens = 0;
    Parser(std::string s, uint64_t a, uint64_t b) : scheme(s), A(a), Bp(b) {
        for (uint64_t i = 0; i + 1 < A; ++i) powA = mulmod(powA, BASE);   // BASE^(w-1)
    }
    static uint64_t mulmod(uint64_t x, uint64_t y) { __uint128_t z = (__uint128_t)x * y; uint64_t lo = (uint64_t)(z & ((1ULL << 61) - 1)), hi = (uint64_t)(z >> 61); uint64_t r = lo + hi; if (r >= (1ULL << 61) - 1) r -= (1ULL << 61) - 1; return r; }
    void emit(size_t keep) {   // the current phrase ends here; the next starts with its last `keep` tokens
        Hash128 h{0x12345, 0x6789}; for (uint64_t t : cur) { h.a = mix(h.a ^ t); h.b = mix(h.b + t * 0x9e3779b97f4a7c15ULL); }
        h.a ^= cur.size();
        auto it = dict.find(h); uint32_t id;
        if (it == dict.end()) { id = dlen.size(); dict.emplace(h, id); dlen.push_back(cur.size()); dict_tokens += cur.size(); } else id = it->second;
        parse.push_back(id);
        std::vector<uint64_t> nxt(cur.end() - std::min(keep, cur.size()), cur.end()); cur.swap(nxt);
    }
    bool boundary(uint64_t t) {   // after pushing token t, is the window ending here a trigger?
        if (scheme == "pfp") {
            uint64_t v = (t * 0x9e3779b97f4a7c15ULL) % MOD + 1;
            if (win.size() == A) { uint64_t o = win.front(); win.pop_front(); kr = (kr + MOD - mulmod(o, powA)) % MOD; }
            kr = (mulmod(kr, BASE) + v) % MOD; win.push_back(v);
            return win.size() == A && kr % Bp == 0;
        } else {   // closed syncmer of k = A tokens with s = Bp
            toks.push_back(mix(t + 0x51));
            if (toks.size() > A) toks.pop_front();
            if (toks.size() < A) return false;
            uint64_t best = UINT64_MAX; size_t bi = 0;
            for (size_t i = 0; i + Bp <= A; ++i) {
                uint64_t h = 0x9ae16a3b2f90404fULL; for (size_t j = 0; j < Bp; ++j) h = mix(h ^ toks[i + j]);
                if (h < best) { best = h; bi = i; }
            }
            return bi == 0 || bi + Bp == A;
        }
    }
    void push(uint64_t t) {
        cur.push_back(t); ++ntok;
        if (boundary(t) && cur.size() > A) emit(A);
    }
    void finish() { if (!cur.empty()) emit(0); }
};

static void report(int level, const Parser &P, double bits_per_token) {
    uint64_t D = P.dlen.size(), Pn = P.parse.size();
    uint64_t overl = P.dict_tokens - std::min<uint64_t>(P.dict_tokens, (uint64_t)D * P.A);   // without the w-token overlaps
    double idbits = std::max(1.0, std::ceil(std::log2((double)std::max<uint64_t>(D, 2))));
    double bytes = (P.dict_tokens * bits_per_token + Pn * idbits) / 8;
    printf("  level %d: tokens %lu, phrases |P| %lu (%.1f tokens each), distinct |D| %lu, tokens in D %lu (%lu without overlaps)"
           ", estimate %.3f GB (D %.3f + P %.3f at %.0f bits/ID)\n",
           level, P.ntok, Pn, (double)P.ntok / Pn, D, P.dict_tokens, overl, bytes / 1e9,
           P.dict_tokens * bits_per_token / 8e9, Pn * idbits / 8e9, idbits);
}

int main(int argc, char **argv) {
    if (argc < 6) die("usage: pfpsize prefix.tagruns both|genus|len pfp|sync a b [levels]");
    std::string mode = argv[2], scheme = argv[3]; uint64_t a = atoll(argv[4]), b = atoll(argv[5]); int levels = argc > 6 ? atoi(argv[6]) : 1;
    if (scheme == "sync" && (b == 0 || b > a)) die("sync needs 1 <= s <= k");
    Runs R(argv[1]);
    Parser L1(scheme, a, b);
    double bpt = mode == "genus" ? 14 : mode == "len" ? 12 : 14;   // bits per token if D is stored plainly
    R.each([&](uint64_t len, uint16_t g) {
        if (mode == "both") { L1.push(len); L1.push((1ULL << 40) + g); }
        else if (mode == "genus") L1.push(g);
        else L1.push(len);
    });
    L1.finish();
    printf("%s, %s a=%lu b=%lu:\n", mode.c_str(), scheme.c_str(), a, b);
    report(1, L1, bpt);
    if (levels >= 2) {
        Parser L2(scheme, a, b);
        for (uint32_t id : L1.parse) L2.push(id);
        L2.finish();
        double idbits1 = std::ceil(std::log2((double)std::max<size_t>(L1.dlen.size(), 2)));
        report(2, L2, idbits1);
        double tot = (L1.dict_tokens * bpt + L2.dict_tokens * idbits1 + L2.parse.size() * std::ceil(std::log2((double)std::max<size_t>(L2.dlen.size(), 2)))) / 8e9;
        printf("  two levels: %.3f GB (D1 + D2 + P2)\n", tot);
    }
    return 0;
}
