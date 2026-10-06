// Copyright (C) 2026 Travis Gagie
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See the LICENSE file for details.
//
// Build the two-level parse index (parse.hpp) of a text S built by rz-prep (undigested).
// usage: rz-parsebuild <prefix> <k> <s> <out.pix>
//   reads prefix.S; suffix-sorts it with divsufsort (32-bit, so |S| < 2^31; about 5|S| bytes of memory).
// Checks that the marked suffixes, in BWT order, start with non-decreasing phrase ids, i.e. that the
// parse's suffix order agrees with the text's.
#include "parse.hpp"
#include <divsufsort.h>
#include <unordered_map>
#include <chrono>
#include <cstring>
using u64 = uint64_t;
static double secs(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: rz-parsebuild prefix k s out.pix\n"); return 1; }
    std::string p = argv[1]; int k = std::stoi(argv[2]), s = std::stoi(argv[3]);
    if (s < 1 || s >= k || s > 16) { fprintf(stderr, "need 1 <= s < k and s <= 16\n"); return 1; }
    auto T0 = std::chrono::steady_clock::now();
    FILE *f = fopen((p + ".S").c_str(), "rb"); if (!f) { perror("S"); return 1; }
    fseek(f, 0, SEEK_END); long N = ftell(f); fseek(f, 0, SEEK_SET);
    if (N >= (1L << 31) - 1) { fprintf(stderr, "S too long for 32-bit divsufsort\n"); return 1; }
    std::vector<unsigned char> S(N); if (fread(S.data(), 1, N, f) != (size_t)N) return 1; fclose(f);

    // 1. closed syncmers (phrase starts), run by run.  Phrase i runs from its syncmer start to the end of the
    // next syncmer's k-mer if that is in the same A/C/G/T run, and otherwise to the separator after the run.
    auto forEachPhrase = [&](auto cb) {           // cb(i, start, len, firstInRun)
        std::vector<char> sy; std::vector<u64> st; u64 i = 0, idx = 0;
        while (i < (u64)N) {
            while (i < (u64)N && rz::dg_b2((char)S[i]) < 0) ++i;
            u64 a = i; while (i < (u64)N && rz::dg_b2((char)S[i]) >= 0) ++i;
            if (i - a < (u64)k) continue;
            rz::closed_syncmers(S.data() + a, i - a, k, s, sy);
            st.clear(); for (u64 q = 0; q + k <= i - a; ++q) if (sy[q]) st.push_back(a + q);
            for (u64 j = 0; j < st.size(); ++j) {
                u64 len = (j + 1 < st.size()) ? st[j + 1] + k - st[j] : std::min<u64>(i + 1, N) - st[j];
                cb(idx++, st[j], len, j == 0);
            }
        }
    };
    sdsl::bit_vector sb(N, 0); u64 np = 0;
    forEachPhrase([&](u64, u64 st, u64, bool) { sb[st] = 1; ++np; });
    sdsl::rank_support_v5<> sbr(&sb);
    fprintf(stderr, "k=%d s=%d: %lu phrases (one per %.2f characters) (%.0f s)\n", k, s, np, (double)N / np, secs(T0));
    // 2. dictionary
    std::unordered_map<u64, std::pair<u64, uint32_t>> dict; dict.reserve(np / 8 + 16);
    u64 coll = 0;
    forEachPhrase([&](u64, u64 st, u64 len, bool) {
        u64 h = rz::ph_hash(S.data() + st, len);
        auto it = dict.find(h);
        if (it == dict.end()) dict.emplace(h, std::make_pair(st, (uint32_t)len));
        else if (it->second.second != len || memcmp(S.data() + it->second.first, S.data() + st, len)) ++coll;
    });
    if (coll) { fprintf(stderr, "%lu hash collisions between distinct phrases; aborting\n", coll); return 1; }
    std::vector<std::pair<u64, std::pair<u64, uint32_t>>> D(dict.begin(), dict.end()); dict.clear(); dict.rehash(0);
    std::sort(D.begin(), D.end(), [&](const auto &x, const auto &y) {
        u64 lx = x.second.second, ly = y.second.second; int c = memcmp(S.data() + x.second.first, S.data() + y.second.first, std::min(lx, ly));
        return c ? c < 0 : lx < ly; });
    u64 d = D.size();
    rz::parse_index X; X.n = N + 1; X.np = np; X.d = d; X.k = k; X.s = s;
    {
        std::vector<std::pair<u64, uint32_t>> hv(d);
        for (u64 j = 0; j < d; ++j) hv[j] = {D[j].first, (uint32_t)(j + 1)};
        std::sort(hv.begin(), hv.end());
        X.dh.resize(d); X.did.resize(d);
        for (u64 j = 0; j < d; ++j) { X.dh[j] = hv[j].first; X.did[j] = hv[j].second; }
    }
    u64 dictchars = 0; for (auto &e : D) dictchars += e.second.second;
    D.clear(); D.shrink_to_fit();
    fprintf(stderr, "%lu distinct phrases, %lu characters (%.0f s)\n", d, dictchars, secs(T0));
    // 3. phrase ids and run starts
    std::vector<uint32_t> pid(np); sdsl::bit_vector fir(np, 0);
    X.C.assign(d + 2, 0);
    forEachPhrase([&](u64 i, u64 st, u64 len, bool fr) {
        pid[i] = X.id(rz::ph_hash(S.data() + st, len)); fir[i] = fr;
        X.C[pid[i] + 1]++;
    });
    for (u64 c = 1; c < d + 2; ++c) X.C[c] += X.C[c - 1];
    // 4. suffix array; walk the rows
    std::vector<int32_t> SA(N); divsufsort(S.data(), SA.data(), (int32_t)N);
    fprintf(stderr, "suffix array (%.0f s)\n", secs(T0));
    std::vector<uint32_t> seq; seq.reserve(np);
    sdsl::bit_vector bb(N + 1, 0); u64 bad = 0; uint32_t last = 0;
    for (u64 row = 1; row <= (u64)N; ++row) {
        u64 q = SA[row - 1];
        if (!sb[q]) continue;
        bb[row] = 1;
        u64 i = sbr(q);                                // phrase index
        if (pid[i] < last) ++bad;
        last = pid[i];
        seq.push_back(fir[i] ? 0 : pid[i - 1]);
    }
    SA.clear(); SA.shrink_to_fit(); S.clear(); S.shrink_to_fit();
    if (bad) fprintf(stderr, "WARNING: %lu marked suffixes out of phrase order\n", bad);
    else fprintf(stderr, "marked suffixes are in phrase order\n");
    X.B = std::move(bb);
    sdsl::util::init_support(X.Br, &X.B); sdsl::util::init_support(X.Bs, &X.B);
    X.rr.build(seq, (uint32_t)(d + 1));
    fprintf(stderr, "parse BWT: %lu symbols, %lu runs; index %.3f GB (B %.3f, RLCSA %.3f, dictionary %.3f) (%.0f s)\n",
            np, X.rr.nruns, X.bytes() / 1e9, sdsl::size_in_bytes(X.B) / 1e9, (X.rr.bytes() + X.C.size() * 8) / 1e9,
            (X.dh.size() * 12) / 1e9, secs(T0));
    X.save(argv[4]);
    return bad ? 2 : 0;
}
