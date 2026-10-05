// Windowed greedy LZ77 (KKP / PSV-NSV) on S or on S^rev, in blocks, in O(W + B) space.
//
// usage: rz-lz77 [-r] [-w window] [-b block] [-t threads] <S> <out.ends>
//   -r         parse S^rev instead of S (read backwards; S^rev is never written to disk)
//   -w window  context kept before each block (default 128M characters)
//   -b block   characters parsed per block   (default 128M characters)
//   -t threads blocks parsed in parallel      (default 1; memory grows linearly)
// writes the phrase ends (0-based, in the coordinates of the parsed string) as uint64s.
//
// Every phrase is either a single character or copies a source starting earlier within
// the window, so the parse is left-referencing: that is all the rz-index needs (the
// leftmost occurrence of any pattern then contains a phrase end).  With the genomes of a
// species consecutive in S, a window of a few genomes loses very little against
// unrestricted LZ77.  Phrases never cross block boundaries.  Memory is about
// 13 (W + B) bytes per thread.

#include <divsufsort.h>
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <unistd.h>
#include <thread>
#include <mutex>

typedef uint64_t u64;

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static u64 parse_size(const char *s) {
    std::string a(s);
    u64 mul = 1;
    char last = a.back();
    if (last == 'K' || last == 'k') mul = 1ULL << 10;
    if (last == 'M' || last == 'm') mul = 1ULL << 20;
    if (last == 'G' || last == 'g') mul = 1ULL << 30;
    if (mul > 1) a.pop_back();
    return std::stoull(a) * mul;
}

// parse block [i0, end) of U with context [c0, i0); phrases never cross `end`
// (fixed block boundaries make blocks independent, so they can be parsed in parallel)
static void parse_block(const std::vector<unsigned char> &X, u64 c0, u64 i0, u64 end, std::vector<u64> &ends) {
    u64 nx = X.size();
    std::vector<int32_t> SA(nx);
    divsufsort(X.data(), SA.data(), (int32_t)nx);
    const uint32_t NIL = UINT32_MAX;
    std::vector<uint32_t> psv(nx, NIL), nsv(nx, NIL);
    {
        std::vector<int32_t> st;
        for (u64 j = 0; j < nx; ++j) {
            while (!st.empty() && st.back() > SA[j]) st.pop_back();
            psv[SA[j]] = st.empty() ? NIL : (uint32_t)st.back();
            st.push_back(SA[j]);
        }
        st.clear();
        for (u64 j = nx; j-- > 0;) {
            while (!st.empty() && st.back() > SA[j]) st.pop_back();
            nsv[SA[j]] = st.empty() ? NIL : (uint32_t)st.back();
            st.push_back(SA[j]);
        }
    }
    std::vector<int32_t>().swap(SA);
    u64 i = i0;
    while (i < end) {
        u64 p = i - c0, best = 0;
        for (uint32_t c : {psv[p], nsv[p]}) {
            if (c == NIL) continue;
            u64 l = 0;
            while (p + l < nx && X[c + l] == X[p + l]) ++l;
            if (l > best) best = l;
        }
        u64 len = best ? best : 1;
        ends.push_back(i + len - 1);
        i += len;
    }
}

int main(int argc, char **argv) {
    bool rev = false;
    u64 W = 128ULL << 20, B = 128ULL << 20, T = 1;
    int opt;
    const char *usage = "usage: rz-lz77 [-r] [-w window] [-b block] [-t threads] <S> <out.ends>\n";
    while ((opt = getopt(argc, argv, "rw:b:t:")) != -1) {
        if (opt == 'r') rev = true;
        else if (opt == 'w') W = parse_size(optarg);
        else if (opt == 'b') B = parse_size(optarg);
        else if (opt == 't') T = std::stoull(optarg);
        else { std::cerr << usage; return 1; }
    }
    if (argc - optind != 2) { std::cerr << usage; return 1; }
    if (W + B >= (1ULL << 31) - 1) { std::cerr << "window + block must be < 2^31\n"; return 1; }
    std::string in_name = argv[optind];
    FILE *probe = fopen(in_name.c_str(), "rb");
    if (!probe) { std::cerr << "cannot open " << in_name << "\n"; return 1; }
    fseeko(probe, 0, SEEK_END);
    u64 N = ftello(probe);
    fclose(probe);
    u64 nb = (N + B - 1) / B;
    std::vector<std::vector<u64>> res(nb);
    std::mutex mtx;
    u64 next = 0, done = 0;
    double t0 = now();
    auto worker = [&]() {
        FILE *in = fopen(in_name.c_str(), "rb");
        std::vector<unsigned char> X;
        while (true) {
            u64 b;
            { std::lock_guard<std::mutex> g(mtx); if (next >= nb) break; b = next++; }
            u64 i0 = b * B, end = std::min(N, i0 + B), c0 = i0 > W ? i0 - W : 0;
            X.resize(end - c0);
            // U = S or S^rev; read U[c0, end)
            if (!rev) { fseeko(in, c0, SEEK_SET); if (fread(X.data(), 1, end - c0, in) != end - c0) exit(2); }
            else {
                fseeko(in, N - end, SEEK_SET);
                if (fread(X.data(), 1, end - c0, in) != end - c0) exit(2);
                std::reverse(X.begin(), X.end());
            }
            parse_block(X, c0, i0, end, res[b]);
            std::lock_guard<std::mutex> g(mtx);
            ++done;
            std::cerr << "\r  blocks " << done << " / " << nb << " (" << now() - t0 << "s)" << std::flush;
        }
        fclose(in);
    };
    std::vector<std::thread> th;
    for (u64 t = 0; t < T; ++t) th.emplace_back(worker);
    for (auto &t : th) t.join();
    std::cerr << "\n";
    FILE *out = fopen(argv[optind + 1], "wb");
    u64 z = 0;
    for (auto &v : res) { fwrite(v.data(), 8, v.size(), out); z += v.size(); }
    fclose(out);
    std::cout << "N=" << N << " z=" << z << " direction=" << (rev ? "right-to-left" : "left-to-right")
              << " window=" << W << " block=" << B << " threads=" << T << " time=" << now() - t0 << "s\n";
    return 0;
}
