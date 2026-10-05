// Count the runs in the tag array (document = genus of each BWT row) of S, for document listing by
// Sadakane's method on runs.  usage: rz-tagruns <prefix>   (reads prefix.S and prefix.tbl)
// Suffixes of S are sorted with divsufsort; the order equals the BWT's (S$ with $ smallest), plus row 0 ($).
#include <divsufsort.h>
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: rz-tagruns prefix\n"); return 1; }
    std::string p = argv[1];
    auto t0 = std::chrono::steady_clock::now();
    FILE *f = fopen((p + ".S").c_str(), "rb"); fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    std::vector<unsigned char> S(n); if (fread(S.data(), 1, n, f) != (size_t)n) return 1; fclose(f);
    std::vector<uint64_t> start;                    // start of each document's pair of strings in S
    std::ifstream tb(p + ".tbl"); std::string line;
    while (std::getline(tb, line)) { std::istringstream ss(line); std::string x; for (int i = 0; i < 5; ++i) std::getline(ss, x, '\t'); start.push_back(std::stoull(x)); }
    std::vector<int32_t> SA(n);
    divsufsort(S.data(), SA.data(), (int32_t)n);
    uint64_t runs = 1, prev = ~0ull, bwtruns = 1; int last = -1;  // row 0 ($) is its own run
    uint64_t changes_within_bwt_runs = 0;
    for (long i = 0; i < n; ++i) {
        uint64_t q = SA[i];
        uint64_t d = std::upper_bound(start.begin(), start.end(), q) - start.begin() - 1;
        int c = q ? S[q - 1] : -2;                   // BWT character
        if (d != prev) { ++runs; if (c == last) ++changes_within_bwt_runs; }
        if (c != last) ++bwtruns;
        prev = d; last = c;
    }
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("n=%ld documents=%zu tag runs rho=%lu BWT runs r~%lu rho/r=%.2f (tag changes inside BWT runs: %lu) | %.0f s\n",
           n + 1, start.size(), runs, bwtruns, (double)runs / bwtruns, changes_within_bwt_runs, sec);
    return 0;
}
