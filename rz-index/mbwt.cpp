// Test helper: multi-string BWT of a list of dataset files, one terminator (0x00) per
// dataset, rotations sorted as in an eBWT with the strings' terminators ordered by content.
// Uses a suffix array of the concatenation (small inputs only, < 2 GB).
// usage: mbwt <datasets-list> <out.bwt>
#include <divsufsort.h>
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <numeric>
int main(int argc, char **argv) {
    if (argc != 3) { std::cerr << "usage: mbwt <datasets-list> <out.bwt>\n"; return 1; }
    std::ifstream l(argv[1]);
    std::vector<std::string> T;
    std::string f;
    while (std::getline(l, f)) { if (f.empty()) continue; std::ifstream in(f, std::ios::binary); std::stringstream ss; ss << in.rdbuf(); T.push_back(ss.str()); }
    size_t k = T.size();
    if (k > 60) { std::cerr << "too many datasets for this helper\n"; return 1; }
    std::vector<size_t> ord(k); std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return T[a] < T[b]; });   // $T_j order
    std::vector<unsigned char> rank(k);
    for (size_t i = 0; i < k; ++i) rank[ord[i]] = 1 + i;
    std::string D; std::vector<size_t> owner, start;
    for (size_t j = 0; j < k; ++j) { start.push_back(D.size()); D += T[j]; D.push_back((char)rank[j]); }
    start.push_back(D.size());
    std::vector<int32_t> SA(D.size());
    divsufsort((const unsigned char *)D.data(), SA.data(), (int32_t)D.size());
    std::string L(D.size(), 0);
    for (size_t i = 0; i < D.size(); ++i) {
        size_t p = SA[i];
        size_t j = std::upper_bound(start.begin(), start.end(), p) - start.begin() - 1;
        size_t prev = p == start[j] ? start[j + 1] - 1 : p - 1;           // cyclic within string j
        unsigned char c = D[prev];
        L[i] = c < 65 ? 0 : c;
    }
    std::ofstream o(argv[2], std::ios::binary); o.write(L.data(), L.size());
    return 0;
}
