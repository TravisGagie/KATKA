// Sample random patterns of length m from S that contain no X (Pizza&Chili format).
// usage: rz-genpat <S> <m> <N> <out> [seed]
#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
int main(int argc, char **argv) {
    if (argc < 5) { std::cerr << "usage: rz-genpat <S> <m> <N> <out> [seed]\n"; return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { std::cerr << "cannot open " << argv[1] << "\n"; return 1; }
    fseeko(f, 0, SEEK_END);
    uint64_t n = ftello(f), m = std::stoull(argv[2]), N = std::stoull(argv[3]);
    std::mt19937_64 rng(argc > 5 ? std::stoull(argv[5]) : 42);
    std::uniform_int_distribution<uint64_t> d(0, n - m);
    std::ofstream out(argv[4], std::ios::binary);
    out << "# number=" << N << " length=" << m << " file=" << argv[1] << " forbidden=X\n";
    std::string p(m, 0);
    for (uint64_t i = 0, tries = 0; i < N; ) {
        if (++tries > 1000 * N + 1000000) { std::cerr << "too many rejections\n"; return 1; }
        fseeko(f, d(rng), SEEK_SET);
        if (fread(&p[0], 1, m, f) != m) continue;
        if (p.find('X') != std::string::npos) continue;
        out.write(p.data(), m); ++i;
    }
    return 0;
}
