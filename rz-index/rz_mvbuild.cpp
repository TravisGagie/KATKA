// Build a move structure (Movi-style, one row per BWT run) from a BWT file.
// usage: rz-mvbuild <bwt> <out.mv>
#include "mv.hpp"
#include <fstream>
#include <chrono>
int main(int argc, char **argv) {
    if (argc != 3) { std::cerr << "usage: rz-mvbuild <bwt> <out.mv>\n"; return 1; }
    auto t0 = std::chrono::steady_clock::now();
    rz::move_bwt M;
    M.build(argv[1]);
    std::ofstream o(argv[2], std::ios::binary);
    uint64_t b = M.serialize(o);
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "n=" << M.n << " r=" << M.R << " bytes=" << b << " (" << (double)b / M.R << " per run) time=" << dt << "s\n";
}
