// Build an RLCSA (explicit or Elias-Fano run lists) from a BWT file, for rz-bench -C / rz-classify -C.
// usage: rz-csabuild <bwt> <out.csa> [explicit|ef]
#include "rz_index.hpp"
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: rz-csabuild <bwt> <out.csa> [explicit|ef]\n"); return 1; }
    bool ef = argc > 3 && std::string(argv[3]) == "ef";
    double t0 = rz::now();
    rz::rlcsa_bwt X;
    X.build(argv[1], ef, [](const std::string &f, auto g) { rz::rlbwt::runs(f, g); });
    X.save(argv[2]);
    printf("n=%lu r=%lu %s RLCSA bytes=%lu time=%.1fs\n", X.n, X.R, ef ? "Elias-Fano" : "explicit", X.bytes(), rz::now() - t0);
}
