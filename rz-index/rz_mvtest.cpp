// Check the move-structure backward search against the RLBWT: same intervals, same answers.
// usage: rz-mvtest <index.rz> <index.mv> <patterns>
#include "rz_index.hpp"
#include <fstream>
using namespace rz; using std::string; using std::vector;
int main(int argc, char **argv) {
    std::ifstream pin(argv[3], std::ios::binary);
    string header; std::getline(pin, header);
    auto get = [&](const string &key) { size_t p = header.find(key + "="); return std::stoull(header.substr(p + key.size() + 1)); };
    u64 N = get("number"), m = get("length");
    rz::index Z; { std::ifstream in(argv[1], std::ios::binary); Z.load(in); }
    move_bwt M; { std::ifstream in(argv[2], std::ios::binary); M.load(in); }
    if (M.n != Z.bwt.n || M.R != Z.bwt.R) { printf("size mismatch: mv n=%lu r=%lu, rlbwt n=%lu r=%lu\n", M.n, M.R, Z.bwt.n, Z.bwt.R); return 1; }
    u64 bad_iv = 0, bad_ans = 0;
    string P(m, 0);
    for (u64 q = 0; q < N; ++q) {
        pin.read(&P[0], m);
        auto a = Z.query(P);
        vector<u64> fs = Z.fs, fe = Z.fe, xs = Z.xs, xe = Z.xe;
        auto b = Z.query(P, M);
        if (a != b) ++bad_ans;
        if (a.first != NONE) {
            bool ok = fs == Z.fs && fe == Z.fe;
            for (u64 i = 0; i <= m && ok; ++i) {
                bool e1 = xs[i] >= xe[i], e2 = Z.xs[i] >= Z.xe[i];
                if (e1 != e2 || (!e1 && (xs[i] != Z.xs[i] || xe[i] != Z.xe[i]))) ok = false;
            }
            if (!ok) ++bad_iv;
        }
    }
    printf("patterns=%lu interval_mismatches=%lu answer_mismatches=%lu ff_steps/pattern=%.1f\n", N, bad_iv, bad_ans, (double)M.ff_steps / N);
    return bad_iv || bad_ans;
}
