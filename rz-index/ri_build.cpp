// Build the baseline r-index (Prezza's r-index, PFP fork by Mun et al.) for S from
// BWT(S$) and the SA samples written by `rz-build -a`, so both indexes share the same BWT.
//
// usage: ri-build-from <S.bwt> <samples-prefix> <out-prefix>
//   reads <samples-prefix>.ssa and <samples-prefix>.esa, writes <out-prefix>.ri

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#define private public
#include "r_index.hpp"
#undef private
#include <fstream>

int main(int argc, char **argv) {
    if (argc != 4) { std::cerr << "usage: ri-build-from <S.bwt> <samples-prefix> <out-prefix>\n"; return 1; }
    std::string bwt = argv[1], sp = argv[2], out = argv[3];
    ri::r_index<> idx;
    std::ifstream ifs(bwt, std::ios::binary);
    idx.bwt = ri::rle_string_sd(ifs);
    idx.r = idx.bwt.number_of_runs();
    ri::ulint n = idx.bwt.size();
    idx.build_F(ifs);
    std::vector<std::pair<ri::ulint, ri::ulint>> first;
    std::vector<ri::ulint> last;
    idx.read_run_starts(sp + ".ssa", n, first);
    idx.read_run_ends(sp + ".esa", n, last);
    if (first.size() != idx.r || last.size() != idx.r) { std::cerr << "sample count != r\n"; return 1; }
    idx.build_phi(first, last);
    std::ofstream o(out + ".ri", std::ios::binary);
    bool fast = false;
    o.write((char *)&fast, sizeof(fast));
    ri::ulint bytes = idx.serialize(o);
    std::cout << "r-index: n=" << n << " r=" << idx.r << " bytes=" << bytes << "\n";
    return 0;
}
