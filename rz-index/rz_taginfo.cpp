// Sizes of the parts of a tag index (tag.hpp).  usage: rz-taginfo <file.tag>...
// Prints, in bytes: everything counting needs (run starts, sampling marks, tags) and the RMQ that listing also needs.
#include "tag.hpp"
using namespace rz;
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: rz-taginfo file.tag...\n"); return 1; }
    printf("file\ts\truns\tcounting_bytes\trmq_bytes\n");
    for (int a = 1; a < argc; ++a) {
        tag_index X; bool rmq = X.load(argv[a], true);
        if (!rmq && !X.load(argv[a], false)) { fprintf(stderr, "cannot load %s\n", argv[a]); return 1; }
        u64 all = X.bytes(), r = X.has_rmq ? sdsl::size_in_bytes(X.rmq) : 0;
        printf("%s\t%lu\t%lu\t%lu\t%lu\n", argv[a], X.s, X.rho, all - r, r);
    }
    return 0;
}
