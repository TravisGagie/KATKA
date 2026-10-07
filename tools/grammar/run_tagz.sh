#!/bin/bash
# LZ77 parse sizes of the run-length encoded tag array (tagz.cpp).  Fetches and compiles libsais, then runs
# tagz on silva/work/rz/bac.tagruns in three modes.  Needs about 9 GB of free memory for "both".
#   cd tools/grammar && ./run_tagz.sh       (results also in tagz.out)
set -e
cd "$(dirname "$0")"
[ -d libsais ] || git clone -q --depth 1 https://github.com/IlyaGrebnov/libsais.git
[ -x tagz ] || { gcc -O3 -c -Ilibsais/include libsais/src/libsais.c -o libsais.o && g++ -O3 -std=c++17 -Ilibsais/include tagz.cpp libsais.o -o tagz; }
F=${1:-../../silva/work/rz/bac.tagruns}
for m in genus len both; do /usr/bin/time -f "  (%e s, peak %M KB)" ./tagz $F $m; done 2>&1 | tee tagz.out
