#!/bin/bash
# Fetches and builds the dependencies at the versions used in the paper, and downloads SILVA 138.1.
#   deps/sdsl-lite, deps/Big-BWT, pfp-merge (only needed for collections split into several datasets),
#   cliffy-experiments and MicrobeMixer (Cliffy's scripts and read simulator), silva/exp1_data (SILVA).
# Also needed: g++ (C++17), cmake, zlib, seqtk, ART (art_illumina), GNU time, Python 3 with
# regex, requests, aiohttp, pandas and multiprocess.
set -e
export CMAKE_POLICY_VERSION_MINIMUM=3.5   # sdsl-lite and Big-BWT ask for CMake versions that CMake 4 rejects
cd "$(dirname "$0")"; ROOT=$PWD
get() {  # get <url> <dir> <commit>
  [ -d "$2" ] || git clone "$1" "$2"
  git -C "$2" checkout -q "$3"
}
mkdir -p deps
get https://github.com/simongog/sdsl-lite deps/sdsl-lite c32874c
[ -d deps/sdsl/lib ] || (cd deps/sdsl-lite && ./install.sh "$ROOT/deps/sdsl")
get https://github.com/alshai/Big-BWT deps/Big-BWT 5b00650
make -C deps/Big-BWT
get https://gitlab.com/manzai/pfp-merge.git pfp-merge 5510113
get https://github.com/oma219/cliffy-experiments.git cliffy-experiments 566f647
get https://github.com/oma219/MicrobeMixer.git MicrobeMixer 6409f5d
# Cliffy itself is only needed to build Cliffy's indexes (silva/run_silva.sh with INDEXES=...):
#   get https://github.com/oma219/cliffy.git silva/cliffy 3763529 && silva/build_cliffy.sh
mkdir -p silva/exp1_data; cd silva/exp1_data
U=https://www.arb-silva.de/fileadmin/silva_databases/release_138.1
for f in Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz Exports/taxonomy/tax_slv_ssu_138.1.txt.gz Exports/taxonomy/tax_slv_ssu_138.1.tre.gz; do
  b=$(basename $f .gz); [ -s $b ] || { wget -q "$U/$f" && gunzip -f $(basename $f); }
done
ls -l
