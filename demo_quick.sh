#!/bin/bash
# Quick demo: KATKA with its default parameters on 5,000 simulated read pairs per 16S region (from the repository).
#  1. fetches and builds the dependencies and SILVA 138.1 (setup_deps.sh) and KATKA (make);
#  2. builds the reference (one document per genus) and the default index (silva/build_katka.sh);
#  3. classifies the reads with MEMs of at least 30 bases, an explicit RLCSA with a 10-base lookup table, an
#     unsampled tag array and proportional credit, first with one thread and then with all cores, and scores them.
# Prints a summary (silva/work/compare/katka-<host>/summary.txt).  Needs about 20 GB of disk and 8 GB of RAM;
# building takes about half an hour.  On a CPU with performance and efficiency cores, PIN=0 runs the
# single-thread timing on CPU 0.  Resumable.
set -e
ROOT=$(cd "$(dirname "$0")" && pwd); S=$ROOT/silva; W=$S/work
cd $ROOT
./setup_deps.sh
make -C rz-index -j4 DEPS=$ROOT/deps
silva/build_katka.sh
mkdir -p $W/reads/sample2000
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  mkdir -p $W/reads/aquatic/$r
  [ -e $W/reads/aquatic/$r/${r}_seqtax.txt ] || cp $S/data/sample2000/${r}_seqtax.txt $W/reads/aquatic/$r/
  for m in 1 2; do f=$S/data/sample2000/${r}_mate_$m.fq.gz; [ -e $W/reads/sample2000/$(basename $f) ] || cp $f $W/reads/sample2000/; done
done
OURS="ud_L30" EVERY=2000 silva/run_katka.sh
