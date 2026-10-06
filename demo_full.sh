#!/bin/bash
# Full demo: every KATKA experiment in the paper, on the same reads (every 50th simulated read pair, 200,000 per
# region, downloaded from the GitHub release), with the same parameters.
#  1. dependencies, SILVA 138.1 and KATKA (setup_deps.sh, make);
#  2. the reference and every index (silva/build_katka.sh with ALL=1);
#  3. the experiments, in this order (each script says what it measures and where its results go):
#       run_bml.sh, run_bml_ud.sh      MEM threshold L, LCA vs listing, equal credit (digested, undigested)
#       run_freq.sh                    equal, proportional and normalized credit (L = 15, 30)
#       run_2l.sh                      two-level parse indexes (also builds them)
#       run_counts_all.sh              proportional credit for every L, times of every configuration with
#                                      occurrence counts, and the Pareto set (work/pareto_counts.txt)
#       run_ftab.sh, run_fb.sh, run_tag.sh   lookup tables, forward-backward, tag-array sampling
#       run_kebab.sh, run_trim.sh, run_hyb.sh, run_vfy_acc.sh   the negative results of the appendix
#       run_threads.sh, run_katka.sh   throughput with several threads; the default configuration
# Timings assume an otherwise idle machine.  Needs about 60 GB of disk and 16 GB of RAM, and takes a day or more.
# Resumable: rerun it after an interruption.  STEPS="..." runs only some of the experiment scripts.
set -e
ROOT=$(cd "$(dirname "$0")" && pwd); S=$ROOT/silva; W=$S/work
STEPS=${STEPS:-"run_bml.sh run_bml_ud.sh run_freq.sh run_2l.sh run_counts_all.sh run_ftab.sh run_fb.sh run_tag.sh run_kebab.sh run_trim.sh run_hyb.sh run_vfy_acc.sh run_threads.sh run_katka.sh"}
cd $ROOT
./setup_deps.sh
make -C rz-index -j4 DEPS=$ROOT/deps
ALL=1 silva/build_katka.sh
silva/get_reads.sh
mkdir -p $W/logs
for x in $STEPS; do
  [ -s $W/logs/$x.done ] && { echo "=== $x: already done"; continue; }
  echo "=== $(date '+%F %T') $x (log: silva/work/logs/$x.log)"
  if [ $x = run_katka.sh ]; then OURS="ud_L30 dg_L30" EVERY=50 silva/$x > $W/logs/$x.log 2>&1
  else silva/$x > $W/logs/$x.log 2>&1; fi
  tail -3 $W/logs/$x.log; date > $W/logs/$x.done
done
echo "=== all done; the Pareto set is in silva/work/pareto_counts.txt"
