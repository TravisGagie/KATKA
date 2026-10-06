#!/bin/bash
# Proportional credit end to end: accuracy for every L (run_freq.sh, freq scoring), times of every listing
# configuration with occurrence counts (run_pareto_ft.sh with COUNTS=1), and the resulting Pareto set.
# About 2.5 hours; run with the Claude app closed.  Output: work/pareto_counts.txt
S=$(cd "$(dirname "$0")" && pwd)
LS="15 20 30 40 50" MODES="freq" $S/run_freq.sh
COUNTS=1 $S/run_pareto_ft.sh
python3 $S/pareto_ft.py --counts > $S/work/pareto_counts.txt
head -40 $S/work/pareto_counts.txt
