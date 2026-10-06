#!/bin/bash
# KATKA's half of the same-machine comparison with Cliffy (see compare_lib.sh for what it does and its options).
#   cd silva && ./run_katka.sh      (resumable; summary: work/compare/katka-<host>/summary.txt)
# Defaults: OURS="ud_L15 ud_L30 dg_L15 dg_L30" (ud = undigested: explicit RLCSA, tag array, 10-base table;
# dg = minimizer digests: explicit RLCSA, tag array, 2-symbol table; _L = minimum MEM length), MODE=freq.
# Needs the reference and reads (run_silva.sh) and our indexes (run_silva.sh, run_silva_dg.sh, build_indexes.sh).
# On the same machine as run_cliffy.sh, run this after it, so that both use the same number of threads.
TOOL=katka exec "$(dirname "$0")/compare_lib.sh" "$@"
