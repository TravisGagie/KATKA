#!/bin/bash
# Cliffy's half of the same-machine comparison with KATKA (see compare_lib.sh for what it does and its options).
#   cd silva && ./run_cliffy.sh      (resumable; summary: work/compare/cliffy-<host>/summary.txt)
# Defaults: CLIFFY="full_text minimizers", THEIRS=1.  Needs the reference and reads (run_silva.sh) and Cliffy
# (build_cliffy.sh); missing Cliffy indexes are built here with the settings of run_silva.sh (over a terabyte of
# temporary files for full_text: set TMPDIR_CLIFFY to a directory on a large disk).  Run before run_katka.sh.
TOOL=cliffy exec "$(dirname "$0")/compare_lib.sh" "$@"
