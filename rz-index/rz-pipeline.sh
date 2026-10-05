#!/usr/bin/env bash
# End-to-end rz-index construction with pfp-merge, plus a small benchmark.
#
# usage: rz-pipeline.sh <list.tsv> <out-prefix>
# environment (defaults in brackets):
#   RZ       directory with the rz-index binaries              [directory of this script]
#   BIGBWT   Big-BWT's bigbwt script                           [bigbwt]
#   PFPMERGE pfp-merge checkout (with build/ made)              [../pfp-merge]
#   PREP_OPTS extra rz-prep options, e.g. -s for multi-genome FASTA files [none]
#   THREADS  threads for rz-lz77 and pfp-merge                  [4]
#   WINDOW   rz-lz77 window = block size                        [256M]
#   LENGTHS  pattern lengths to benchmark                       ["10 20 50 100 200"]
#   NPAT     patterns per length                                [1000]
set -euo pipefail
list=$1; out=$2
RZ=${RZ:-$(cd "$(dirname "$0")" && pwd)}
BIGBWT=${BIGBWT:-bigbwt}
PFPMERGE=${PFPMERGE:-$RZ/../pfp-merge}
PREP_OPTS=${PREP_OPTS:-}
THREADS=${THREADS:-4}
WINDOW=${WINDOW:-256M}
LENGTHS=${LENGTHS:-"10 20 50 100 200"}
NPAT=${NPAT:-1000}
outabs=$(cd "$(dirname "$out")" && pwd)/$(basename "$out")

echo "== S and one dataset per species"
"$RZ/rz-prep" -d $PREP_OPTS "$list" "$outabs"
mapfile -t ds < "$outabs.datasets"

echo "== BWT: Big-BWT per dataset, then pfp-merge"
for f in "${ds[@]}"; do "$BIGBWT" -t "$THREADS" "$f" > "$f.bigbwt.log" 2>&1; done
if [ "${#ds[@]}" -ge 2 ]; then
  (cd "$PFPMERGE" && ./PFPmerge.py -p "$THREADS" -o "$outabs.bwt" "${ds[@]}")
else
  cp "${ds[0]}.bwt" "$outabs.bwt"
fi

echo "== parses"
"$RZ/rz-lz77" -t "$THREADS" -w "$WINDOW" -b "$WINDOW" "$outabs.S" "$outabs.left"
"$RZ/rz-lz77" -r -t "$THREADS" -w "$WINDOW" -b "$WINDOW" "$outabs.S" "$outabs.right"

echo "== rz-index and baseline r-index"
"$RZ/rz-build" -a "$outabs.S" "$outabs.bwt" "$outabs.left" "$outabs.right" "$outabs.tbl" "$outabs"

echo "== benchmark"
for m in $LENGTHS; do
  "$RZ/rz-genpat" "$outabs.S" "$m" "$NPAT" "$outabs.pat$m" "$m"
  "$RZ/rz-bench" "$outabs.rz" "$outabs.rix" "$outabs.pat$m"
done
