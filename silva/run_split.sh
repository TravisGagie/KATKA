#!/bin/bash
# Classify a read file with P rz-classify processes, each on a contiguous slice of the reads, and concatenate
# their outputs in order, so the result is identical to a single run.  Usage:
#   run_split.sh P [rz-classify options] index.rz reads.fq out.listings
# e.g. run_split.sh 16 -L 20 -l -C work/rz/bac.csa -T work/rz/bac.s1.tag -F 12 work/rz/bac.rz reads.fq out.txt
# Each process loads its own copy of the index, so memory grows P-fold (use threads, rz-classify -j, to share
# one copy).  Prints each slice's statistics, the wall-clock time and the throughput in reads per second.
set -e
S=$(cd "$(dirname "$0")" && pwd); B=${RZCLASSIFY:-$(dirname "$S")/rz-index/rz-classify}
P=$1; shift
[ "$P" -ge 1 ] 2>/dev/null && [ $# -ge 3 ] || { echo "usage: run_split.sh P [rz-classify options] index.rz reads.fq out" >&2; exit 1; }
args=("$@"); n=${#args[@]}
OUT=${args[n-1]}; FQ=${args[n-2]}; opts=("${args[@]:0:n-2}")
T=$(mktemp -d "${TMPDIR:-/tmp}/rzsplit.XXXXXX"); trap 'rm -rf "$T"' EXIT
nr=$(( $(wc -l < "$FQ") / 4 )); per=$(( (nr + P - 1) / P ))
awk -v per=$per -v T="$T" '{ print > sprintf("%s/in.%04d.fq", T, int((NR-1)/(4*per))) }' "$FQ"
t0=$(date +%s.%N); pids=()
for f in "$T"/in.*.fq; do
  i=${f##*/in.}; i=${i%.fq}
  "$B" "${opts[@]}" "$f" "$T/out.$i" 2> "$T/err.$i" & pids+=($!)
done
rc=0; for p in "${pids[@]}"; do wait $p || rc=$?; done
t1=$(date +%s.%N)
cat "$T"/out.* > "$OUT"
for e in "$T"/err.*; do echo "slice ${e##*.}: $(tail -n 3 "$e" | tr '\n' ' ')"; done
awk -v n=$nr -v a=$t0 -v b=$t1 -v p=$P 'BEGIN { w = b - a; printf "P=%d reads=%d wall=%.2f s throughput=%.0f reads/s (%.1f us/read wall)\n", p, n, w, n / w, 1e6 * w / n }'
exit $rc
