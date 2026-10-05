#!/usr/bin/env bash
# Build everything from the FASTA files and benchmark (rz-index, r-index, sr-index; RLFM-index and
# move-structure backends), optionally over strand-symmetric minimizer digests.
#     DIGEST=3,11 OUT=$HOME/rz/run30m11 nohup ./runrz.sh >> run30m11.log 2>&1 &
# Environment: DIGEST (k,w; empty = no digestion), NSP (first NSP species; 0 = all), OUT, DATA, GLOB,
# THREADS, WINDOW, BWTPAR, LENGTHS (pattern lengths, in symbols of the indexed text), SVALS, NPAT.
# Resumable: finished steps are skipped.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$PWD
RZ=$ROOT/rz-index
BIGBWT=$ROOT/deps/Big-BWT/bigbwt
PFPMERGE=$ROOT/pfp-merge
DIGEST=${DIGEST:-}
NSP=${NSP:-0}
OUT=${OUT:-$ROOT/run$( [ -n "$DIGEST" ] && echo "m${DIGEST/,/_}" )}
DATA=${DATA:-$ROOT/data/30Bacteria}
GLOB=${GLOB:-*.fa.gz}
THREADS=${THREADS:-6}
WINDOW=${WINDOW:-32M}
BWTPAR=${BWTPAR:-4}
LENGTHS=${LENGTHS-"4 5 6 8 10 12 15 20"}
SVALS=${SVALS-"0 4 8 16 32 64"}
NPAT=${NPAT:-1000}
PREPFLAGS=${PREPFLAGS--s}   # rz-prep flags besides -d (default -s: split FASTA files into genomes by sample ID)
MT="python3 $RZ/test/memtime.py"
mkdir -p "$OUT"; cd "$OUT"
step() { echo; echo "== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }

step "list (digest: ${DIGEST:-none})"
if [ ! -s list.tsv ]; then
  files=$(ls "$DATA"/$GLOB); [ "$NSP" -gt 0 ] && files=$(echo "$files" | head -n "$NSP")
  for f in $files; do printf '%s\t%s\n' "$(basename "$f" | cut -d. -f1 | sed 's/_g[0-9]*$//')" "$f"; done > list.tsv
fi
wc -l < list.tsv

step "S and datasets"
if [ ! -s bac.datasets ]; then
  $MT "$RZ/rz-prep" -d ${PREPFLAGS--s} ${DIGEST:+-M $DIGEST -V} list.tsv "$OUT/bac" || fail "rz-prep"
fi
mapfile -t ds < bac.datasets

step "Big-BWT per dataset"
todo=(); for f in "${ds[@]}"; do [ -s "$f.bwt" ] || todo+=("$f"); done
if [ ${#todo[@]} -gt 0 ]; then
  printf '%s\n' "${todo[@]}" | xargs -P "$BWTPAR" -I{} sh -c "$MT $BIGBWT {} > {}.bigbwt.log 2>&1 || { echo 'Big-BWT failed on {}'; exit 255; }" || fail "Big-BWT"
fi

step "pfp-merge"
if [ ! -s bac.bwt ]; then
  if [ "${#ds[@]}" -ge 2 ]; then (cd "$PFPMERGE" && $MT ./PFPmerge.py -p 4 -o "$OUT/bac.bwt" "${ds[@]}") || fail "pfp-merge"
  else cp "${ds[0]}.bwt" bac.bwt; fi
fi
# the merged BWT must have one row per character of S plus one terminator per dataset
want=$(( $(stat -c %s bac.S) + ${#ds[@]} )); got=$(stat -c %s bac.bwt 2>/dev/null || echo 0)
[ "$got" = "$want" ] || fail "merged BWT has $got bytes, expected $want (see bac.ds0.S.mrg.log)"

step "parses"
[ -s bac.left ]  || $MT "$RZ/rz-lz77"    -t "$THREADS" -w "$WINDOW" -b "$WINDOW" bac.S bac.left  || fail "left parse"
[ -s bac.right ] || $MT "$RZ/rz-lz77" -r -t "$THREADS" -w "$WINDOW" -b "$WINDOW" bac.S bac.right || fail "right parse"

step "rz-index";       [ -s bac.rz ]  || { $MT "$RZ/rz-build" bac.S bac.bwt bac.left bac.right bac.tbl bac | tee rz-build.out; [ -s bac.rz ] || fail "rz-build"; }
step "r-index";        [ -s bac.rix ] || { $MT "$RZ/rz-build" -A bac.S bac.bwt - - bac.tbl bac | tee ri-build.out; [ -s bac.rix ] || fail "rz-build -A"; }
step "move structure"; [ -n "${NOMV:-}" ] || [ -s bac.mv ] || $MT "$RZ/rz-mvbuild" bac.bwt bac.mv | tee mv-build.out
step "grid fixes";     [ -s bac.aux ] || $MT "$RZ/rz-aux" bac.rz bac.aux | tee aux-build.out
step "sr-indexes"
need=""; for s in $SVALS; do [ -s bac.s$s.sri ] || need="$need,$s"; done; need=${need#,}
if [ -n "$need" ]; then
  if [ "${need/,/}" = "$need" ]; then $MT "$RZ/sr-build" bac.rix - "$need" bac.s$need.sri | tee -a sr-build.out
  else $MT "$RZ/sr-build" bac.rix - "$need" bac | tee -a sr-build.out; fi
fi

step "benchmark (one index per process, sequentially)"
OUTF=bench.out
: > $OUTF
for m in $LENGTHS; do
  [ -s pat$m ] || "$RZ/rz-genpat" bac.S "$m" "$NPAT" pat$m "$m"
  $MT "$RZ/rz-bench" -o ans$m.rz                          bac.rz pat$m | tee -a $OUTF
  $MT "$RZ/rz-bench" -o ans$m.rzx   -x bac.aux -T 16        bac.rz pat$m | tee -a $OUTF
  $MT "$RZ/rz-bench" -o ans$m.rzmvx -x bac.aux -T 16 -m bac.mv bac.rz pat$m | tee -a $OUTF
  $MT "$RZ/rz-bench" -o ans$m.rix                         bac.rix pat$m | tee -a $OUTF
  for s in $SVALS; do
    $MT "$RZ/rz-bench" -o ans$m.s$s     -m bac.mv  bac.s$s.sri pat$m | tee -a $OUTF
    $MT "$RZ/rz-bench" -o ans$m.s${s}rl -b bac.rix bac.s$s.sri pat$m | tee -a $OUTF
  done
  for a in rzx rzmvx rix $(for s in $SVALS; do echo s$s s${s}rl; done); do
    "$RZ/rz-bench" -c ans$m.rz ans$m.$a | sed "s/^/m=$m rz vs $a: /" | tee -a $OUTF
  done
done
step "ALL DONE"
