#!/bin/bash
# SILVA with minimizer digests, k = 4, w = 11 (Cliffy's and SPUMONI 2's default), strand-symmetric and
# byte-mapped (rz-prep -B; the 6 rarest of the 256 minimizers share one byte: 94 of 261 M occurrences).
#  1. rz-index, r-index, sr-index (s = 4, 16), explicit and Elias-Fano RLCSAs over the digested reference
#  2. accuracy: all four aquatic regions, every 10th read pair (1 M), rz-index on the explicit RLCSA,
#     once with Ziv-Merhav phrases (Cliffy's features) and once with MEMs (rz-classify -M),
#     scored with rz_score.py (their method0 = Cliffy's LCA query)
#  3. speed: 2,000 reads of each region, both feature types,, rz-index and sr-index on each backend (answers must be identical)
# Log: silva/run_silva_dg.log;  speed table: silva/work/rzdg/speed.out
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work; D=$W/rzdg
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
T() { /usr/bin/time --format='user= %U system= %S elapsed= %e CPU= %P MemMax= %M' "$@"; }
mkdir -p $D; exec 9>$D/.lock; flock -n 9 || fail "another run_silva_dg.sh is running"; cp -n $W/rz/list.tsv $D/list.tsv
step "indexes over the digested reference"
OUT=$D PREPFLAGS="-B 4,11" NOMV=1 SVALS="4 16" LENGTHS="" THREADS=6 $RZ/runrz.sh || fail "runrz.sh"
[ -s $D/bac.csa ]   || $B/rz-csabuild $D/bac.bwt $D/bac.csa || fail csa
[ -s $D/bac.csaef ] || $B/rz-csabuild $D/bac.bwt $D/bac.csaef ef || fail csaef
cat $D/rz-build.out; ls -l $D/bac.rz $D/bac.aux $D/bac.rix $D/bac.s*.sri $D/bac.csa $D/bac.csaef $D/bac.map

for F in zm mem; do
  MF=""; [ $F = mem ] && MF=-M
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do
    step "accuracy ($F): aquatic $r"
    d=$W/results/rzdg_$F/aquatic/$r; [ -s $d/done ] && { cat $d/output.classification_results.csv; continue; }
    mkdir -p $d
    for m in 1 2; do
      grep -q "^reads=" $d/mate_$m.log 2>/dev/null && continue
      awk 'int((NR-1)/4) % 10 == 0' $W/reads/aquatic/$r/${r}_mate_$m.fq > $d/reads_$m.fq
      T --output=$d/mate_$m.time $B/rz-classify $MF -B $D/bac.map -C $D/bac.csa -x $D/bac.aux $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log \
        || fail "rz-classify $F $r $m"
      rm $d/reads_$m.fq; cat $d/mate_$m.log
    done
    python3 $S/rz_score.py --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
        --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
        --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
        --output $d/output.classification_results.csv --exp-src $RZ/cliffy-experiments/src > $d/classify.log 2>&1 || fail "scoring $F $r"
    cat $d/output.classification_results.csv; date > $d/done
  done
done

step "speed: 2,000 reads per region, both feature types, every index and backend"
SP=$D/speed.out; : > $SP
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  awk 'int((NR-1)/4) % 5000 == 0' $W/reads/aquatic/$r/${r}_mate_1.fq > /tmp/dg_reads.fq      # 2,000 reads
  for F in zm mem; do
    MF=""; [ $F = mem ] && MF=-M
    for v in rl csa csaef; do
      c=""; [ $v != rl ] && c="-C $D/bac.$v"
      echo -n "$r $F rz $v: " | tee -a $SP
      $B/rz-classify $MF -B $D/bac.map $c -x $D/bac.aux $D/bac.rz /tmp/dg_reads.fq /tmp/dg_rz_$v 2>&1 | tee -a $SP
      for s in 4 16; do
        echo -n "$r $F sr s=$s $v: " | tee -a $SP
        $B/rz-classify $MF -B $D/bac.map $c -S $D/bac.s$s.sri -R $D/bac.rix $D/bac.rz /tmp/dg_reads.fq /tmp/dg_s${s}_$v 2>&1 | tee -a $SP
      done
    done
    for f in /tmp/dg_rz_csa /tmp/dg_rz_csaef /tmp/dg_s4_rl /tmp/dg_s4_csa /tmp/dg_s4_csaef /tmp/dg_s16_rl /tmp/dg_s16_csa /tmp/dg_s16_csaef; do
      cmp -s /tmp/dg_rz_rl $f && echo "  $r $F $(basename $f): identical" | tee -a $SP || echo "  $r $F $(basename $f): DIFFERS" | tee -a $SP
    done
    rm -f /tmp/dg_rz_* /tmp/dg_s*
  done
done
step "ALL DONE"
