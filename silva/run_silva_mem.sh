#!/bin/bash
# Undigested SILVA with MEMs, to complete the 2x2 (digested/undigested x Ziv-Merhav/MEMs).
# Waits for run_silva_dg.sh to finish (shares its lock), then:
#  1. accuracy: four aquatic regions, every 10th read pair, rz-index + explicit RLCSA, rz-classify -M
#  2. speed: the same 2,000 reads per region as run_silva_dg.sh, ZM and MEMs, rz and sr (s = 4, 16)
#     on the RLBWT, explicit RLCSA and Elias-Fano RLCSA, with identity checks
# Log: silva/run_silva_mem.log;  speed table: silva/work/rz/speed.out
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work; D=$W/rz
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
T() { /usr/bin/time --format='user= %U system= %S elapsed= %e CPU= %P MemMax= %M' "$@"; }
step "waiting for run_silva_dg.sh"
exec 9>$W/rzdg/.lock; flock 9
step "accuracy (mem, undigested)"
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  step "accuracy (mem): aquatic $r"
  d=$W/results/rz_mem/aquatic/$r; [ -s $d/done ] && { cat $d/output.classification_results.csv; continue; }
  mkdir -p $d
  for m in 1 2; do
    grep -q "^reads=" $d/mate_$m.log 2>/dev/null && continue
    awk 'int((NR-1)/4) % 10 == 0' $W/reads/aquatic/$r/${r}_mate_$m.fq > $d/reads_$m.fq
    T --output=$d/mate_$m.time $B/rz-classify -M -C $D/bac.csa -x $D/bac.aux $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log \
      || fail "rz-classify mem $r $m"
    rm $d/reads_$m.fq; cat $d/mate_$m.log
  done
  python3 $S/rz_score.py --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $d/output.classification_results.csv --exp-src $RZ/cliffy-experiments/src > $d/classify.log 2>&1 || fail "scoring mem $r"
  cat $d/output.classification_results.csv; date > $d/done
done
step "speed (undigested): 2,000 reads per region"
SP=$D/speed.out; : > $SP
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  awk 'int((NR-1)/4) % 5000 == 0' $W/reads/aquatic/$r/${r}_mate_1.fq > /tmp/ud_reads.fq
  for F in zm mem; do
    MF=""; [ $F = mem ] && MF=-M
    for v in rl csa csaef; do
      c=""; [ $v != rl ] && c="-C $D/bac.$v"
      echo -n "$r $F rz $v: " | tee -a $SP
      $B/rz-classify $MF $c -x $D/bac.aux $D/bac.rz /tmp/ud_reads.fq /tmp/ud_rz_$v 2>&1 | tee -a $SP
      for s in 4 16; do
        echo -n "$r $F sr s=$s $v: " | tee -a $SP
        $B/rz-classify $MF $c -S $D/bac.s$s.sri -R $D/bac.rix $D/bac.rz /tmp/ud_reads.fq /tmp/ud_s${s}_$v 2>&1 | tee -a $SP
      done
    done
    for f in /tmp/ud_rz_csa /tmp/ud_rz_csaef /tmp/ud_s4_rl /tmp/ud_s4_csa /tmp/ud_s4_csaef /tmp/ud_s16_rl /tmp/ud_s16_csa /tmp/ud_s16_csaef; do
      cmp -s /tmp/ud_rz_rl $f && echo "  $r $F $(basename $f): identical" | tee -a $SP || echo "  $r $F $(basename $f): DIFFERS" | tee -a $SP
    done
    rm -f /tmp/ud_rz_* /tmp/ud_s*
  done
done
step "ALL DONE"
