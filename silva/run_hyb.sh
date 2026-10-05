#!/bin/bash
# Hybrid LCA/listing on the digested SILVA: BML with threshold L; a MEM with more than T occurrences is
# answered by its LCA (rz-index grids), the others by listing (sr-index).  Accuracy on the same 200 K pairs per
# region as run_bml.sh (sr s = 0, explicit RLCSA, with aux); speed on the same 2,000 reads, sr with
# s = 0, 4, 16 on RLBWT / explicit / EF RLCSA, with and without the aux grid fixes.
# Log: silva/run_hyb.log; results under work/results/bml_hyb and work/rzdg/hyb_speed.out
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rzdg
LS=${LS-"15 20 30"}; TS=${TS-"100 300 1000"}
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
exec 9>$D/.lock; flock -n 9 || fail "another run is using the digest index"
R=$W/results/bml_hyb; mkdir -p $R
acc_region() {   # $1 = region, $2 = L, $3 = T
  local r=$1 L=$2 T=$3 d=$R/hyb_L$2_T$3/$1
  [ -s $d/done ] && return 0
  mkdir -p $d
  for m in 1 2; do
    awk 'int((NR-1)/4) % 50 == 0' $W/reads/aquatic/$r/${r}_mate_$m.fq > $d/reads_$m.fq
    $B -L $L -l -H $T -B $D/bac.map -C $D/bac.csa -x $D/bac.aux -S $D/bac.s0.sri -R $D/bac.rix $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log \
      || { echo "rz-classify failed: $r L=$L T=$T mate $m"; return 1; }
    rm $d/reads_$m.fq
  done
  python3 $S/rz_score2.py --mode list --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $d/scores.csv --exp-src $RZ/cliffy-experiments/src --threads 2 > $d/classify.log 2>&1 || { echo "scoring failed: $r L=$L T=$T"; return 1; }
  date > $d/done
}
for L in $LS; do for T in $TS; do
  step "accuracy L=$L T=$T"
  pids=""; for r in V1_V2 V3_V4 V4_V4 V4_V5; do acc_region $r $L $T & pids="$pids $!"; done
  for p in $pids; do wait $p || fail "accuracy L=$L T=$T"; done
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do echo "$r $(grep genus $R/hyb_L${L}_T$T/$r/scores.csv)"; done
done; done

step "speed: 2,000 reads per region"
SP=$D/hyb_speed.out; [ -s $SP ] || : > $SP
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  awk 'int((NR-1)/4) % 5000 == 0' $W/reads/aquatic/$r/${r}_mate_1.fq > $D/hyb_reads.fq
  for L in $LS; do for T in $TS; do
    grep -q "^$r L=$L T=$T done" $SP && continue
    O=$D/hyb_out; mkdir -p $O; rm -f $O/*
    for v in rl csa csaef; do
      c=""; [ $v != rl ] && c="-C $D/bac.$v"
      for x in aux noaux; do
        a=""; [ $x = aux ] && a="-x $D/bac.aux"
        for s in 0 4 16; do
          echo -n "$r L=$L T=$T sr s=$s $v $x: " >> $SP
          $B -L $L -l -H $T -B $D/bac.map $c $a -S $D/bac.s$s.sri -R $D/bac.rix $D/bac.rz $D/hyb_reads.fq $O/h_${s}_${v}_$x 2>> $SP || echo "EXIT $?" >> $SP
        done
      done
    done
    bad=0; for f in $O/h_*; do cmp -s $f $O/h_0_rl_aux || { echo "  $r L=$L T=$T $(basename $f): DIFFERS" >> $SP; bad=1; }; done
    echo "$r L=$L T=$T done (answers $([ $bad = 0 ] && echo identical || echo DIFFER))" >> $SP
  done; done
done
rm -rf $D/hyb_out $D/hyb_reads.fq
step "ALL DONE"
