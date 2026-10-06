#!/bin/bash
# Boyer-Moore-Li on the undigested SILVA (plain DNA; same design as run_bml.sh): MEMs of at least L bases, then
# either the LCA of their leftmost/rightmost genera, or (list) the genera of all their occurrences.
#  1. accuracy: L in LS, LCA (rz-index) and list (sr-index), every 50th read pair (200 K pairs per region),
#     four regions in parallel; scored by rz_score2.py (unclassified reads count as wrong)
#  2. speed: 2,000 reads per region (as before), every combination: rz on RLBWT / explicit / EF RLCSA,
#     with and without the aux grid fixes; sr with s = 0, 4, 16 on the same backends, LCA and list;
#     answers must agree across backends
# Log: ~/rz/silva/run_bml_ud.log; results under work/results/bml_ud and work/rz/bml_speed.out
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rz
LS=${LS-"15 20 30 40 50 75 100"}
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
exec 9>$D/.lock; flock -n 9 || fail "another run is using the digest index"
R=$W/results/bml_ud; mkdir -p $R
acc_region() {   # $1 = region, $2 = L, $3 = mode
  local r=$1 L=$2 M=$3 d=$R/$3_L$2/$1
  [ -s $d/done ] && return 0
  mkdir -p $d
  for m in 1 2; do
    $S/reads.sh $r $m 50 > $d/reads_$m.fq
    if [ $M = lca ]; then $B -L $L -C $D/bac.csa -x $D/bac.aux $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log
    else $B -L $L -l -C $D/bac.csa -S $D/bac.s0.sri -R $D/bac.rix $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log; fi \
      || { echo "rz-classify failed: $r L=$L $M mate $m"; return 1; }
    rm $d/reads_$m.fq
  done
  python3 $S/rz_score2.py --mode $M --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $d/scores.csv --exp-src $RZ/cliffy-experiments/src --threads 2 > $d/classify.log 2>&1 || { echo "scoring failed: $r L=$L $M"; return 1; }
  rm $d/mate_?.listings; date > $d/done
}
for L in $LS; do for M in lca list; do
  step "accuracy L=$L $M"
  pids=""; for r in V1_V2 V3_V4 V4_V4 V4_V5; do acc_region $r $L $M & pids="$pids $!"; done
  for p in $pids; do wait $p || fail "accuracy L=$L $M"; done
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do echo "$r $(grep genus $R/${M}_L$L/$r/scores.csv)"; cat $R/${M}_L$L/$r/mate_1.log; done
done; done

step "speed: 2,000 reads per region, every combination"
SP=$D/bml_speed.out; [ -s $SP ] || : > $SP
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  $S/reads.sh $r 1 5000 > $D/bml_reads.fq
  for L in $LS; do
    grep -q "^$r L=$L done" $SP && continue
    O=$D/bml_out; mkdir -p $O; rm -f $O/*
    for v in rl csa csaef; do
      c=""; [ $v != rl ] && c="-C $D/bac.$v"
      for x in aux noaux; do
        a=""; [ $x = aux ] && a="-x $D/bac.aux"
        echo -n "$r L=$L rz $v $x: " >> $SP
        $B -L $L $c $a $D/bac.rz $D/bml_reads.fq $O/rz_${v}_$x 2>> $SP || echo "EXIT $?" >> $SP
      done
      for s in 0 4 16; do
        echo -n "$r L=$L sr s=$s $v lca: " >> $SP
        $B -L $L $c -S $D/bac.s$s.sri -R $D/bac.rix $D/bac.rz $D/bml_reads.fq $O/srlca_${s}_$v 2>> $SP || echo "EXIT $?" >> $SP
        echo -n "$r L=$L sr s=$s $v list: " >> $SP
        $B -L $L -l $c -S $D/bac.s$s.sri -R $D/bac.rix $D/bac.rz $D/bml_reads.fq $O/srlist_${s}_$v 2>> $SP || echo "EXIT $?" >> $SP
      done
    done
    bad=0
    for f in $O/rz_* $O/srlca_*; do cmp -s $f $O/rz_rl_aux || { echo "  $r L=$L $(basename $f): DIFFERS" >> $SP; bad=1; }; done
    for f in $O/srlist_*; do cmp -s $f $O/srlist_0_rl || { echo "  $r L=$L $(basename $f): DIFFERS" >> $SP; bad=1; }; done
    echo "$r L=$L done (answers $([ $bad = 0 ] && echo identical || echo DIFFER))" >> $SP
  done
done
rm -rf $D/bml_out $D/bml_reads.fq
step "ALL DONE"
