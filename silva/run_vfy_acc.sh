#!/bin/bash
# accuracy of verified listing (rz-classify -V) at L, digested SILVA, same 200 K pairs per region as run_bml.sh
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rzdg; L=${L:-30}
R=$W/results/bml_vfy; mkdir -p $R
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  d=$R/vlist_L$L/$r; [ -s $d/done ] && continue; mkdir -p $d
  for m in 1 2; do
    awk 'int((NR-1)/4) % 50 == 0' $W/reads/aquatic/$r/${r}_mate_$m.fq > $d/reads_$m.fq
    $B -L $L -l -B $D/bac.map -C $D/bac.csa -S $D/bac.s0.sri -R $D/bac.rix -V $D/bac.vfy $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log || exit 1
    rm $d/reads_$m.fq
  done
  python3 $S/rz_score2.py --mode list --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $d/scores.csv --exp-src $RZ/cliffy-experiments/src --threads 2 > $d/classify.log 2>&1 || exit 1
  date > $d/done; echo "$r $(grep genus $d/scores.csv)"
done
echo ALL DONE
