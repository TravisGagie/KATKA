#!/bin/bash
# Estimates what a one-level phrase index (minimizer phrases, k = 4, w = 11) would lose: each undigested MEM
# (BML, L = 30, listing with the unsampled tag array) is trimmed to the part such an index would certainly match
# (RZ_TRIM=4,11) or by a fixed number of bases at each end (RZ_TRIMFIX=t), and the genera of the trimmed strings are
# listed and scored.  Every 500th read pair of the V1-V2 and V4 regions.  Output: silva/work/results/trim/
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rz; O=$W/results/trim
mkdir -p $O
for r in V1_V2 V4_V4; do
  for m in 1 2; do awk 'int((NR-1)/4) % 500 == 0' $W/reads/aquatic/$r/${r}_mate_$m.fq > $O/${r}_$m.fq; done
  for v in base trim fix1 fix2 fix4 fix6; do
    case $v in base) E="";; trim) E="RZ_TRIM=4,11";; fix*) E="RZ_TRIMFIX=${v#fix}";; esac
    for m in 1 2; do env $E $B -L 30 -l -C $D/bac.csa -T $D/bac.s1.tag $D/bac.rz $O/${r}_$m.fq $O/${r}_${v}_$m.listings 2> $O/${r}_${v}_$m.log; done
    python3 $S/rz_score2.py --mode list --mate1-listings $O/${r}_${v}_1.listings --mate2-listings $O/${r}_${v}_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $O/${r}_$v.csv --exp-src $RZ/cliffy-experiments/src --threads 2 > /dev/null 2>&1
    echo "$r $v $(grep genus $O/${r}_$v.csv | grep -o 'strict=.*')  $(grep -h '^trim' $O/${r}_${v}_1.log)"
    rm -f $O/${r}_${v}_?.listings
  done
  rm -f $O/${r}_?.fq
done
