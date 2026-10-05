#!/bin/bash
# Ideal KeBaB before undigested BML: each read is first cut into pseudo-MEMs with exact k-mer membership
# (RZ_KEBAB=k; neither timed nor counted as steps), then BML lists the genera of each MEM with the unsampled
# tag array on the explicit RLCSA.  Same 2,000 reads per region as the speed runs; the answers must be
# identical to those without KeBaB.  Output: silva/work/rz/kebab_speed.out
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rz
OUT=$D/kebab_speed.out; : > $OUT
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  awk 'int((NR-1)/4) % 5000 == 0' $W/reads/aquatic/$r/${r}_mate_1.fq > $D/kebab_reads.fq
  for c in "30 0" "30 20" "30 25" "15 0" "15 12"; do
    set -- $c
    echo -n "$r L=$1 k=$2: " >> $OUT
    RZ_KEBAB=$2 $B -L $1 -l -C $D/bac.csa -T $D/bac.s1.tag $D/bac.rz $D/kebab_reads.fq $D/kebab_out_${1}_$2 2>&1 | tr '\n' ' ' >> $OUT; echo >> $OUT
  done
  cmp -s $D/kebab_out_30_0 $D/kebab_out_30_20 && cmp -s $D/kebab_out_30_0 $D/kebab_out_30_25 && cmp -s $D/kebab_out_15_0 $D/kebab_out_15_12 && echo "  $r: answers identical" >> $OUT || echo "  $r: answers DIFFER" >> $OUT
done
rm -f $D/kebab_reads.fq $D/kebab_out_*
cat $OUT
