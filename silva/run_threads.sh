#!/bin/bash
# Throughput with rz-classify -j p (p = 1, 2, 4, ..., JS): undigested (tag array s = 1, explicit RLCSA, -F 12) and
# digested (tag array s = 1, explicit RLCSA, -F 2) BML with listing, L = 20, on 20,000 reads per region.  Checks that
# every thread count gives the same output as one thread.  Run with the Claude app closed.  Output: work/threads.out
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work
JS=${JS:-"1 2 4 6 8 12"}; OUT=$W/threads.out; : > $OUT
for ds in rz rzdg; do
  D=$W/$ds; if [ $ds = rzdg ]; then M="-B $D/bac.map"; F=2; else M=""; F=12; fi
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do
    $S/reads.sh $r 1 500 > $D/th_reads.fq
    for j in $JS; do
      echo -n "$ds $r j=$j: " >> $OUT
      $B -j $j -L 20 -l $M -C $D/bac.csa -T $D/bac.s1.tag -F $F $D/bac.rz $D/th_reads.fq $D/th_out_$j 2>&1 \
        | grep -o "[0-9.]* us/read\|wall [0-9.]* s, [0-9]* reads/s" | tr '\n' ' ' >> $OUT
      if [ $j = 1 ]; then echo >> $OUT; else cmp -s $D/th_out_1 $D/th_out_$j && echo "[identical]" >> $OUT || echo "[ANSWERS DIFFER]" >> $OUT; fi
    done
    rm -f $D/th_reads.fq $D/th_out_*
  done
done
echo "done: $(grep -c identical $OUT) identical, $(grep -c DIFFER $OUT) differ"
