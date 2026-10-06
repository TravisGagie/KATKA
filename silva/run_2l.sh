#!/bin/bash
# Two-level indexing with closed syncmers (rz-parsebuild, rz-classify -P) on the undigested SILVA:
#  1. builds the parse indexes for each (k, s) in KS (about 9 GB of memory and a few minutes each);
#  2. times BML + tag-array listing on the explicit RLCSA, with and without the parse index, on the usual
#     2,000 reads per region, for L in LS; the answers must be identical.
# Output: silva/work/rz/2l_speed.out
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work; D=$W/rz
KS=${KS-"12,4 16,5 20,6"}; LS=${LS-"15 20 30 40"}
OUT=$D/2l_speed.out; : > $OUT
for ks in $KS; do
  k=${ks%,*}; s=${ks#*,}
  [ -s $D/bac.k${k}s$s.pix ] || { echo "== building k=$k s=$s"; /usr/bin/time -f "build: %e s, %M KB" $B/rz-parsebuild $D/bac $k $s $D/bac.k${k}s$s.pix 2>&1 | tee -a $OUT; }
  ls -l $D/bac.k${k}s$s.pix $D/bac.k${k}s$s.pix.B | awk '{s+=$5} END {print "files: " s/1e9 " GB"}' | tee -a $OUT
done
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  $S/reads.sh $r 1 5000 > $D/2l_reads.fq
  for L in $LS; do
    echo -n "$r L=$L base: " >> $OUT
    $B/rz-classify -L $L -l -C $D/bac.csa -T $D/bac.s1.tag $D/bac.rz $D/2l_reads.fq $D/2l_out_base 2>&1 | tail -1 >> $OUT
    for ks in $KS; do
      k=${ks%,*}; s=${ks#*,}
      echo -n "$r L=$L k=$k s=$s: " >> $OUT
      $B/rz-classify -L $L -l -C $D/bac.csa -T $D/bac.s1.tag -P $D/bac.k${k}s$s.pix $D/bac.rz $D/2l_reads.fq $D/2l_out_2l 2>&1 | tail -2 | tr '\n' ' ' >> $OUT
      cmp -s $D/2l_out_base $D/2l_out_2l && echo " [identical]" >> $OUT || echo " [ANSWERS DIFFER]" >> $OUT
    done
  done
done
rm -f $D/2l_reads.fq $D/2l_out_*
grep -o "^.*: \|[0-9.]* us/read\|two-level.*per read\|\[.*\]\|build.*\|files.*" $OUT | paste -sd' ' | sed 's/ \(V[0-9]_V[0-9] L=\)/\n\1/g; s/ \(files\)/\n\1/g; s/ \(build\)/\n\1/g'
