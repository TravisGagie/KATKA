#!/bin/bash
# Lookup table of t-mer intervals (rz-classify -F t) for undigested BML with tag-array listing:
# plain BML on the RLBWT, explicit and Elias-Fano RLCSAs with t = 0 (no table), 8, 10, 12, and two-level BML
# (-P, explicit RLCSA) with t = 10 and 12; 2,000 reads per region, L in LS; the answers must be identical.
# Output: silva/work/rz/ftab_speed.out
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rz
LS=${LS-"15 20 30 40"}
OUT=$D/ftab_speed.out; : > $OUT
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  $S/reads.sh $r 1 5000 > $D/ft_reads.fq
  for L in $LS; do
    $B -L $L -l -C $D/bac.csa -T $D/bac.s1.tag $D/bac.rz $D/ft_reads.fq $D/ft_base 2>/dev/null
    while read name opts; do
      echo -n "$r L=$L $name: " >> $OUT
      $B -L $L -l -T $D/bac.s1.tag $opts $D/bac.rz $D/ft_reads.fq $D/ft_out 2>&1 | grep -o "steps/read=[0-9.]*\|[0-9.]* us/read\|lookups: [0-9.]*\|lookup table.*" | tr '\n' ' ' >> $OUT
      cmp -s $D/ft_base $D/ft_out && echo "[identical]" >> $OUT || echo "[ANSWERS DIFFER]" >> $OUT
    done <<LIST
rl
rl-F10      -F 10
csa         -C $D/bac.csa
csa-F8      -C $D/bac.csa -F 8
csa-F10     -C $D/bac.csa -F 10
csa-F12     -C $D/bac.csa -F 12
csaef       -C $D/bac.csaef
csaef-F10   -C $D/bac.csaef -F 10
LIST
  done
done
rm -f $D/ft_reads.fq $D/ft_base $D/ft_out
cat $OUT
