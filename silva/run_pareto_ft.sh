#!/bin/bash
# Timings for the Pareto surface with lookup tables (rz-classify -F t): BML with listing, on the digested
# (work/rzdg, t = 0 or 2 digest symbols) and undigested (work/rz, t = 0, 10 or 12 bases) references; tag arrays
# (s = 1, 4, 16) and sr-indexes (s = 0, 4, 16) on the RLBWT, explicit and Elias-Fano RLCSAs; plus two-level BML
# (undigested, k = 12, s = 4, t = 12).  2,000 reads per region, L = 15, 20, 30, 40 (and 50 for digests); the
# answers must equal those without a table.  Run with the Claude app closed.  Output: work/pareto_speed.out
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work
# COUNTS=1: list each genus with its number of occurrences (RZ_COUNTS), for proportional credit
if [ "${COUNTS-0}" = 1 ]; then export RZ_COUNTS=1; OUT=$W/pareto_speed_counts.out; else OUT=$W/pareto_speed.out; fi
: > $OUT
for ds in rzdg rz; do
  D=$W/$ds
  if [ $ds = rzdg ]; then M="-B $D/bac.map"; FS="0 2"; LS="15 20 30 40 50"; else M=""; FS="0 10 12"; LS="15 20 30 40"; fi
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do
    $S/reads.sh $r 1 5000 > $D/pf_reads.fq
    for L in $LS; do
      for v in rl csa csaef; do
        c=""; [ $v != rl ] && c="-C $D/bac.$v"
        for kind in tag sr; do
          for s in $([ $kind = tag ] && echo "1 4 16" || echo "0 4 16"); do
            if [ $kind = tag ]; then X="-T $D/bac.s$s.tag"; else X="-S $D/bac.s$s.sri -R $D/bac.rix"; fi
            for F in $FS; do
              f=""; [ $F != 0 ] && f="-F $F"
              echo -n "$ds $r L=$L $kind s=$s $v F=$F: " >> $OUT
              $B -L $L -l $M $c $X $f $D/bac.rz $D/pf_reads.fq $D/pf_out_$F 2>&1 | grep -o "[0-9.]* us/read\|table [0-9.]* MB" | tr '\n' ' ' >> $OUT
              if [ $F = 0 ]; then echo >> $OUT; else cmp -s $D/pf_out_0 $D/pf_out_$F && echo "[identical]" >> $OUT || echo "[ANSWERS DIFFER]" >> $OUT; fi
            done
          done
        done
      done
    done
    rm -f $D/pf_reads.fq $D/pf_out_*
  done
done
echo "done: $(grep -c identical $OUT) identical, $(grep -c DIFFER $OUT) differ"
