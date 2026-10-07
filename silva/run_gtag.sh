#!/bin/bash
# Grammar-compressed tag array (rz-index/gtag.hpp): builds work/rz/bac.s1.gtag if missing (the genus sequence of the
# runs of work/rz/bac.s1.tag, compressed with RePair: Navarro's irepair, which avoids deep rules, from BigRePair's
# largeb_repair; about 5 minutes and 9 GB on SILVA), checks that classification gives identical answers, and
# times it against the run-length tag array as run_pareto_ft.sh does (COUNTS=1, 2,000 reads per region, undigested,
# L = 30 and 40, RLBWT / explicit / Elias-Fano RLCSA, lookup tables of 0, 10 and 12 bases).  Run with the Claude app
# closed.  Output: work/pareto_speed_gtag.out (read by pareto_ft.py --counts).
# Needs irepair (IREPAIR=path; default ../deps/bigrepair/largeb_repair/irepair; REPAIR_MB=12000) unless the grammar exists.
# GTAG=name: time work/rz/<name> (an existing .gtag file, e.g. from another grammar) instead; output pareto_speed_<name>.out
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work; D=$W/rz; T=$RZ/tools/grammar
BR=${IREPAIR:-$RZ/deps/bigrepair/largeb_repair/irepair}
fail() { echo "FAILED: $*"; exit 1; }
make -C $B rz-classify rz-gtagbuild > /dev/null || fail make
GT=${GTAG:-bac.s1.gtag}
if [ ! -s $D/$GT ]; then
  [ -n "$GTAG" ] && fail "no $D/$GTAG"
  G=$D/tags_genus.int32; [ -s $RZ/tools/tagz/rp_genus.int32.R ] && G=$RZ/tools/tagz/rp_genus.int32
  if [ ! -s $G.R ]; then
    [ -x $T/tokwrite ] || g++ -O3 -std=c++17 $T/tokwrite.cpp -o $T/tokwrite || fail "tokwrite"
    [ -s $G ] || $T/tokwrite $D/bac.tagruns genus $G || fail "tokwrite"
    [ -x $BR ] || fail "no irepair at $BR (git clone https://gitlab.com/manzai/bigrepair.git ../deps/bigrepair && make -C ../deps/bigrepair)"
    /usr/bin/time -f "irepair: %e s, peak %M KB" $BR $G ${REPAIR_MB:-12000} || fail irepair
  fi
  /usr/bin/time -f "rz-gtagbuild: %e s, peak %M KB" $B/rz-gtagbuild $D/bac.s1.tag $G $D/bac.s1.gtag || fail rz-gtagbuild
fi
ls -l $D/bac.s1.tag $D/$GT
export RZ_COUNTS=1; OUT=$W/pareto_speed_gtag.out; [ -n "$GTAG" ] && OUT=$W/pareto_speed_$GTAG.out; : > $OUT
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  $S/reads.sh $r 1 5000 > $D/pg_reads.fq
  for L in 30 40; do for v in rl csa csaef; do
    c=""; [ $v != rl ] && c="-C $D/bac.$v"
    for F in 0 10 12; do
      f=""; [ $F != 0 ] && f="-F $F"
      for kind in tag gtag; do
        X=$D/bac.s1.tag; [ $kind = gtag ] && X=$D/$GT
        echo -n "rz $r L=$L $kind s=1 $v F=$F: " >> $OUT
        $B/rz-classify -L $L -l $c -T $X $f $D/bac.rz $D/pg_reads.fq $D/pg_out_$kind 2>&1 | grep -o "[0-9.]* us/read\|table [0-9.]* MB" | tr '\n' ' ' >> $OUT
        [ $kind = tag ] && echo >> $OUT
      done
      cmp -s $D/pg_out_tag $D/pg_out_gtag && echo "[identical]" >> $OUT || echo "[ANSWERS DIFFER]" >> $OUT
    done
  done; done
  rm -f $D/pg_reads.fq $D/pg_out_*
done
echo "done: $(grep -c identical $OUT) identical, $(grep -c DIFFER $OUT) differ"
awk '/ L=30 / && / csa F=10:/' $OUT
