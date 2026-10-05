#!/bin/bash
# Tag-array listing (rz-tagbuild / rz-classify -T): build the undigested tag indexes (s = 1, 4, 16; the
# digested ones exist), then time BML + tags on the usual 2,000 reads per region, both reference sets,
# L in 15 20 30 40 50 75 100, backends RLBWT / explicit / EF RLCSA, s = 1, 4, 16, list and LCA.
# Answers equal the sr-index listings (accuracy: work/results/bml*/list_L*), so only time is measured;
# the answers of all tag configurations are checked to agree.
# Log: silva/run_tag.log; times in work/rzdg/tag_speed.out and work/rz/tag_speed.out
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; TB=$RZ/rz-index/rz-tagbuild; W=$S/work
LS=${LS-"15 20 30 40 50 75 100"}
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
for s in 1 4 16; do
  [ -s $W/rz/bac.s$s.tag ] || { step "undigested tags s=$s"; (cd $W/rz && $TB bac $s bac.s$s.tag) || fail "tagbuild s=$s"; }
done
for ds in rzdg rz; do
  D=$W/$ds; M=""; [ $ds = rzdg ] && M="-B $D/bac.map"
  exec 9>$D/.lock; flock -n 9 || fail "another run is using $D"
  step "speed: $ds"
  SP=$D/tag_speed.out; [ -s $SP ] || : > $SP
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do
    awk 'int((NR-1)/4) % 5000 == 0' $W/reads/aquatic/$r/${r}_mate_1.fq > $D/tag_reads.fq
    for L in $LS; do
      grep -q "^$r L=$L done" $SP && continue
      O=$D/tag_out; mkdir -p $O; rm -f $O/*
      for v in rl csa csaef; do
        c=""; [ $v != rl ] && c="-C $D/bac.$v"
        for s in 1 4 16; do
          for mode in list lca; do
            l=""; [ $mode = list ] && l="-l"
            echo -n "$r L=$L tag s=$s $v $mode: " >> $SP
            $B -L $L $l -T $D/bac.s$s.tag $M $c $D/bac.rz $D/tag_reads.fq $O/${mode}_${s}_$v 2>> $SP || echo "EXIT $?" >> $SP
          done
        done
      done
      bad=0
      for f in $O/list_*; do cmp -s $f $O/list_1_rl || { echo "  $r L=$L $(basename $f): DIFFERS" >> $SP; bad=1; }; done
      for f in $O/lca_*; do cmp -s $f $O/lca_1_rl || { echo "  $r L=$L $(basename $f): DIFFERS" >> $SP; bad=1; }; done
      echo "$r L=$L done (answers $([ $bad = 0 ] && echo identical || echo DIFFER))" >> $SP
    done
  done
  rm -rf $D/tag_out $D/tag_reads.fq
done
step "ALL DONE"
