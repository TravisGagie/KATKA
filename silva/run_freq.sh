#!/bin/bash
# Credit each genus in proportion to the MEM's occurrences in it, instead of evenly among the genera containing
# it.  Listings with occurrence counts (RZ_COUNTS=1; sr-index s = 0 on the explicit RLCSA) for every 50th read
# pair (200 K per region), digested and undigested, L = 15 and 30; scored by rz_score2.py in three modes:
#   list (even split, as in the paper), freq (proportional to occurrences),
#   norm (proportional to the fraction of the genus's sequences, both strands, that contain the MEM).
# Classification is not timed, so two regions run at a time.  Output: work/results/freq/summary.txt
LS=${LS-"15 30"}; MODES=${MODES-"list freq norm"}
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; R=$W/results/freq
mkdir -p $R
one() {   # ds L region
  local ds=$1 L=$2 r=$3 D=$W/$1 d=$R/$1_L$2/$3; mkdir -p $d
  local M=""; [ $ds = rzdg ] && M="-B $D/bac.map"
  for m in 1 2; do
    [ -s $d/mate_$m.listings ] && continue
    $S/reads.sh $r $m 50 > $d/reads_$m.fq
    RZ_COUNTS=1 $B -L $L -l $M -C $D/bac.csa -S $D/bac.s0.sri -R $D/bac.rix $D/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log
    rm -f $d/reads_$m.fq
  done
  for mode in $MODES; do
    [ -s $d/$mode.csv ] && continue
    python3 $S/rz_score2.py --mode $mode --nseq $D/bac.nseq --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $d/$mode.csv --exp-src $RZ/cliffy-experiments/src --threads 3 > $d/$mode.log 2>&1
  done
}
for ds in rzdg rz; do for L in $LS; do
  [ $ds = rz ] && [ $L = 50 ] && continue
  echo "== $(date '+%T') $ds L=$L"
  one $ds $L V1_V2 & one $ds $L V3_V4 & wait
  one $ds $L V4_V4 & one $ds $L V4_V5 & wait
done; done
python3 - $R <<'PY' | tee $R/summary.txt
import sys, os
R = sys.argv[1]
print("genus level, mean over the four regions: strict / Cliffy's scoring")
for ds in ("rzdg", "rz"):
    for L in (15, 20, 30, 40, 50):
        row = []
        for mode in ("list", "freq", "norm"):
            if not all(os.path.exists(f"{R}/{ds}_L{L}/{r}/{mode}.csv") for r in ("V1_V2", "V3_V4", "V4_V4", "V4_V5")): continue
            st = th = 0
            for r in ("V1_V2", "V3_V4", "V4_V4", "V4_V5"):
                for l in open(f"{R}/{ds}_L{L}/{r}/{mode}.csv"):
                    f = l.strip().split(",")
                    if f[2] == "genus": dd = dict(x.split("=") for x in f[4:]); st += float(dd["strict"]); th += float(dd["theirs"])
            row.append(f"{mode} {100*st/4:.2f} / {100*th/4:.2f}")
        print(f"{'digest' if ds == 'rzdg' else 'undigested':<11} L={L:<3}", "   ".join(row))
PY
