#!/bin/bash
# Builds the indexes used by run_bml*.sh, run_tag.sh, run_hyb.sh and run_vfy_acc.sh that run_silva.sh and
# run_silva_dg.sh do not: for both the undigested (work/rz) and digested (work/rzdg) references,
#   explicit and Elias-Fano RLCSAs, sr-indexes with s = 0 (plain r-index), 4 and 16,
#   tag arrays with s = 1 (unsampled), 4 and 16, and the number of sequences per genus (for run_hyb.sh);
# and, for the digested reference only, the verification structures (for run_vfy_acc.sh).
# Run after run_silva.sh and run_silva_dg.sh.  Resumable: existing files are kept.
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work
fail() { echo "FAILED: $*"; exit 1; }
n=$(cut -d' ' -f1 $W/ref/done)
for ds in rz rzdg; do
  D=$W/$ds; [ -s $D/bac.rz ] || fail "no index in $D (run run_silva.sh / run_silva_dg.sh first)"
  echo "== $D"
  [ -s $D/bac.csa ]   || $B/rz-csabuild $D/bac.bwt $D/bac.csa explicit || fail "csa $ds"
  [ -s $D/bac.csaef ] || $B/rz-csabuild $D/bac.bwt $D/bac.csaef ef     || fail "csaef $ds"
  for s in 0 4 16; do [ -s $D/bac.s$s.sri ] || $B/sr-build $D/bac.rix - $s $D/bac.s$s.sri || fail "sr s=$s $ds"; done
  for s in 1 4 16; do [ -s $D/bac.s$s.tag ] || (cd $D && $B/rz-tagbuild bac $s bac.s$s.tag) || fail "tags s=$s $ds"; done
  [ -s $D/bac.nseq ] || for i in $(seq 1 $n); do grep -c '>' $W/ref/rna/doc_${i}_seq.fa; done > $D/bac.nseq
done
[ -s $W/rzdg/bac.vfy ] || $B/rz-vfybuild $W/rzdg/bac.tbl $W/rzdg/bac.map $W/rzdg/bac.vfy || fail vfy
echo "ALL DONE"
