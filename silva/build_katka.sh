#!/bin/bash
# Builds the SILVA 138.1 reference (one document per genus, in tree order, as in Cliffy's experiments) and KATKA's
# indexes over it.  No reads are needed.  Resumable: finished steps are skipped.
#   ./build_katka.sh          the default configuration only: undigested index (work/rz) with an explicit RLCSA and
#                             an unsampled tag array (about 1.6 GB in memory when classifying)
#   ALL=1 ./build_katka.sh    every index in the paper: also the Elias-Fano RLCSA, sr-indexes (s = 0, 4, 16), sampled
#                             tag arrays (s = 4, 16), the rz-index's grid fixes and the number of sequences per genus,
#                             and the same over minimizer digests (work/rzdg; k = 4, w = 11, byte-mapped), with the
#                             verification structures.  (run_2l.sh builds the two-level parse indexes.)
# Needs ../setup_deps.sh (sdsl-lite, Big-BWT, cliffy-experiments, SILVA), make -C ../rz-index, and seqtk.
# THREADS (default: all cores) for the parses.
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work; EXP=$RZ/cliffy-experiments/src
FASTA=$S/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta
TAXTXT=$S/exp1_data/tax_slv_ssu_138.1.txt; TAXTRE=$S/exp1_data/tax_slv_ssu_138.1.tre
ALL=${ALL:-0}; THREADS=${THREADS:-$(nproc)}
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
command -v python3 >/dev/null || fail "missing python3"
[ -s $W/ref/done ] || command -v seqtk >/dev/null || fail "missing seqtk"
for f in $FASTA $TAXTXT $TAXTRE $EXP/write_silva_genera.py $RZ/deps/Big-BWT/bigbwt $B/rz-prep $B/rz-classify; do
  [ -e $f ] || fail "missing $f (run ../setup_deps.sh and make -C ../rz-index)"; done
mkdir -p $W; cd $W

step "reference: one file per genus, in SILVA tree order"
if [ ! -s ref/done ]; then
  rm -rf ref; mkdir -p ref/rna ref/dna
  seqtk seq -U $FASTA > ref/silva_database.fa || fail seqtk
  python3 $EXP/write_silva_genera.py -i ref/silva_database.fa -o $W/ref/rna/ --tree $TAXTRE --tree-rank $TAXTXT -n 10000 || fail write_silva_genera
  n=$(ls ref/rna/doc_*_seq.fa | wc -l); : > ref/dna/filelist.txt
  for i in $(seq 1 $n); do        # as in their rule: seqtk seq -r (reverse complement; also turns U into A)
    seqtk seq -r ref/rna/doc_${i}_seq.fa > ref/dna/doc_${i}_seq.fa || fail "seqtk -r $i"
    echo "$W/ref/dna/doc_${i}_seq.fa $i" >> ref/dna/filelist.txt
  done
  python3 - <<'PY' || fail trav_to_length
out = open("ref/trav_to_length.txt", "w")     # genus traversal, then total reference length of that genus
for line in open("ref/rna/doc_to_traversal.txt"):
    doc, trav = line.split(None, 1)
    L = sum(len(l.strip()) for l in open(f"ref/rna/doc_{doc}_seq.fa") if not l.startswith(">"))
    out.write(f"{trav.strip()} {L}\n")
PY
  rm ref/silva_database.fa; echo "$n genera" > ref/done
fi
cat ref/done; n=$(cut -d' ' -f1 ref/done)

index() {   # $1 = rz (undigested) or rzdg (digested)
  local D=$W/$1 prep="" sv=""; [ $1 = rzdg ] && prep="-B 4,11"; [ $ALL = 1 ] && sv="0 4 16"
  step "$1: text, BWT, rz-index, r-index${sv:+, sr-indexes $sv}"
  mkdir -p $D
  [ -s $D/list.tsv ] || for i in $(seq 1 $n); do printf 'silva\t%s\n' "$W/ref/dna/doc_${i}_seq.fa"; done > $D/list.tsv
  OUT=$D PREPFLAGS="$prep" NOMV=1 SVALS="$sv" LENGTHS="" THREADS=$THREADS $RZ/runrz.sh > $D/build.log 2>&1 || { tail $D/build.log; fail "runrz.sh $1 (see $D/build.log)"; }
  step "$1: RLCSA and tag arrays"
  [ -s $D/bac.csa ] || $B/rz-csabuild $D/bac.bwt $D/bac.csa explicit || fail "csa $1"
  for s in 1 $([ $ALL = 1 ] && echo 4 16); do [ -s $D/bac.s$s.tag ] || (cd $D && $B/rz-tagbuild bac $s bac.s$s.tag) || fail "tags s=$s $1"; done
  if [ $ALL = 1 ]; then
    [ -s $D/bac.csaef ] || $B/rz-csabuild $D/bac.bwt $D/bac.csaef ef || fail "csaef $1"
    [ -s $D/bac.nseq ] || for i in $(seq 1 $n); do grep -c '>' $W/ref/rna/doc_${i}_seq.fa; done > $D/bac.nseq
  fi
  ls -l $D/bac.rz $D/bac.csa $D/bac.s*.tag $D/bac.s*.sri $D/bac.csaef 2>/dev/null | awk '{print "  " $5 "\t" $NF}'
}
index rz
if [ $ALL = 1 ]; then
  index rzdg
  [ -s $W/rzdg/bac.vfy ] || $B/rz-vfybuild $W/rzdg/bac.tbl $W/rzdg/bac.map $W/rzdg/bac.vfy || fail vfy
fi
step "done"
