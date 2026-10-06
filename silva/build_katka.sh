#!/bin/bash
# Builds the SILVA 138.1 reference (one document per genus, in tree order, as in Cliffy's experiments) and the
# KATKA indexes you ask for; by default only those of the default configuration.  No reads are needed.
# Resumable: finished steps are skipped.
#   ./build_katka.sh                       default configuration: undigested text (work/rz), explicit RLCSA,
#                                          unsampled tag array
#   INDEXES="rz rzdg" ./build_katka.sh     also the same over minimizer digests (work/rzdg; k = 4, w = 11)
#   ALL=1 ./build_katka.sh                 every index in the paper (both texts, all of the below)
# What to build, for each text in INDEXES (defaults in brackets):
#   CSA="explicit ef"   RLCSAs: explicit and/or Elias-Fano                                    [explicit]
#   TAGS="1 4 16"       tag arrays with these sampling parameters (1 = unsampled)              [1]
#   SR="0 4 16"         r-index and sr-indexes with these sampling parameters (0 = r-index)    [none]
#   GRIDS=1             the rz-index's LZ77 grids and their auxiliary structures, for LCA
#                       queries without tags (rz-classify without -T or -S, -H/-A)             [0]
#   NSEQ=1              number of sequences per genus (for normalized credit and -A)           [0]
#   VFY=1               verification structures for digested matches (rzdg only)              [0]
# (run_2l.sh builds the two-level parse indexes.)  THREADS (default: all cores) for the LZ77 parses.
# Needs ../setup_deps.sh (sdsl-lite, Big-BWT, cliffy-experiments, SILVA), make -C ../rz-index, and seqtk.
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index; W=$S/work; EXP=$RZ/cliffy-experiments/src
FASTA=$S/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta
TAXTXT=$S/exp1_data/tax_slv_ssu_138.1.txt; TAXTRE=$S/exp1_data/tax_slv_ssu_138.1.tre
ALL=${ALL:-0}; THREADS=${THREADS:-$(nproc)}
if [ $ALL = 1 ]; then INDEXES="rz rzdg"; CSA="explicit ef"; TAGS="1 4 16"; SR="0 4 16"; GRIDS=1; NSEQ=1; VFY=1; fi
INDEXES=${INDEXES:-rz}; CSA=${CSA:-explicit}; TAGS=${TAGS-1}; SR=${SR-}; GRIDS=${GRIDS:-0}; NSEQ=${NSEQ:-0}; VFY=${VFY:-0}
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
  local D=$W/$1 prep="" ; [ $1 = rzdg ] && prep="-B 4,11"
  mkdir -p $D; cd $D
  [ -s list.tsv ] || for i in $(seq 1 $n); do printf 'silva\t%s\n' "$W/ref/dna/doc_${i}_seq.fa"; done > list.tsv
  step "$1: text and BWT"
  [ -s bac.datasets ] || $B/rz-prep -d $prep list.tsv $D/bac > prep.log 2>&1 || fail "rz-prep $1 (see $D/prep.log)"
  if [ ! -s bac.bwt ]; then
    [ $(wc -l < bac.datasets) = 1 ] || fail "expected one dataset in $D/bac.datasets"
    [ -s bac.ds0.S.bwt ] || $RZ/deps/Big-BWT/bigbwt $D/bac.ds0.S > bac.ds0.S.bigbwt.log 2>&1 || fail "Big-BWT $1 (see $D/bac.ds0.S.bigbwt.log)"
    cp bac.ds0.S.bwt bac.bwt
  fi
  [ $(stat -c %s bac.bwt) = $(( $(stat -c %s bac.S) + 1 )) ] || fail "$D/bac.bwt has the wrong length"
  if [ $GRIDS = 1 ]; then
    if [ ! -s bac.rz ] || [ -e bac.rz.nogrids ]; then
      step "$1: LZ77 parses and the rz-index's grids"
      [ -s bac.left ]  || $B/rz-lz77    -t $THREADS -w 32M -b 32M bac.S bac.left  > lz77.log 2>&1 || fail "left parse $1"
      [ -s bac.right ] || $B/rz-lz77 -r -t $THREADS -w 32M -b 32M bac.S bac.right >> lz77.log 2>&1 || fail "right parse $1"
      $B/rz-build bac.S bac.bwt bac.left bac.right bac.tbl bac > rz-build.out 2>&1 || fail "rz-build $1 (see $D/rz-build.out)"
      rm -f bac.rz.nogrids bac.aux
    fi
    [ -s bac.aux ] || $B/rz-aux bac.rz bac.aux > aux-build.out 2>&1 || fail "rz-aux $1"
  elif [ ! -s bac.rz ]; then
    step "$1: RLBWT (no grids)"
    $B/rz-build -G bac.S bac.bwt - - bac.tbl bac > rz-build.out 2>&1 || fail "rz-build -G $1 (see $D/rz-build.out)"
    : > bac.rz.nogrids
  fi
  for c in $CSA; do
    f=bac.csa; [ $c = ef ] && f=bac.csaef
    [ -s $f ] || { step "$1: $c RLCSA"; $B/rz-csabuild bac.bwt $f $c || fail "csa $c $1"; }
  done
  for s in $TAGS; do [ -s bac.s$s.tag ] || { step "$1: tag array, s = $s"; $B/rz-tagbuild bac $s bac.s$s.tag || fail "tags s=$s $1"; }; done
  if [ -n "$SR" ]; then
    [ -s bac.rix ] || { step "$1: r-index"; $B/rz-build -A bac.S bac.bwt - - bac.tbl bac > ri-build.out 2>&1 || fail "r-index $1 (see $D/ri-build.out)"; }
    for s in $SR; do [ -s bac.s$s.sri ] || { step "$1: sr-index, s = $s"; $B/sr-build bac.rix - $s bac.s$s.sri >> sr-build.out 2>&1 || fail "sr-index s=$s $1"; }; done
  fi
  [ $NSEQ = 1 ] && { [ -s bac.nseq ] || for i in $(seq 1 $n); do grep -c '>' $W/ref/rna/doc_${i}_seq.fa; done > bac.nseq; }
  [ $VFY = 1 ] && [ $1 = rzdg ] && { [ -s bac.vfy ] || $B/rz-vfybuild bac.tbl bac.map bac.vfy || fail vfy; }
  ls -l bac.rz bac.csa* bac.s*.tag bac.s*.sri bac.rix 2>/dev/null | awk '{print "  " $5 "\t" $NF}'
  cd $W
}
for x in $INDEXES; do index $x; done
step "done"
