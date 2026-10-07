#!/bin/bash
# Mini test of Cliffy's installation and execution (temporary): Cliffy on a small subset of SILVA.
#  1. builds Cliffy (silva/build_cliffy.sh) at the commit used in the paper, if it is not built yet;
#  2. takes every MINI_EVERY-th genus (default 30: about 300 genera) of the reference that demo_quick.sh builds,
#     renumbered 1..m, and the quick demo's read pairs (silva/data/sample2000) whose true genus is among them;
#  3. builds Cliffy's index(es) on the subset (MINI_INDEXES, default "full_text minimizers", with the settings of
#     run_silva.sh), classifies both mates of the reads of each region and scores them with Cliffy's own script.
# Run ./demo_quick.sh first (it fetches the dependencies and SILVA and builds the reference).  Needs cmake and a
# C++ compiler for Cliffy, and up to 2 x MINI_TMPSIZE (default 16GB) of temporary disk space.  Accuracy on a subset means little: this only checks that everything runs.  Resumable;
# results in silva/work/mini/ (summary.txt).
set -o pipefail
ROOT=$(cd "$(dirname "$0")" && pwd); S=$ROOT/silva; W=$S/work; M=$W/mini
EVERY=${MINI_EVERY:-30}; TMPSIZE=${MINI_TMPSIZE:-16GB}; IXS=${MINI_INDEXES:-"full_text minimizers"}; REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}
EXP=$ROOT/cliffy-experiments/src; TAXTXT=$S/exp1_data/tax_slv_ssu_138.1.txt
CL=$S/cliffy/build/cliffy; export PFPDOC_BUILD_DIR=$S/cliffy/build/install
mkdir -p $M; exec > >(tee -a $M/mini.log) 2>&1
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
T() { /usr/bin/time --format='elapsed %e s, peak %M KB' "$@"; }
[ -s $W/ref/done ] || fail "no reference: run ./demo_quick.sh first"
[ -s $EXP/classify_silva_readset.py ] || fail "no cliffy-experiments: run ./setup_deps.sh"
for t in cmake git /usr/bin/time python3; do command -v $t > /dev/null || fail "missing $t"; done

if [ ! -x $CL ]; then
  step "building Cliffy (log: silva/build_cliffy.log)"
  [ -d $S/cliffy ] || git clone https://github.com/oma219/cliffy.git $S/cliffy || fail "git clone cliffy"
  git -C $S/cliffy checkout -q 3763529 || fail "checkout"
  $S/build_cliffy.sh; [ -x $CL ] || fail "building Cliffy (see silva/build_cliffy.log)"
fi

if [ ! -s $M/ref/done ]; then
  step "subset: every ${EVERY}th genus"
  rm -rf $M/ref; mkdir -p $M/ref/dna $M/ref/rna
  n=$(cut -d' ' -f1 $W/ref/done); j=0
  for i in $(seq 1 $EVERY $n); do
    j=$((j + 1)); cp $W/ref/dna/doc_${i}_seq.fa $M/ref/dna/doc_${j}_seq.fa
    sed -n "${i}s/^[0-9]* /$j /p" $W/ref/rna/doc_to_traversal.txt >> $M/ref/rna/doc_to_traversal.txt
  done
  cp $W/ref/trav_to_length.txt $M/ref/
  echo "$j genera, $(cat $M/ref/dna/doc_*_seq.fa | grep -c '>') sequences, $(cat $M/ref/dna/doc_*_seq.fa | grep -v '>' | tr -d '\n' | wc -c) bases" | tee $M/ref/done
fi
m=$(ls $M/ref/dna/doc_*_seq.fa | wc -l)   # Cliffy's file list (absolute paths, so written on every run)
for j in $(seq 1 $m); do echo "$M/ref/dna/doc_${j}_seq.fa $j"; done > $M/ref/dna/filelist.txt

for r in $REGIONS; do
  d=$M/reads/$r; [ -s $d/done ] && continue; mkdir -p $d
  cp $S/data/sample2000/${r}_seqtax.txt $d/
  for m in 1 2; do   # keep the reads whose true genus (column 2 of the truth file = N in genus_N_read_...) is in the subset
    python3 - $M/ref/rna/doc_to_traversal.txt $d/${r}_seqtax.txt $S/data/sample2000/${r}_mate_$m.fq.gz > $d/${r}_mate_$m.fq <<'PY' || fail "reads $r"
import sys, gzip
keep = set(l.split(" ", 1)[1].strip() for l in open(sys.argv[1]))
ids = set(l.rstrip("\n").split("\t")[1] for l in open(sys.argv[2]) if l.rstrip("\n").split("\t")[0].strip() in keep)
with gzip.open(sys.argv[3], "rt") as f:
    while True:
        rec = [f.readline() for _ in range(4)]
        if not rec[0]: break
        if rec[0][1:].split("_")[1] in ids: sys.stdout.write("".join(rec))
PY
  done
  echo "$r: $(( $(wc -l < $d/${r}_mate_1.fq) / 4 )) pairs" | tee $d/done
done

for ix in $IXS; do
  d=$M/index/$ix; [ -s $d/done ] && continue
  step "Cliffy index $ix"
  rm -rf $d; mkdir -p $d
  case $ix in
    full_text)  opt="--num-col 7" ;;
    minimizers) opt="--num-col 5 --minimizers --small-window 4 --large-window 11" ;;
    *) fail "unknown index $ix" ;;
  esac
  T -o $d/time.txt $CL build --filelist $M/ref/dna/filelist.txt --output $d/output --revcomp --taxcomp $opt \
     --two-pass $d/temp --tmp-size $TMPSIZE 2> $d/build.log || fail "cliffy build $ix (see $d/build.log)"
  rm -rf $d/temp $d/temp.tmp_*; date > $d/done; cat $d/time.txt; du -sh $d
done

for ix in $IXS; do for r in $REGIONS; do
  d=$M/runs/$ix/$r; [ -s $d/output.classification_results.csv ] && continue; mkdir -p $d
  step "Cliffy $ix on $r"
  case $ix in
    full_text)  opt="--num-col 7" ;;
    minimizers) opt="--num-col 5 --minimizers --small-window 4 --large-window 11" ;;
  esac
  for m in 1 2; do
    T -o $d/mate_$m.time $CL run --ref $M/index/$ix/output --pattern $M/reads/$r/${r}_mate_$m.fq --output $d/mate_$m \
       --taxcomp --ftab $opt 2> $d/mate_$m.log || fail "cliffy run $ix $r mate $m (see $d/mate_$m.log)"
  done
  python3 $EXP/classify_silva_readset.py --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
    --doc-id-to-traversal $M/ref/rna/doc_to_traversal.txt --silva-tax-ranks $TAXTXT \
    --readset-truthset $M/reads/$r/${r}_seqtax.txt --output-dir $d/ --trav-to-length $M/ref/trav_to_length.txt \
    --sample-name output --classify-docprof > $d/classify.log 2>&1 || fail "scoring $ix $r (see $d/classify.log)"
done; done

step "summary"
{ echo "host: $(hostname); $(cat $M/ref/done)"; for r in $REGIONS; do cat $M/reads/$r/done; done
  for ix in $IXS; do echo "--- $ix: index $(du -sh $M/index/$ix | cut -f1), build $(cat $M/index/$ix/time.txt)"
    for r in $REGIONS; do echo "  $r: $(cat $M/runs/$ix/$r/mate_1.time); genus: $(grep -h ',genus,' $M/runs/$ix/$r/output.classification_results.csv | tr '\n' ' ')"; done
  done; } | tee $M/summary.txt
echo "Cliffy ran end to end (accuracy on a subset is not meaningful)."
