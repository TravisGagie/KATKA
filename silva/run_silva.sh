#!/bin/bash
# Reproduces the Cliffy paper's 16S experiment (Ahmed, Boucher, Langmead, Genome Res 2025) on this
# desktop, using their scripts (./cliffy-experiments), MicrobeMixer (./MicrobeMixer) and
# Cliffy (silva/cliffy/build/cliffy) with the settings in their Snakemake rules.
# Resumable: every step is skipped if its output exists.  Log: silva/run_silva.log
#
# Environment (defaults): BIOMES="aquatic"  REGIONS="V1_V2 V3_V4 V4_V4 V4_V5"  NREADS=10000000
#                         INDEXES=""  TMPSIZE=10GB
# What differs from their pipeline, and why:
#   - abundance profiles come from MicrobeMixer/data/<biome>.tsv (their MGnify snapshot) instead of a
#     live MGnify query, which would give different genera today; "Candidatus_" is respelled "Candidatus "
#     to match SILVA 138.1;
#   - one query trial instead of two;
#   - trav_to_length.txt (not produced by their repository) = total reference length per genus;
#   - Kraken 2 / Bracken not yet included.
#   - Cliffy indexes are NOT built by default (INDEXES=""): Cliffy's two-pass build first writes the
#     uncompressed document-array profiles to temporary files, 1,172 GB for the minimizer index and
#     2,265 GB for full text (their Table 1), far more disk than this desktop has.
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S")
EXP=$RZ/cliffy-experiments/src; MM=$RZ/MicrobeMixer; CL=$S/cliffy/build/cliffy
export PFPDOC_BUILD_DIR=$S/cliffy/build/install   # Cliffy finds its helper programs (newscan.x, pfbwt.x, ...) in $PFPDOC_BUILD_DIR/bin
W=$S/work; mkdir -p $W; cd $W || exit 1
BIOMES=${BIOMES:-aquatic}; REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; NREADS=${NREADS:-10000000}
INDEXES=${INDEXES-""}   # Cliffy indexes: off by default, see the note above
TMPSIZE=${TMPSIZE:-10GB}  # Cliffy temporary-file size, as in their config
FASTA=$S/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta
TAXTXT=$S/exp1_data/tax_slv_ssu_138.1.txt; TAXTRE=$S/exp1_data/tax_slv_ssu_138.1.tre
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
T() { /usr/bin/time --format='user= %U system= %S elapsed= %e CPU= %P MemMax= %M' "$@"; }
for t in seqtk art_illumina /usr/bin/time $([ -n "$INDEXES" ] && echo $CL); do command -v $t >/dev/null || fail "missing $t"; done
python3 -c "import regex, requests, aiohttp, pandas, multiprocess" || fail "missing Python modules"

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
  # trav_to_length.txt: genus traversal, then total reference length of that genus
  python3 - <<'PY' || fail trav_to_length
import os
out = open("ref/trav_to_length.txt", "w")
for line in open("ref/rna/doc_to_traversal.txt"):
    doc, trav = line.split(None, 1)
    L = sum(len(l.strip()) for l in open(f"ref/rna/doc_{doc}_seq.fa") if not l.startswith(">"))
    out.write(f"{trav.strip()} {L}\n")
PY
  rm ref/silva_database.fa; echo "$n genera" > ref/done
fi
cat ref/done

for ix in $INDEXES; do
  step "Cliffy index: $ix"
  d=index/$ix; [ -s $d/done ] && continue
  rm -rf $d; mkdir -p $d; cp ref/dna/filelist.txt $d/filelist.txt
  case $ix in
    full_text)  opt="--num-col 7" ;;
    minimizers) opt="--num-col 5 --minimizers --small-window 4 --large-window 11" ;;
    dna_minimizers) opt="--num-col 7 --dna-minimizers --small-window 4 --large-window 11" ;;
  esac
  T --output=$d/time_and_mem.log $CL build --filelist $d/filelist.txt --output $d/output --revcomp --taxcomp $opt \
      --two-pass $d/temp --tmp-size $TMPSIZE 2> $d/build.log || fail "cliffy build $ix (see $W/$d/build.log)"
  rm -rf $d/temp $d/temp.tmp_*; date > $d/done; cat $d/time_and_mem.log
done

mkdir -p abundance
for b in $BIOMES; do   # their MGnify snapshots spell "Candidatus_X" where SILVA 138.1 has "Candidatus X"
  sed 's/Candidatus_/Candidatus /g' $MM/data/$b.tsv > abundance/$b.tsv
  tail -n +2 abundance/$b.tsv | cut -f1 | while read g; do grep -qF "$g	" $TAXTXT || fail "genus not in SILVA 138.1: $g"; done || exit 1
done
for b in $BIOMES; do for r in $REGIONS; do
  step "reads: $b $r ($NREADS pairs)"
  d=reads/$b/$r; [ -s $d/done ] && continue
  rm -rf $d; mkdir -p $d
  p=${r/V4_V4/V4}; P=$MM/data/${p}_primers.txt; [ -s $P ] || fail "no primers $P"
  python3 $MM/src/sim_16s_reads.py simulate --biome-abundance $W/abundance/$b.tsv --primers $P \
      --silva-taxonomy $TAXTXT --silva-ref $FASTA --num-reads $NREADS --output-name $r --temp-dir $W/$d/ > $d/stdout.log 2>&1 \
      || fail "MicrobeMixer $b $r (see $W/$d/stdout.log)"
  rm -f $d/genus_*_reads_1.fq $d/genus_*_reads_2.fq $d/genus_*_seqs.fna $d/genus_*.aln
  [ -s $d/${r}_mate_1.fq ] && [ -s $d/${r}_seqtax.txt ] || fail "no reads for $b $r"
  date > $d/done
done; done

for ix in $INDEXES; do for b in $BIOMES; do for r in $REGIONS; do
  step "Cliffy query + classification: $ix $b $r"
  d=results/$ix/$b/$r; [ -s $d/done ] && continue
  rm -rf $d; mkdir -p $d
  case $ix in
    full_text)  opt="--num-col 7" ;;
    minimizers) opt="--num-col 5 --minimizers --small-window 4 --large-window 11" ;;
    dna_minimizers) opt="--num-col 7 --dna-minimizers --small-window 4 --large-window 11" ;;
  esac
  for m in 1 2; do
    T --output=$d/mate_$m.time $CL run --ref index/$ix/output --pattern reads/$b/$r/${r}_mate_$m.fq \
        --output $d/mate_$m --taxcomp --ftab $opt 2> $d/mate_$m.log || fail "cliffy run $ix $b $r mate $m"
  done
  python3 $EXP/classify_silva_readset.py --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal ref/rna/doc_to_traversal.txt --silva-tax-ranks $TAXTXT \
      --readset-truthset reads/$b/$r/${r}_seqtax.txt --output-dir $W/$d/ --trav-to-length ref/trav_to_length.txt \
      --sample-name output --classify-docprof > $d/classify.log 2>&1 || fail "classify $ix $b $r (see $W/$d/classify.log)"
  date > $d/done
done; done; done
# ---------------------------------------------------------------------------------------------
# rz-index over the same reference: S = doc_1 X rc(doc_1) X ... in tree order (one "genome" per genus,
# its sequences joined by X), then Cliffy-style classification (greedy maximal left-extensions from the
# end of each read; for each match, the leftmost and rightmost genus containing it or its reverse
# complement), scored with their classify_silva_readset.py.  Its "method0" is Cliffy's LCA query.
# EVERY=k classifies every k-th read pair (default 10: 1M of the 10M pairs; the read files are grouped
# by genus, so taking the first reads would be biased).
EVERY=${EVERY:-10}
step "rz-index over the reference"
if [ ! -s rz/bac.aux ]; then
  mkdir -p rz; n=$(cut -d' ' -f1 ref/done)
  [ -s rz/list.tsv ] || for i in $(seq 1 $n); do printf 'silva\t%s\n' "$W/ref/dna/doc_${i}_seq.fa"; done > rz/list.tsv
  OUT=$W/rz PREPFLAGS="" NOMV=1 SVALS="${RZSVALS-}" LENGTHS="" THREADS=6 $RZ/runrz.sh || fail "runrz.sh (rz-index build)"
  [ -s rz/bac.aux ] || fail "rz-index build"
fi
cat rz/rz-build.out
# explicit RLCSA for the classification's backward searches (identical answers, faster than the RLBWT)
[ -s rz/bac.csa ] || $RZ/rz-index/rz-csabuild rz/bac.bwt rz/bac.csa || fail "rz-csabuild"
for b in $BIOMES; do for r in $REGIONS; do
  step "rz-index classification: $b $r (every $EVERY-th pair)"
  d=results/rz/$b/$r; [ -s $d/done ] && continue
  mkdir -p $d
  for m in 1 2; do
    grep -q "^reads=" $d/mate_$m.log 2>/dev/null && { cat $d/mate_$m.log; continue; }   # already classified
    awk -v k=$EVERY 'int((NR-1)/4) % k == 0' reads/$b/$r/${r}_mate_$m.fq > $d/reads_$m.fq
    T --output=$d/mate_$m.time $RZ/rz-index/rz-classify -C rz/bac.csa -x rz/bac.aux rz/bac.rz $d/reads_$m.fq $d/mate_$m.listings 2> $d/mate_$m.log \
        || fail "rz-classify $b $r mate $m (see $W/$d/mate_$m.log)"
    rm $d/reads_$m.fq; cat $d/mate_$m.log
  done
  # scoring: rz_score.py = their method0 (LCA query) without their per-read table copies; it gave the
  # same numbers as classify_silva_readset.py on a 1,500-pair sample, at every level
  python3 $S/rz_score.py --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal ref/rna/doc_to_traversal.txt --silva-tax-ranks $TAXTXT \
      --readset-truthset reads/$b/$r/${r}_seqtax.txt --trav-to-length ref/trav_to_length.txt \
      --output $d/output.classification_results.csv --exp-src $EXP > $d/classify.log 2>&1 || fail "scoring rz $b $r (see $W/$d/classify.log)"
  cat $d/output.classification_results.csv
  date > $d/done
done; done
step "ALL DONE"
