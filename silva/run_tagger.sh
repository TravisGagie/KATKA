#!/bin/bash
# Tagger (Depuydt et al., iScience 2025; https://github.com/biointec/tagger) on the paper's reads, on this machine,
# for comparison with KATKA (run_katka.sh), Cliffy and Kraken 2.
#   cd silva && ./run_tagger.sh           (resumable; summary: work/compare/tagger-<host>/summary.txt, or
#   tagger-<host>_L<L>/ with TAGGER_L other than 25; the index does not depend on L and is reused)
# Needs cmake and a C++ compiler; the reference (build_katka.sh) and the reads (get_reads.sh).
#  1. clones and builds Tagger (TAGGER_COMMIT, default: the latest);
#  2. writes the reference as one FASTA file whose headers carry the genus as a tag G00000..G09117 (the same
#     genus numbers as KATKA's documents), with a categories file listing the tags, and builds Tagger's index
#     (its PFP-based build script, with Big-BWT, unless TAGGER_PFP=0; Tagger handles reverse complements at
#     query time, so the reference is indexed one strand only);
#  3. classifies each mate of every 50th read pair of each region with Tagger's defaults (L = 25, 10-mer table),
#     with one thread and then with P threads (-t P -R, output in input order); Tagger gives each read one genus;
#  4. scores the genus level per pair with our scorer (rz_score2.py, listing mode): each mate votes for its genus
#     with its length, so mates that disagree give a tie, which is wrong under the strict measure and correct
#     under Cliffy's rule (we also report the expected accuracy if such ties are broken at random).
#     Tagger itself breaks ties between genera at random, so its outputs with several threads differ slightly.  The time per read is the total time minus the time on a single read (loading).
# Environment: REGIONS, EVERY=50, P, TAGGER_L (default 25), TAGGER_PFP=1, TAGGER_COMMIT.
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); W=$S/work; EXP=$RZ/cliffy-experiments/src
REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; EVERY=${EVERY:-50}; TL=${TAGGER_L:-25}; TAGGER_PFP=${TAGGER_PFP:-1}
TG=$RZ/deps/tagger; TW=$W/tagger
R=$W/compare/tagger-$(hostname -s)$( [ $TL = 25 ] || echo _L$TL ); RD=$W/compare/reads; mkdir -p $R/runs $R/scores $RD $TW
exec > >(tee -a $R/tagger.log) 2>&1
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
for t in cmake git /usr/bin/time python3; do command -v $t >/dev/null || fail "missing $t"; done
[ -s $W/ref/done ] || fail "no reference (run ./build_katka.sh)"

if [ ! -x $TG/build/tagger ]; then
  step "building Tagger"
  [ -d $TG ] || git clone https://github.com/biointec/tagger.git $TG || fail "git clone"
  [ -n "$TAGGER_COMMIT" ] && git -C $TG checkout -q $TAGGER_COMMIT
  (mkdir -p $TG/build && cd $TG/build && cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DSDSL_INCLUDE_DIR=$RZ/deps/sdsl/include -DSDSL_LIBRARY=$RZ/deps/sdsl/lib/libsdsl.a .. > cmake.log 2>&1 && make -j$(nproc) > make.log 2>&1) || fail "building Tagger (see $TG/build/*.log)"
fi
git -C $TG rev-parse --short HEAD > $R/tagger_commit

n=$(cut -d' ' -f1 $W/ref/done)
if [ ! -s $TW/ref_tags.fa ] || head -1 $TW/ref_tags.fa | grep -qv _s; then
  step "reference with genus tags"
  for i in $(seq 1 $n); do
    t=$(printf 'G%05d' $((i - 1)))
    awk -v t=$t '/^>/ { print ">" t "_s" ++k; next } { print }' $W/ref/dna/doc_${i}_seq.fa   # tag + unique number (accessions can contain tags)
  done > $TW/ref_tags.fa.tmp && mv $TW/ref_tags.fa.tmp $TW/ref_tags.fa
  for i in $(seq 1 $n); do printf 'G%05d\n' $((i - 1)); done > $TW/categories.txt
fi
if [ ! -s $TW/build.done ]; then
  step "Tagger index"
  if [ $TAGGER_PFP = 1 ]; then
    (cd $TG/build && /usr/bin/time -f "%e s, %M KB" -o $TW/build.time \
       bash ./tagger_build_pfp.sh -r $TW/silva -t $TW/categories.txt -f $TW/ref_tags.fa > $TW/build.log 2>&1)   # uses paths relative to build/ || fail "Tagger PFP build (see $TW/build.log; try TAGGER_PFP=0)"
  else
    (cd $TW && /usr/bin/time -f "%e s, %M KB" -o build.time $TG/build/tagger_build -f $TW/ref_tags.fa -r $TW/silva -t $TW/categories.txt > build.log 2>&1) \
      || fail "tagger_build (see $TW/build.log)"
  fi
  date > $TW/build.done; cat $TW/build.time
fi

for r in $REGIONS; do for m in 1 2; do
  f=$RD/${r}_$m.fq; [ -s $f ] || $S/reads.sh $r $m $EVERY > $f || fail "reads $r $m"
  [ -s $RD/${r}_$m.one.fq ] || head -4 $f > $RD/${r}_$m.one.fq
done; done
P=${P:-$( [ -s $W/compare/P_cliffy ] && cat $W/compare/P_cliffy || nproc )}; echo $P > $R/P
{ echo "host: $(hostname)"; lscpu | grep -E "Model name|^CPU\(s\)"; grep MemTotal /proc/meminfo; echo "tagger: $(cat $R/tagger_commit), L = $TL"; echo "P = $P"; } > $R/machine.txt

run() {   # $1 region_mate, $2 threads, $3 all|one
  local o=$R/runs/$1.p$2.$3 in=$RD/$1.fq; [ $3 = one ] && in=$RD/$1.one.fq
  [ -s $o.time ] && [ -e $o.out ] && return 0
  local R_=""; [ $2 -gt 1 ] && R_="-R"
  /usr/bin/time -f "%e %U %S %M" -o $o.time $TG/build/tagger -f $in -r $TW/silva -L $TL -t $2 $R_ -o $o.out > $o.log 2>&1 \
    || fail "tagger $1 p=$2 (see $o.log)"
}
step "single thread"; for r in $REGIONS; do for m in 1 2; do run ${r}_$m 1 one; run ${r}_$m 1 all; done; done
step "P = $P threads"; for r in $REGIONS; do for m in 1 2; do run ${r}_$m $P all
  cmp -s $R/runs/${r}_$m.p1.all.out $R/runs/${r}_$m.p$P.all.out || echo "note: ${r}_$m outputs differ with $P threads"; done; done

step "accuracy"
for r in $REGIONS; do
  d=$R/scores/$r; [ -s $d/scores.csv ] && continue; mkdir -p $d
  for m in 1 2; do   # Tagger's "ReadID<tab>Tag" lines -> one listing per read: the whole read, voting for its genus
    python3 - $RD/${r}_$m.fq $R/runs/${r}_$m.p1.all.out > $d/mate_$m.listings <<'PY' || fail "converting $r $m"
import sys, re
lens = {}
with open(sys.argv[1]) as f:
    for i, l in enumerate(f):
        if i % 4 == 0: name = l[1:].split()[0]
        elif i % 4 == 1: lens[name] = len(l.strip())
tag = {}
for l in open(sys.argv[2]):
    t = l.split()
    if len(t) >= 2:
        if t[1].isdigit() and int(t[1]) < 65535:   # index into the categories file (= KATKA's genus); 65535 = no tag
            tag[t[0].lstrip("@>")] = int(t[1])
for name, n in lens.items():
    g = tag.get(name)
    print(">" + name); print(f"[0,{n - 1}] {{{g}}}" if g is not None else "")
PY
  done
  python3 $S/rz_score2.py --mode list --mate1-listings $d/mate_1.listings --mate2-listings $d/mate_2.listings \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt \
      --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $d/scores.csv --exp-src $EXP > $d/score.log 2>&1 || fail "scoring $r (see $d/score.log)"
done

step "summary"
python3 - "$R" "$P" "$REGIONS" "$RD" "$TW" <<'PY' | tee $R/summary.txt
import sys, re, glob, os
R, P, regions, RD, TW = sys.argv[1], sys.argv[2], sys.argv[3].split(), sys.argv[4], sys.argv[5]
t = lambda f: [float(x) for x in open(f).read().split()[-4:]]
nreads = lambda f: sum(1 for _ in open(f)) // 4
size = sum(os.path.getsize(f) for f in glob.glob(TW + "/silva*") if os.path.isfile(f)) / 1e9
us, t1, tP, pk, st, th, tagged = [], [], [], [], [], [], []
for r in regions:
    for m in (1, 2):
        b = f"{R}/runs/{r}_{m}"; n = nreads(f"{RD}/{r}_{m}.fq")
        e1 = t(b + ".p1.one.time")[0]; eA, _, _, mA = t(b + ".p1.all.time"); eP = t(b + f".p{P}.all.time")[0]
        us.append(1e6 * (eA - e1) / n); t1.append(n / max(eA - e1, .01)); tP.append(n / max(eP - e1, .01)); pk.append(mA / 2**20)
    for l in open(f"{R}/scores/{r}/scores.csv"):
        if ",genus," in l:
            st.append(float(re.search(r"strict=([\d.]+)", l).group(1))); th.append(float(l.split(",")[3]))
mean = lambda v: sum(v) / len(v)
print(open(f"{R}/machine.txt").read())
print(f"Tagger index: {size:.2f} GB (files {TW}/silva*); peak memory {max(pk):.2f} GB")
print(f"time per read, one thread, excluding loading: {mean(us):.1f} us; reads/s: 1 thread {mean(t1):.0f}, {P} threads {mean(tP):.0f}")
print(f"genus level, mean over the regions: strict {100*mean(st):.2f}% (range {100*min(st):.2f}-{100*max(st):.2f}); "
      f"mates' disagreements counted correct {100*mean(th):.2f}%; broken at random {50*mean(st)+50*mean(th):.2f}%")
PY
