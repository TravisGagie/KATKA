#!/bin/bash
# Shared code of run_katka.sh and run_cliffy.sh (run those, not this): the same reads, timed with one thread and
# with P threads, outputs checked, scored and summarized.  TOOL=katka or TOOL=cliffy.
#
# For every region, every EVERY-th read pair (default 50: 200,000 pairs per region, as in the paper), both mates:
#  1. single thread: each configuration classifies the reads with one thread, and once more a file with a single
#     read, so that loading the index (and building lookup tables) can be subtracted:
#       us/read = (time on all reads - time on one read) / reads.
#     PIN=c runs the single-thread runs on logical CPU c (taskset), e.g. on a performance core of a hybrid CPU;
#  2. parallel: the same reads with P threads.  KATKA: rz-classify -j P, one shared copy of the index
#     (KATKA_SPLIT=1: P processes).  Cliffy, whose queries are single-threaded: P processes on P contiguous slices,
#     outputs concatenated.  Every parallel output must equal the single-thread one.
#     Default P: for Cliffy, the number of logical cores, reduced if P copies of its largest index would not fit in
#     90% of the available memory (written to work/compare/P_cliffy); for KATKA, P_cliffy if run_cliffy.sh has run
#     on this machine (so that both use the same number of threads), else the number of logical cores.  P=... overrides;
#  3. accuracy at the genus level: KATKA with rz_score2.py --mode $MODE (default freq: proportional credit, with
#     RZ_COUNTS=1); Cliffy with rz_score2.py --mode lca (its LCA query, their "method0") and, if THEIRS=1 (default),
#     with their classify_silva_readset.py --classify-docprof (all their methods, including approximate listing;
#     slow).  "strict" counts ties for the maximum as wrong; "theirs" is Cliffy's rule (ties and reads without
#     matches count as correct, except "uncultured").
# Results: work/compare/<tool>-<host>/ (summary.txt, machine.txt, <tool>.log); reads: work/compare/reads/.
# Run on an otherwise idle machine (on a laptop: plugged in, performance power mode).
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); W=$S/work
B=${RZCLASSIFY:-$RZ/rz-index/rz-classify}; CL=${CL:-$S/cliffy/build/cliffy}; EXP=$RZ/cliffy-experiments/src
export PFPDOC_BUILD_DIR=$S/cliffy/build/install
case $TOOL in
  katka)  CLIFFY=""; OURS=${OURS-"ud_L15 ud_L30 dg_L15 dg_L30"} ;;
  cliffy) OURS=""; CLIFFY=${CLIFFY-"full_text minimizers"} ;;
  *) echo "TOOL must be katka or cliffy (run run_katka.sh or run_cliffy.sh)"; exit 1 ;;
esac
REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; EVERY=${EVERY:-50}; MODE=${MODE:-freq}; THEIRS=${THEIRS:-1}
KATKA_SPLIT=${KATKA_SPLIT:-0}; TMPSIZE=${TMPSIZE:-10GB}
RD=$W/compare/reads; R=$W/compare/$TOOL-$(hostname -s); mkdir -p $RD $R/runs $R/scores
exec > >(tee -a $R/$TOOL.log) 2>&1
PIN=${PIN:-}; TS=""; [ -n "$PIN" ] && TS="taskset -c $PIN"
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
TAX=$S/exp1_data/tax_slv_ssu_138.1.txt

# ---- checks -----------------------------------------------------------------------------------------------
for t in /usr/bin/time python3 awk; do command -v $t >/dev/null || fail "missing $t"; done
[ -n "$OURS" ] && { [ -x $B ] || fail "no $B (cd rz-index && make)"; }
[ -s $W/ref/done ] || fail "no reference (run ./run_silva.sh first)"
for r in $REGIONS; do [ -s $W/reads/aquatic/$r/${r}_seqtax.txt ] || fail "no reads for $r (run ./run_silva.sh, or get the paper sample)"; done
ours_opts() {   # options for configuration $1, ending with the index
  local L=${1##*_L}
  case $1 in
    ud_L*) echo "-L $L -l -C $W/rz/bac.csa -T $W/rz/bac.s1.tag -F 10 $W/rz/bac.rz" ;;
    dg_L*) echo "-L $L -l -B $W/rzdg/bac.map -C $W/rzdg/bac.csa -T $W/rzdg/bac.s1.tag -F 2 $W/rzdg/bac.rz" ;;
    *) fail "unknown configuration $1" ;;
  esac
}
ours_files() { case $1 in ud_*) echo "$W/rz/bac.csa $W/rz/bac.s1.tag";; dg_*) echo "$W/rzdg/bac.csa $W/rzdg/bac.s1.tag $W/rzdg/bac.map";; esac; }
for c in $OURS; do for f in $(ours_files $c) $(ours_opts $c | awk '{print $NF}'); do
  [ -s $f ] || fail "missing $f (see the README: run_silva.sh, run_silva_dg.sh, build_indexes.sh)"; done; done
cl_opts() {
  case $1 in
    full_text)  echo "--num-col 7" ;;
    minimizers) echo "--num-col 5 --minimizers --small-window 4 --large-window 11" ;;
    *) fail "unknown Cliffy index $1" ;;
  esac
}
if [ -n "$CLIFFY" ]; then
  [ -x $CL ] || fail "no Cliffy binary $CL (run ./build_cliffy.sh)"
  for ix in $CLIFFY; do
    d=$W/index/$ix; [ -s $d/done ] && continue
    step "building Cliffy index $ix"
    mkdir -p $d; cp $W/ref/dna/filelist.txt $d/filelist.txt; tmp=${TMPDIR_CLIFFY:-$d}/temp_$ix
    /usr/bin/time --output=$d/time_and_mem.log $CL build --filelist $d/filelist.txt --output $d/output --revcomp --taxcomp \
        $(cl_opts $ix) --two-pass $tmp --tmp-size $TMPSIZE 2> $d/build.log || fail "cliffy build $ix (see $d/build.log)"
    rm -rf $tmp $tmp.tmp_*; date > $d/done
  done
fi
if [ ! -s $R/machine.txt ]; then
  { echo "host: $(hostname)"; lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|MHz|L3"; grep -E "MemTotal|MemAvailable" /proc/meminfo
    echo "logical CPUs (core, max MHz; different max MHz = hybrid CPU):"; lscpu --extended=CPU,CORE,MAXMHZ 2>/dev/null | tail -n +2 | tr '\n' ';'; echo
    echo "single-thread runs pinned to CPU: ${PIN:-no}"
    echo "compiler: $(${CXX:-g++} --version | head -1)"; echo "katka: $(git -C $RZ rev-parse --short HEAD 2>/dev/null)"
    echo "cliffy: $(git -C $S/cliffy rev-parse --short HEAD 2>/dev/null)"; } > $R/machine.txt
fi
cat $R/machine.txt

# ---- reads ------------------------------------------------------------------------------------------------
for r in $REGIONS; do for m in 1 2; do
  f=$RD/${r}_$m.fq; [ -s $f ] || $S/reads.sh $r $m $EVERY > $f
  [ -s $RD/${r}_$m.one.fq ] || head -4 $f > $RD/${r}_$m.one.fq
done; done

# ---- running ----------------------------------------------------------------------------------------------
TF='%e %U %S %M'   # elapsed, user, system (s), peak memory (KB)
ours_run() {   # $1 config, $2 threads, $3 reads, $4 output; time to $4.time
  local env=""; [ $MODE = freq ] && env="RZ_COUNTS=1"
  if [ $2 -gt 1 ] && [ $KATKA_SPLIT = 1 ]; then
    env $env /usr/bin/time -f "$TF" -o $4.time $S/run_split.sh $2 $(ours_opts $1 | sed 's/ [^ ]*$//') $(ours_opts $1 | awk '{print $NF}') $3 $4 > $4.log 2>&1
  else
    local ts=""; [ $2 = 1 ] && ts=$TS
    env $env /usr/bin/time -f "$TF" -o $4.time $ts $B -j $2 $(ours_opts $1) $3 $4 2> $4.log
  fi
}
cliffy_run() {   # $1 index, $2 processes, $3 reads, $4 output (writes $4 = the listings); time to $4.time
  local o=$4.cl
  if [ $2 = 1 ]; then
    /usr/bin/time -f "$TF" -o $4.time $TS $CL run --ref $W/index/$1/output --pattern $3 --output $o --taxcomp --ftab $(cl_opts $1) > $4.log 2>&1 \
      && mv $o.listings $4
  else
    local T=$4.slices; rm -rf $T; mkdir -p $T
    local nr=$(( $(wc -l < $3) / 4 )); local per=$(( (nr + $2 - 1) / $2 ))
    awk -v per=$per -v T=$T '{ print > sprintf("%s/in.%04d.fq", T, int((NR-1)/(4*per))) }' $3
    cat > $T/run.sh <<EOF
#!/bin/bash
pids=(); for f in $T/in.*.fq; do i=\${f##*/in.}; i=\${i%.fq}
  $CL run --ref $W/index/$1/output --pattern \$f --output $T/out.\$i --taxcomp --ftab $(cl_opts $1) > $T/log.\$i 2>&1 & pids+=(\$!); done
rc=0; for p in "\${pids[@]}"; do wait \$p || rc=1; done; exit \$rc
EOF
    /usr/bin/time -f "$TF" -o $4.time bash $T/run.sh && cat $T/out.*.listings > $4 && rm -rf $T
  fi
}
# one run, resumable: $1 tool (ours|cliffy), $2 config, $3 threads, $4 region_mate, $5 suffix (all|one)
run1() {
  local o=$R/runs/$2.$4.p$3.$5 in=$RD/$4.fq; [ $5 = one ] && in=$RD/$4.one.fq
  [ -s $o.time ] && [ -s $o ] && return 0
  if [ $1 = ours ]; then ours_run $2 $3 $in $o; else cliffy_run $2 $3 $in $o; fi || fail "$1 $2 p=$3 $4 $5 (see $o.log)"
}

step "single thread"
for r in $REGIONS; do for m in 1 2; do
  for c in $OURS;   do run1 ours $c 1 ${r}_$m one;   run1 ours $c 1 ${r}_$m all; done
  for c in $CLIFFY; do run1 cliffy $c 1 ${r}_$m one; run1 cliffy $c 1 ${r}_$m all; done
done; done

if [ -z "$P" ] && [ $TOOL = katka ] && [ -s $W/compare/P_cliffy ]; then P=$(cat $W/compare/P_cliffy); echo "P = $P (as for Cliffy on this machine)"; fi
if [ -z "$P" ]; then
  P=$(nproc); avail=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
  big=0; for c in $CLIFFY; do for f in $R/runs/$c.*.p1.all.time; do k=$(awk '{print $4}' $f); [ $k -gt $big ] && big=$k; done; done
  if [ $KATKA_SPLIT = 1 ]; then for c in $OURS; do for f in $R/runs/$c.*.p1.all.time; do k=$(awk '{print $4}' $f); [ $k -gt $big ] && big=$k; done; done; fi
  if [ $big -gt 0 ]; then pm=$(( avail * 9 / 10 / big )); [ $pm -lt $P ] && P=$pm; fi
  [ $P -lt 1 ] && P=1
  echo "P = $P (logical cores $(nproc), available memory $((avail / 1024)) MB, largest per-process peak $((big / 1024)) MB)"
fi
echo $P > $R/P; [ $TOOL = cliffy ] && echo $P > $W/compare/P_cliffy

step "parallel, P = $P"
for r in $REGIONS; do for m in 1 2; do
  for c in $OURS;   do run1 ours $c $P ${r}_$m all; done
  for c in $CLIFFY; do run1 cliffy $c $P ${r}_$m all; done
done; done
bad=0
for f in $R/runs/*.p$P.all; do g=${f%.p$P.all}.p1.all; cmp -s $f $g || { echo "OUTPUTS DIFFER: $f $g"; bad=1; }; done
[ $bad = 0 ] && echo "all parallel outputs identical to the single-thread ones"

# ---- accuracy ---------------------------------------------------------------------------------------------
step "accuracy"
score() {   # $1 config, $2 region, $3 mode, $4 output
  [ -s $4 ] && return 0
  python3 $S/rz_score2.py --mode $3 --mate1-listings $R/runs/$1.${2}_1.p1.all --mate2-listings $R/runs/$1.${2}_2.p1.all \
      --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $TAX \
      --readset-truthset $W/reads/aquatic/$2/${2}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
      --output $4 --exp-src $EXP > $4.log 2>&1 || fail "scoring $1 $2 (see $4.log)"
}
for r in $REGIONS; do
  for c in $OURS;   do score $c $r $MODE $R/scores/$c.$r.csv; done
  for c in $CLIFFY; do score $c $r lca $R/scores/$c.$r.csv; done
  if [ $THEIRS = 1 ]; then for c in $CLIFFY; do
    d=$R/scores/$c.$r.theirs; [ -s $d/done ] && continue; mkdir -p $d
    python3 $EXP/classify_silva_readset.py --mate1-listings $R/runs/$c.${r}_1.p1.all --mate2-listings $R/runs/$c.${r}_2.p1.all \
        --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks $TAX \
        --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --output-dir $d/ --trav-to-length $W/ref/trav_to_length.txt \
        --sample-name output --classify-docprof > $d/classify.log 2>&1 || fail "their scoring $c $r (see $d/classify.log)"
    date > $d/done
  done; fi
done

# ---- summary ----------------------------------------------------------------------------------------------
step "summary"
{
  cat $R/machine.txt; echo "reads: every ${EVERY}th pair of $REGIONS (both mates); P = $P"; echo
  printf "%-22s %9s %9s %10s %9s %11s %8s %9s %9s %8s %8s\n" tool/config index_GB load_s us/read peak_GB "P:reads/s" speedup "P:peakGB" "1:reads/s" strict theirs
  for c in $OURS $CLIFFY; do
    if [[ " $OURS " == *" $c "* ]]; then tool=KATKA; sz=$(du -cb $(ours_files $c) | tail -1 | cut -f1); else tool=Cliffy; sz=$(du -cb $W/index/$c/output* | tail -1 | cut -f1); fi
    python3 - "$R" "$c" "$P" "$REGIONS" "$sz" "$tool" "$RD" <<'PY'
import sys, glob, re
R, c, P, regions, sz, tool, RD = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(), int(sys.argv[5]), sys.argv[6], sys.argv[7]
def t(f): e, u, s, m = open(f).read().split()[-4:]; return max(float(e), 0.01), int(m)
def nreads(f): return sum(1 for l in open(f)) // 4
def iw(f):
    try: x = re.search(r"threads=\d+ wall ([0-9.]+) s", open(f).read())
    except OSError: return None
    return float(x.group(1)) if x and float(x.group(1)) > 0 else None
load = []; us = []; pk = []; ppk = []; thr = []; thr1 = []
for r in regions:
    for m in (1, 2):
        b = f"{R}/runs/{c}.{r}_{m}"; n = nreads(f"{RD}/{r}_{m}.fq")
        e1, _ = t(b + ".p1.one.time"); eA, mA = t(b + ".p1.all.time"); eP, mP = t(b + f".p{P}.all.time")
        load.append(e1); us.append(1e6 * (eA - e1) / max(n - 1, 1)); pk.append(mA); ppk.append(mP)
        wP, wA = iw(b + f".p{P}.all.log"), iw(b + ".p1.all.log")   # KATKA reports its own wall time without loading
        thr.append(n / (wP if wP else max(eP - e1, 0.01))); thr1.append(n / (wA if wA else max(eA - e1, 0.01)))
acc = {"strict": [], "theirs": []}
for r in regions:
    for l in open(f"{R}/scores/{c}.{r}.csv"):
        if ",genus," in l:
            for k in acc:
                x = re.search(k + r"=([0-9.]+)", l)
                if x: acc[k].append(float(x.group(1)))
mean = lambda v: sum(v) / len(v) if v else float("nan")
pc = lambda v: (100 * mean(v) if mean(v) <= 1 else mean(v)) if v else float("nan")
print(f"{tool+' '+c:22s} {sz/1e9:9.2f} {mean(load):9.2f} {mean(us):10.1f} {max(pk)/2**20:9.2f} {mean(thr):11.0f} {mean(thr)/mean(thr1):8.1f} {max(ppk)/2**20:9.2f} {mean(thr1):9.0f} {pc(acc['strict']):7.2f}% {pc(acc['theirs']):7.2f}%")
PY
  done
  echo
  echo "us/read: single thread, excluding loading (time on all reads minus time on one read); P:reads/s, 1:reads/s and speedup:"
  echo "wall clock with P threads (KATKA) or P processes (Cliffy) and with one, also excluding loading (KATKA: as"
  echo "  measured by rz-classify; Cliffy: total time minus the time on one read, so only meaningful on many reads); peak: maximum resident memory (GB = 2^30 bytes),"
  echo "for Cliffy with P processes that of one process (P copies in total); accuracy: genus level, mean over the regions,"
  if [ $TOOL = katka ]; then echo "KATKA with $MODE credit."; else echo "Cliffy with its LCA query (method0); its own scoring: $R/scores/<index>.<region>.theirs/"; fi
} | tee $R/summary.txt
