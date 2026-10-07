#!/bin/bash
# Kraken 2 on the paper's reads, on this machine, for comparison with KATKA (run_katka.sh) and Cliffy.
#   cd silva && ./run_kraken2.sh          (resumable; summary: work/compare/kraken2-<host>/summary.txt)
# Needs kraken2 and kraken2-build in the PATH (e.g. sudo apt install kraken2, or conda install -c bioconda kraken2),
# the reference (build_katka.sh) and the reads (get_reads.sh).
#  1. builds Kraken 2's SILVA database with its own script (kraken2-build --special silva: SILVA 138.1 SSU NR99,
#     the release we use; default k = 35, minimizer length 31), as Cliffy's experiments did;
#  2. classifies every 50th read pair of each region (both mates together, --paired, as in Cliffy's scripts)
#     with one thread and then with P threads (default: the P saved by run_cliffy.sh, else all logical cores);
#     the time per read is Kraken 2's own "processed in" time (excluding loading the database) divided by the
#     number of reads (two per pair);
#  3. scores the genus-level accuracy with Cliffy's own Kraken 2 scoring (classify_silva_readset.py, its
#     KrakenReadSet class), without Bracken, which only estimates abundances.
# Environment: REGIONS, EVERY=50, P, THREADS_BUILD (default: all cores), KDB (default work/kraken2/silva).
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); W=$S/work; EXP=$RZ/cliffy-experiments/src
REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; EVERY=${EVERY:-50}
KDB=${KDB:-$W/kraken2/silva}; THREADS_BUILD=${THREADS_BUILD:-$(nproc)}
R=$W/compare/kraken2-$(hostname -s); RD=$W/compare/reads; mkdir -p $R/runs $R/scores $RD
exec > >(tee -a $R/kraken2.log) 2>&1
step() { echo "=== $(date '+%F %T') $*"; }
fail() { echo "FAILED: $*"; exit 1; }
for t in kraken2 kraken2-build /usr/bin/time python3; do command -v $t >/dev/null || fail "missing $t"; done
[ -s $W/ref/done ] || fail "no reference (run ./build_katka.sh)"
TAX=$S/exp1_data/tax_slv_ssu_138.1.txt

if [ ! -s $KDB/hash.k2d ]; then
  step "building Kraken 2's SILVA database in $KDB"
  mkdir -p $(dirname $KDB)
  /usr/bin/time -f "%e s, %M KB" -o $(dirname $KDB)/build.time \
    kraken2-build --special silva --db $KDB --threads $THREADS_BUILD > $(dirname $KDB)/build.log 2>&1 || fail "kraken2-build (see $(dirname $KDB)/build.log)"
  cat $(dirname $KDB)/build.time
fi
du -sh $KDB/*.k2d

for r in $REGIONS; do for m in 1 2; do
  f=$RD/${r}_$m.fq; [ -s $f ] || $S/reads.sh $r $m $EVERY > $f || fail "reads $r $m"
done; done
P=${P:-$( [ -s $W/compare/P_cliffy ] && cat $W/compare/P_cliffy || nproc )}; echo $P > $R/P
{ echo "host: $(hostname)"; lscpu | grep -E "Model name|^CPU\(s\)"; grep MemTotal /proc/meminfo; kraken2 --version | head -1; echo "P = $P"; } > $R/machine.txt

run() {   # $1 region, $2 threads
  local o=$R/runs/$1.p$2
  [ -s $o.out ] && [ -s $o.log ] && return 0
  /usr/bin/time -f "%e %U %S %M" -o $o.time kraken2 --db $KDB --threads $2 --paired $RD/${1}_1.fq $RD/${1}_2.fq > $o.out 2> $o.log \
    || fail "kraken2 $1 p=$2 (see $o.log)"
}
step "single thread"; for r in $REGIONS; do run $r 1; done
step "P = $P threads"; for r in $REGIONS; do run $r $P; cmp -s $R/runs/$r.p1.out $R/runs/$r.p$P.out || echo "note: $r outputs differ in order or content with $P threads"; done

step "accuracy"
for r in $REGIONS; do
  d=$R/scores/$r; [ -s $d/output.classification_results.csv ] && continue; mkdir -p $d
  python3 - "$EXP" "$TAX" "$W/reads/aquatic/$r/${r}_seqtax.txt" "$R/runs/$r.p1.out" "$d/" > $d/score.log 2>&1 <<'PY' || fail "scoring $r (see $d/score.log)"
import sys, importlib.util
exp, tax, truth, kout, outdir = sys.argv[1:6]
spec = importlib.util.spec_from_file_location("csr", exp + "/classify_silva_readset.py")
csr = importlib.util.module_from_spec(spec); spec.loader.exec_module(csr)
rs = csr.KrakenReadSet(tax, truth, kout, None)     # Bracken is only used for abundances
rs.generate_sensitivity_plot(outdir, "output")
PY
done

step "summary"
python3 - "$R" "$P" "$REGIONS" "$KDB" <<'PY' | tee $R/summary.txt
import sys, re, glob, os
R, P, regions, kdb = sys.argv[1], sys.argv[2], sys.argv[3].split(), sys.argv[4]
def proc(f):            # "N sequences (... Mbp) processed in T s"
    m = re.search(r"(\d+) sequences.*processed in ([\d.]+)s", open(f).read()); return int(m.group(1)), float(m.group(2))
size = sum(os.path.getsize(f) for f in glob.glob(kdb + "/*.k2d")) / 1e9
us, thr1, thrP, pk = [], [], [], []
for r in regions:
    n, t1 = proc(f"{R}/runs/{r}.p1.log"); _, tP = proc(f"{R}/runs/{r}.p{P}.log")
    us.append(1e6 * t1 / (2 * n)); thr1.append(2 * n / t1); thrP.append(2 * n / tP)
    pk.append(int(open(f"{R}/runs/{r}.p1.time").read().split()[-1]) / 2**20)
acc = {}
for r in regions:
    for l in open(f"{R}/scores/{r}/output.classification_results.csv"):
        f = l.strip().split(",")
        if len(f) == 4 and f[2] == "genus": acc.setdefault(f[0], {})[r] = float(f[3])
mean = lambda v: sum(v) / len(v)
print(open(f"{R}/machine.txt").read())
print(f"Kraken 2 SILVA database: {size:.2f} GB; peak memory {max(pk):.2f} GB")
print(f"time per read, one thread, excluding loading: {mean(us):.1f} us; reads/s: 1 thread {mean(thr1):.0f}, {P} threads {mean(thrP):.0f}")
print("genus level, Cliffy's Kraken 2 scoring; without_vp: only pairs assigned to their true genus count as correct;")
print("with_vp: pairs assigned to a correct higher taxon also count (Cliffy's 'vague positives')")
for k, v in sorted(acc.items()):
    vals = [v[r] for r in regions]
    print(f"  {k}: mean {100*mean(vals):.2f}%, range {100*min(vals):.2f}-{100*max(vals):.2f}%  ({', '.join(f'{r} {100*v[r]:.2f}' for r in regions)})")
PY
