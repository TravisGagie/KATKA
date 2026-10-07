#!/bin/bash
# Throughput of KATKA (default configuration), Kraken 2 and Tagger with 1, 6 and 12 threads on every 50th read pair
# (200,000 pairs per region).  The comparison scripts (run_katka.sh, run_kraken2.sh, run_tagger.sh) use 5,000
# pairs per region, on which runs with 12 threads take well under a second, so their throughputs are dominated by
# start-up; this measures it on 40 times as many reads.  No scoring.  Needs the indexes from build_katka.sh,
# run_kraken2.sh and run_tagger.sh.
#   cd silva && ./run_scale.sh          (summary: work/compare/scale-<host>/summary.txt)
# Reads per second exclude loading the index: KATKA and Kraken 2 report their own classification time; for Tagger
# we subtract the time of a run on a single read.  Environment: TOOLS="katka kraken2 tagger", REGIONS, EVERY=50,
# THREADS="1 6 12".
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); W=$S/work
TOOLS=${TOOLS:-"katka kraken2 tagger"}; REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; EVERY=${EVERY:-50}; THREADS=${THREADS:-"1 6 12"}
KDB=$W/kraken2/silva; TG=$RZ/deps/tagger/build/tagger; TI=$W/tagger/silva; B=$RZ/rz-index/rz-classify
R=$W/compare/scale-$(hostname -s); RD=$W/compare/reads_every$EVERY; mkdir -p $R $RD
exec > >(tee -a $R/scale.log) 2>&1
fail() { echo "FAILED: $*"; exit 1; }
for r in $REGIONS; do for m in 1 2; do
  f=$RD/${r}_$m.fq; [ -s $f ] || { $S/reads.sh $r $m $EVERY > $f.tmp && mv $f.tmp $f; } || fail "reads $r $m"
  [ -s $RD/${r}_$m.one.fq ] || head -4 $f > $RD/${r}_$m.one.fq
done; done
cat $RD/*.fq > /dev/null   # warm the page cache
run() {   # $1 tool, $2 region_mate or region, $3 threads, $4 reads (all|one)
  local o=$R/$1.$2.p$3.$4; [ -s $o.time ] && return 0
  local in=$RD/$2.fq; [ $4 = one ] && in=$RD/$2.one.fq
  case $1 in
    katka)   RZ_COUNTS=1 /usr/bin/time -f "%e" -o $o.time $B -j $3 -L 30 -l -C $W/rz/bac.csa -T $W/rz/bac.s1.tag -F 10 $W/rz/bac.rz $in $o.out 2> $o.log ;;
    kraken2) /usr/bin/time -f "%e" -o $o.time kraken2 --db $KDB --threads $3 --paired $RD/${2}_1.fq $RD/${2}_2.fq > $o.out 2> $o.log ;;
    tagger)  local R_=""; [ $3 -gt 1 ] && R_="-R"
             /usr/bin/time -f "%e" -o $o.time $TG -f $in -r $TI -L 25 -t $3 $R_ -o $o.out > $o.log 2>&1 ;;
  esac || fail "$1 $2 p=$3"
  rm -f $o.out
}
for t in $TOOLS; do for p in $THREADS; do echo "$(date '+%T') $t, $p threads"; for r in $REGIONS; do
  if [ $t = kraken2 ]; then run $t $r $p all; else for m in 1 2; do run $t ${r}_$m $p all; [ $t = tagger ] && run $t ${r}_$m 1 one; done; fi
done; done; done
python3 - "$R" "$TOOLS" "$REGIONS" "$THREADS" "$RD" <<'PY' | tee $R/summary.txt
import sys, re, os
R, tools, regions, threads, RD = sys.argv[1], sys.argv[2].split(), sys.argv[3].split(), sys.argv[4].split(), sys.argv[5]
nreads = lambda f: sum(1 for _ in open(f)) // 4
el = lambda f: float(open(f).read().split()[-1])
print("reads per second excluding loading, mean over the regions (and mates); speedup over one thread")
print("tool       " + "".join(f"{p:>10} thr" for p in threads))
for t in tools:
    rates = {}
    for p in threads:
        v = []
        for r in regions:
            if t == "kraken2":
                m = re.search(r"(\d+) sequences.*processed in ([\d.]+)s", open(f"{R}/{t}.{r}.p{p}.all.log").read())
                v.append(2 * int(m.group(1)) / float(m.group(2)))
            else:
                for mate in (1, 2):
                    b = f"{R}/{t}.{r}_{mate}"; n = nreads(f"{RD}/{r}_{mate}.fq")
                    if t == "katka":
                        w = float(re.search(r"threads=\d+ wall ([0-9.]+) s", open(b + f".p{p}.all.log").read()).group(1))
                    else:
                        w = el(b + f".p{p}.all.time") - el(b + ".p1.one.time")
                    v.append(n / w)
        rates[p] = sum(v) / len(v)
    print(f"{t:<10} " + "".join(f"{rates[p]:>14.0f}" for p in threads))
    print(f"{'':<10} " + "".join(f"{rates[p]/rates[threads[0]]:>13.1f}x" for p in threads))
PY
