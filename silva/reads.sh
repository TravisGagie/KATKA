#!/bin/bash
# reads.sh <region> <mate> <k>: writes every k-th read of the region's simulated reads (mate 1 or 2) to stdout.
# Uses the full read files (work/reads/aquatic/<region>/<region>_mate_<m>.fq, from run_silva.sh) if present;
# otherwise the paper's sample, every 50th read pair (work/reads/sample50/<region>_mate_<m>.fq.gz, from the
# release or make_samples.sh), which gives the same reads whenever k is a multiple of 50; or the quick demo's
# sample, every 2,000th pair (work/reads/sample2000), for k a multiple of 2,000.
S=$(cd "$(dirname "$0")" && pwd); W=$S/work; r=$1; m=$2; k=$3
f=$W/reads/aquatic/$r/${r}_mate_$m.fq
if [ -s $f ]; then awk -v k=$k 'int((NR-1)/4) % k == 0' $f; exit; fi
for e in 50 2000; do             # 2000: the quick demo's sample (5,000 pairs per region)
  z=$W/reads/sample$e/${r}_mate_$m.fq.gz
  if [ -s $z ] && [ $((k % e)) = 0 ]; then gzip -dc $z | awk -v k=$((k / e)) 'int((NR-1)/4) % k == 0'; exit; fi
done
echo "reads.sh: no reads for $r mate $m, every ${k}th pair (need $f, or a sample in $W/reads/sample50 or sample2000)" >&2; exit 1
