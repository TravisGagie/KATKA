#!/bin/bash
# Makes the read samples for the GitHub release and the quick demo from the full simulated reads (run_silva.sh):
#   work/reads/sample50/<region>_mate_<m>.fq.gz     every 50th pair (200,000 per region), as used in the paper
#   work/reads/sample2000/<region>_mate_<m>.fq.gz   every 2,000th pair (5,000 per region), for demo_quick.sh
# with the truth files (<region>_seqtax.txt) in both, and the release file work/katka_silva_sample50.tar.gz.
# Checks that reads.sh gives the same reads from the sample as from the full files.  Takes about 15 minutes.
set -e
S=$(cd "$(dirname "$0")" && pwd); W=$S/work; cd $W/reads
mkdir -p sample50 sample2000
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  cp aquatic/$r/${r}_seqtax.txt sample50/; cp aquatic/$r/${r}_seqtax.txt sample2000/
  for m in 1 2; do
    echo "$r mate $m"
    [ -s sample50/${r}_mate_$m.fq.gz ] || awk 'int((NR-1)/4) % 50 == 0' aquatic/$r/${r}_mate_$m.fq | gzip -9 > sample50/${r}_mate_$m.fq.gz
    [ -s sample2000/${r}_mate_$m.fq.gz ] || gzip -dc sample50/${r}_mate_$m.fq.gz | awk 'int((NR-1)/4) % 40 == 0' | gzip -9 > sample2000/${r}_mate_$m.fq.gz
    # check: every 5,000th pair, from the full file and from the sample
    a=$(awk 'int((NR-1)/4) % 5000 == 0' aquatic/$r/${r}_mate_$m.fq | md5sum)
    b=$(gzip -dc sample50/${r}_mate_$m.fq.gz | awk 'int((NR-1)/4) % 100 == 0' | md5sum)
    [ "$a" = "$b" ] || { echo "MISMATCH $r $m"; exit 1; }
  done
done
cd $W && tar czf katka_silva_sample50.tar.gz reads/sample50 && ls -l katka_silva_sample50.tar.gz && du -sh reads/sample50 reads/sample2000
echo DONE
