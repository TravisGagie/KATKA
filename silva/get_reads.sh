#!/bin/bash
# Downloads the paper's read sample (every 50th simulated read pair, 200,000 per region, with truth files; 258 MB)
# from the GitHub release into work/reads/sample50, and puts the truth files where the scripts expect them.
# The experiment scripts read it through reads.sh.  (To simulate all 10 million pairs per region instead, as in
# Cliffy's experiments, run run_silva.sh.)
set -e
S=$(cd "$(dirname "$0")" && pwd); W=$S/work
URL=https://github.com/TravisGagie/KATKA/releases/download/silva-sample/katka_silva_sample50.tar.gz
mkdir -p $W
if [ ! -s $W/reads/sample50/V4_V5_mate_2.fq.gz ]; then
  echo "downloading the read sample (258 MB)"
  wget -q -O $W/katka_silva_sample50.tar.gz $URL
  echo "168ffc65ca1070f730272a16bc6532cb  $W/katka_silva_sample50.tar.gz" | md5sum -c --quiet
  tar xzf $W/katka_silva_sample50.tar.gz -C $W && rm $W/katka_silva_sample50.tar.gz
fi
for r in V1_V2 V3_V4 V4_V4 V4_V5; do mkdir -p $W/reads/aquatic/$r; cp -n $W/reads/sample50/${r}_seqtax.txt $W/reads/aquatic/$r/; done
ls $W/reads/sample50
