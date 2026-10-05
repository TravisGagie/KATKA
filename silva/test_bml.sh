#!/bin/bash
# quick check of rz-classify -L (BML on the digest) and -l (listing): all backends must agree
W=$(cd "$(dirname "$0")" && pwd)/work/rzdg; B=$(cd "$(dirname "$0")/.." && pwd)/rz-index/rz-classify
T=$1; O=$2; L=${3:-20}; mkdir -p $O
for v in rl csa csaef; do c=""; [ $v != rl ] && c="-C $W/bac.$v"
  echo -n "rz $v: "; $B -L $L -B $W/bac.map $c -x $W/bac.aux $W/bac.rz $T $O/rz_$v 2>&1
  for s in 0 4 16; do
    echo -n "sr$s $v lca: "; $B -L $L -B $W/bac.map $c -S $W/bac.s$s.sri -R $W/bac.rix $W/bac.rz $T $O/s${s}_$v 2>&1
    echo -n "sr$s $v list: "; $B -L $L -l -B $W/bac.map $c -S $W/bac.s$s.sri -R $W/bac.rix $W/bac.rz $T $O/l${s}_$v 2>&1
  done
done
for f in $O/rz_* $O/s*_*; do cmp -s $f $O/rz_rl && echo "same $(basename $f)" || echo "DIFF $(basename $f)"; done
for f in $O/l*_*; do cmp -s $f $O/l0_rl && echo "same $(basename $f)" || echo "DIFF $(basename $f)"; done
# listing ends == LCA ends
python3 - $O/rz_rl $O/l0_rl <<'PY'
import sys,re
a=open(sys.argv[1]).read().split('\n'); b=open(sys.argv[2]).read().split('\n'); bad=0; n=0
for x,y in zip(a,b):
    if x.startswith('>'): continue
    ex=re.findall(r'\[(\d+),(\d+)\] \{(\d+),(\d+)\}',x); ey=re.findall(r'\[(\d+),(\d+)\] \{([\d,]+)\}',y)
    if len(ex)!=len(ey): bad+=1; continue
    for (s,e,l,r),(s2,e2,d) in zip(ex,ey):
        d=d.split(','); n+=1
        if (s,e)!=(s2,e2) or d[0]!=l or d[-1]!=r: bad+=1
print("listing vs LCA ends: features", n, "mismatches", bad)
PY
echo DONE
