#!/bin/bash
cd $(dirname $0)/work/dgtest; B=../../../rz-index; T=$1; O=${2:-/tmp/memtest}; mkdir -p $O
for M in "" -M; do for v in rl csa csaef; do c=""; [ $v != rl ] && c="-C bac.$v"
 echo -n "rz $M $v: "; $B/rz-classify $M -B bac.map $c -x bac.aux bac.rz $T $O/rz${M}_$v 2>&1
 for s in 4 16; do echo -n "sr$s $M $v: "; $B/rz-classify $M -B bac.map $c -S bac.s$s.sri -R bac.rix bac.rz $T $O/s$s${M}_$v 2>&1; done; done
 for f in $O/*${M}_rl $O/*${M}_csa $O/*${M}_csaef; do cmp -s $f $O/rz${M}_rl && echo "same $(basename $f)" || echo "DIFF $(basename $f)"; done; done
echo DONE
