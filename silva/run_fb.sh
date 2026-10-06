#!/bin/bash
# Forward-backward (rz-classify -W: every MEM found, only those of at least L bases listed) against BML on the
# undigested SILVA, with unsampled tags on the explicit RLCSA, without and with lookup tables (-F) and two-level
# search (-P, k = 12, s = 4).  2,000 reads per region, L = 15, 20, 30, 40; the answers must be identical.
# Run with the Claude app closed.  Output: work/rz/fb_speed.out
S=$(cd "$(dirname "$0")" && pwd); RZ=$(dirname "$S"); B=$RZ/rz-index/rz-classify; W=$S/work; D=$W/rz
OUT=$D/fb_speed.out; : > $OUT
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  $S/reads.sh $r 1 5000 > $D/fb_reads.fq
  for L in 15 20 30 40; do
    $B -L $L -l -C $D/bac.csa -T $D/bac.s1.tag $D/bac.rz $D/fb_reads.fq $D/fb_base 2>/dev/null
    while read name opts; do
      echo -n "$r L=$L $name: " >> $OUT
      $B -L $L -l -C $D/bac.csa -T $D/bac.s1.tag $opts $D/bac.rz $D/fb_reads.fq $D/fb_out 2>&1 | grep -o "steps/read=[0-9.]*\|[0-9.]* us/read" | tr '\n' ' ' >> $OUT
      cmp -s $D/fb_base $D/fb_out && echo "[identical]" >> $OUT || echo "[ANSWERS DIFFER]" >> $OUT
    done <<LIST
bml
bml-F10     -F 10
bml-F12     -F 12
bml-2l-F12  -P $D/bac.k12s4.pix -F 12
fb          -W
fb-F10      -W -F 10
fb-F12      -W -F 12
fb-2l       -W -P $D/bac.k12s4.pix
fb-2l-F12   -W -P $D/bac.k12s4.pix -F 12
LIST
  done
done
rm -f $D/fb_reads.fq $D/fb_base $D/fb_out
python3 - "$OUT" <<'PY'
import re, sys, collections
d = collections.defaultdict(list)
for l in open(sys.argv[1]):
    m = re.match(r'(\S+) L=(\d+) (\S+): .*?([\d.]+) us/read', l)
    if m: d[(m.group(3), int(m.group(2)))].append(float(m.group(4)))
names = []
for (n, L) in d:
    if n not in names: names.append(n)
print('config'.ljust(12), *[f'L={L}'.rjust(8) for L in (15, 20, 30, 40)])
for n in names: print(n.ljust(12), *[f'{sum(d[(n, L)]) / len(d[(n, L)]):8.1f}' for L in (15, 20, 30, 40)])
PY
echo "$(grep -c identical $OUT) identical, $(grep -c DIFFER $OUT) differ"
