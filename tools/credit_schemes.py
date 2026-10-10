#!/usr/bin/env python3
# credit_schemes.py listings truth tbl [seed]
# Scores one classification run (rz-classify -l with RZ_COUNTS=1, so that each MEM is listed as
# [a,b] {tag:count,...}) under several ways of splitting a MEM's credit (its length) among the tags
# (species or genera) whose genomes contain it, and prints the accuracy of each:
#   prop    in proportion to the number of occurrences in each tag (KATKA)
#   even    evenly among the distinct tags containing the MEM
#   one     all to one tag, drawn at random with probability proportional to its occurrences
#           (a stand-in for Tagger's tag of the MEM's last occurrence in the BWT)
#   cliffs  evenly among a sample shaped like Cliffy's cliffs: give each tag containing the MEM a random
#           value and keep the left-to-right and right-to-left maxima in tag order (it always includes the
#           leftmost and rightmost tags, and has about 2 ln d of the d tags; not weighted by occurrences)
#   range   evenly over every tag from the leftmost to the rightmost one containing the MEM, whether it
#           contains it or not (Cliffy's LCA query when tags are in tree order)
# A read goes to the tag with the most credit; ties and reads without MEMs count as wrong.
import re, sys, random, collections
out, truthf, tbl = sys.argv[1:4]; rnd = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
names = []
for l in open(tbl):
    s = l.split('\t')[1]
    if s not in names: names.append(s)
truth = dict(l.split()[:2] for l in open(truthf))
pat = re.compile(r'\[(\d+),(\d+)\] \{([^}]*)\}')
S = ('prop', 'even', 'one', 'cliffs', 'range')
ok = collections.Counter(); n = unc = 0; name = None; ksum = csum = nm = 0

def cliffs(ts):
    v = [rnd.random() for _ in ts]; keep = set(); m = -1.0
    for i, x in enumerate(v):
        if x > m: keep.add(i); m = x
    m = -1.0
    for i in range(len(v) - 1, -1, -1):
        if v[i] > m: keep.add(i); m = v[i]
    return [ts[i] for i in sorted(keep)]

def best(w):
    if not w: return None
    b = max(w.values()); top = [t for t in w if w[t] == b]
    return top[0] if len(top) == 1 else None

for line in open(out):
    line = line.rstrip('\n')
    if line.startswith('>'): name = line[1:]; continue
    n += 1; W = {s: collections.Counter() for s in S}; seen = False
    for a, b, d in pat.findall(line):
        ln = int(b) - int(a) + 1
        cs = sorted((int(x.split(':')[0]), int(x.split(':')[1]) if ':' in x else 1) for x in d.split(',') if x)
        if not cs: continue
        seen = True; ts = [t for t, _ in cs]; tot = sum(c for _, c in cs); nm += 1; ksum += len(ts)
        for t, c in cs: W['prop'][t] += ln * c / tot; W['even'][t] += ln / len(ts)
        r = rnd.random() * tot
        for t, c in cs:
            r -= c
            if r < 0: W['one'][t] += ln; break
        else: W['one'][cs[-1][0]] += ln
        cl = cliffs(ts); csum += len(cl)
        for t in cl: W['cliffs'][t] += ln / len(cl)
        lo, hi = ts[0], ts[-1]
        for t in range(lo, hi + 1): W['range'][t] += ln / (hi - lo + 1)
    if not seen: unc += 1; continue
    for s in S:
        t = best(W[s]); ok[s] += t is not None and names[t] == truth[name]
print('reads %d, unclassified %d (%.2f%%); %d MEMs, %.2f tags per MEM, %.2f in the cliffs-style sample'
      % (n, unc, 100 * unc / n, nm, ksum / max(1, nm), csum / max(1, nm)))
for s in S: print('  %-7s %7.3f%%' % (s, 100 * ok[s] / n))
