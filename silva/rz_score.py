#!/usr/bin/env python3
# Fast, equivalent re-implementation of "method0" (Cliffy's LCA query) from cliffy-experiments'
# classify_silva_readset.py.  Their version deep-copies a weight table with one entry per clade for
# every read and sends it to worker processes with each read, which takes hours for 1M read pairs on a
# desktop.  This one loads the data with their own classes and applies the same rule with a sparse
# weight table:
#   each match of length L in [left doc, right doc] adds L/(right-left+1) to every document's clade
#   in that range; the read is correct if its true clade is among the clades of maximum weight
#   (so a read with no weight at all counts as correct unless its clade is "uncultured", as in their code); reads whose true clade is
#   missing or ambiguous at a level count as correct, as in their code.
# Output line format is theirs: cliffy,method0,<level>,<fraction correct rounded to 4 places>
import sys, os, argparse, importlib.util
import multiprocessing
# fork, so workers share the loaded read set (Python 3.14 defaults to forkserver on Linux)
Pool = multiprocessing.get_context("fork").Pool

p = argparse.ArgumentParser()
p.add_argument("--mate1-listings", required=True); p.add_argument("--mate2-listings", required=True)
p.add_argument("--doc-id-to-traversal", required=True); p.add_argument("--silva-tax-ranks", required=True)
p.add_argument("--readset-truthset", required=True); p.add_argument("--trav-to-length", required=True)
p.add_argument("--output", required=True); p.add_argument("--exp-src", required=True)
p.add_argument("--threads", type=int, default=os.cpu_count())
a = p.parse_args()
spec = importlib.util.spec_from_file_location("csr", os.path.join(a.exp_src, "classify_silva_readset.py"))
csr = importlib.util.module_from_spec(spec); spec.loader.exec_module(csr)
rs = csr.DocProfReadSet(a.mate1_listings, a.mate2_listings, a.doc_id_to_traversal, a.silva_tax_ranks,
                        a.readset_truthset, a.trav_to_length)

def parse(listing):
    t = listing.split(); out = []
    for i in range(0, len(t), 2):
        s, e = t[i][1:-1].split(","); d = [int(x) for x in t[i + 1][1:-1].split(",")]
        out.append((int(e) - int(s) + 1, d[0], d[-1]))
    return out

G = {}
def work(rng):
    lo, hi = rng; keys = G["keys"]; d2c = G["d2c"]; labels = G["labels"]; level = G["level"]; tp = 0
    for r in rs.read_list[lo:hi]:
        correct = r.correct_taxa.get_certain_level(level)
        n = labels.get(correct, 0)
        if n != 1: tp += 1; continue
        acc = {}
        for listing in (r.mate1_listings, r.mate2_listings):
            for L, left, right in parse(listing):
                w = L / (right + 1 - left)
                for doc in range(left, right + 1):
                    t = d2c.get(doc)
                    if t and t in keys: acc[t] = acc.get(t, 0) + w
        if not acc: tp += correct in keys; continue   # every clade (but "uncultured") ties at weight 0
        m = max(acc.values())
        if m <= 0 or acc.get(correct) == m: tp += 1
    return tp

N = len(rs.read_list)
with open(a.output, "w") as out:
    for level in ["genus", "family", "order", "class", "phylum", "domain"]:
        labels = {}
        for x in rs.silva_nodes_at_each_level[level]:
            c = x.get_certain_level(level); labels[c] = labels.get(c, 0) + 1
        G.update(level=level, labels=labels, keys={c for c in labels if c != "uncultured"},
                 d2c={i: x.get_certain_level(level) for i, x in rs.doc_id_to_full_traversal.items()})
        step = (N + 4 * a.threads - 1) // (4 * a.threads)
        with Pool(a.threads) as pool:           # forked workers see G and rs without pickling them
            tp = sum(pool.map(work, [(i, min(N, i + step)) for i in range(0, N, step)]))
        line = f"cliffy,method0,{level},{round(tp / N, 4)}"
        print(line, flush=True); out.write(line + "\n")
