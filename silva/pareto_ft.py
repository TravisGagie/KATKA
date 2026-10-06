#!/usr/bin/env python3
# Time/space/accuracy Pareto set from run_pareto_ft.sh (work/pareto_speed.out), with lookup tables.
# Space = structures loaded: backend (RLBWT or RLCSA) + tags or sr samples (+ the RLBWT when an RLCSA is used
# with sr samples, or with sampled tags) + the lookup table (+ the parse index for two-level BML).
# Accuracy: work/results/bml{,_ud}/list_L*/<region>/scores.csv (unchanged by the tables).
import re, os, sys, collections
COUNTS = "--counts" in sys.argv     # proportional credit: times from pareto_speed_counts.out, accuracy from results/freq/*/freq.csv
S = os.path.dirname(os.path.abspath(__file__)); W = S + "/work"; R = ["V1_V2", "V3_V4", "V4_V4", "V4_V5"]
RES = {"rzdg": "bml", "rz": "bml_ud"}
def sizes(ds):
    D = f"{W}/{ds}"; sz = lambda f: os.path.getsize(f"{D}/{f}") / 1e9
    rlbwt = int(re.search(r"rlbwt=(\d+)", open(f"{D}/rz-build.out").read()).group(1)) / 1e9
    back = {"rl": rlbwt, "csa": sz("bac.csa"), "csaef": sz("bac.csaef")}
    pix = sz("bac.k12s4.pix") + sz("bac.k12s4.pix.B") if os.path.exists(f"{D}/bac.k12s4.pix") else 0
    def space(kind, s, v, tab):
        if kind == "sr": x = sz(f"bac.s{s}.sri") + (rlbwt if v != "rl" else 0)
        else: x = sz(f"bac.s{s}.tag") + (rlbwt if s != "1" and v != "rl" else 0)
        if kind == "2l": x += pix
        return back[v] + x + tab
    return space
def acc(ds, L, key):
    v = []
    for r in R:
        f = f"{W}/results/freq/{ds}_L{L}/{r}/freq.csv" if COUNTS else f"{W}/results/{RES[ds]}/list_L{L}/{r}/scores.csv"
        for l in open(f):
            f = l.strip().split(",")
            if f[2] == "genus": d = dict(x.split("=") for x in f[4:]); d["acc"] = f[3]; v.append(float(d[key]))
    return sum(v) / len(v)
t = collections.defaultdict(list); tab = {}; bad = 0
for l in open(f"{W}/pareto_speed_counts.out" if COUNTS else f"{W}/pareto_speed.out"):
    m = re.match(r"(\S+) (\S+) L=(\d+) (tag|sr|2l) s=(\d+) (\S+) F=(\d+): ([\d.]+) us/read(?: table ([\d.]+) MB)?", l)
    if not m: continue
    if "DIFFER" in l: bad += 1
    k = (m.group(1), int(m.group(3)), m.group(4), m.group(5), m.group(6), int(m.group(7)))
    t[k].append(float(m.group(8))); tab[k] = float(m.group(9) or 0) / 1e3
spaces = {ds: sizes(ds) for ds in ("rzdg", "rz")}
rows = []
for k, ts in t.items():
    if len(ts) != 4: continue
    ds, L, kind, s, v, F = k
    rows.append(dict(k=k, time=sum(ts) / 4, space=round(spaces[ds](kind, s, v, tab[k]), 2), strict=acc(ds, L, "strict"), ties=acc(ds, L, "acc")))
def pareto(rows):
    return sorted([a for a in rows if not any(b is not a and b["time"] <= a["time"] and b["space"] <= a["space"] and b["strict"] >= a["strict"]
               and (b["time"] < a["time"] or b["space"] < a["space"] or b["strict"] > a["strict"]) for b in rows)], key=lambda r: (r["space"], r["time"]))
nm = lambda k: f"{'dg' if k[0] == 'rzdg' else 'ud'} L={k[1]:<3} {k[2]} s={k[3]:<2} {k[4]:<5} F={k[5]}"
print(f"answers differing from no table: {bad}")
print("Pareto set (strict accuracy; time us/read, mean of 4 regions; space GB)")
for r in pareto(rows): print(f"  {nm(r['k']):<36} space {r['space']:.2f}  time {r['time']:8.1f}  strict {r['strict']*100:.2f}  ties {r['ties']*100:.2f}")
print("\nall rows")
for r in sorted(rows, key=lambda r: r["k"]): print(f"  {nm(r['k']):<36} {r['space']:.3f} {r['time']:9.1f} {r['strict']*100:.2f}")
