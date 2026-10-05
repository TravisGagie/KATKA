#!/usr/bin/env python3
# Time/space/accuracy table and Pareto set for run_bml.sh (digested SILVA) and run_bml_ud.sh (undigested).
# usage: pareto_bml.py [dg] [ud]     (default dg; with both, one combined Pareto set)
import re, os, sys, collections
S = os.path.dirname(os.path.abspath(__file__)); W = S + "/work"
R = ["V1_V2", "V3_V4", "V4_V4", "V4_V5"]
DS = sys.argv[1:] or ["dg"]
INFO = {"dg": ("rzdg", "bml"), "ud": ("rz", "bml_ud")}     # index dir, results dir

def load(ds):
    D = f"{W}/{INFO[ds][0]}"; RES = f"{W}/results/{INFO[ds][1]}"
    sz = lambda f: os.path.getsize(f"{D}/{f}") / 1e9
    rlbwt = int(re.search(r"rlbwt=(\d+)", open(f"{D}/rz-build.out").read()).group(1)) / 1e9   # RLBWT inside bac.rz
    rzrest = sz("bac.rz") - rlbwt
    back = {"rl": rlbwt, "csa": sz("bac.csa"), "csaef": sz("bac.csaef")}
    def space(cfg):
        if cfg[0] == "tag":                  # tags + backend (+ the RLBWT for LF when the tags are sampled)
            s, v = cfg[1], cfg[2]
            return sz(f"bac.s{s}.tag") + back[v] + (rlbwt if s != "1" and v != "rl" else 0)
        if cfg[0] == "rz": return rzrest + back[cfg[1]] + (sz("bac.aux") if cfg[2] == "aux" else 0)
        s, v = cfg[1], cfg[2]                # sr: sampling + backend (+ the RLBWT's run heads/starts for an RLCSA)
        return sz(f"bac.s{s}.sri") + back[v] + (rlbwt if v != "rl" else 0)
    def acc(L, mode, key):
        v = []
        for r in R:
            for l in open(f"{RES}/{mode}_L{L}/{r}/scores.csv"):
                f = l.strip().split(",")
                if f[2] == "genus":
                    d = dict(x.split("=") for x in f[4:]); d["acc"] = f[3]; v.append(float(d[key]))
        return sum(v) / len(v)
    t = collections.defaultdict(list)        # time: mean us/read over the four regions
    for l in open(D + "/bml_speed.out"):
        m = re.match(r"(\S+) L=(\d+) (rz (\S+) (aux|noaux)|sr s=(\d+) (\S+) (lca|list)): .*? ([\d.]+) us/read", l)
        if not m: continue
        cfg = ("rz", m.group(4), m.group(5), "lca") if m.group(4) else ("sr", m.group(6), m.group(7), m.group(8))
        t[(int(m.group(2)), cfg)].append(float(m.group(9)))
    if os.path.exists(D + "/tag_speed.out"):
        for l in open(D + "/tag_speed.out"):
            m = re.match(r"(\S+) L=(\d+) tag s=(\d+) (\S+) (list|lca): .*? ([\d.]+) us/read", l)
            if m: t[(int(m.group(2)), ("tag", m.group(3), m.group(4), m.group(5)))].append(float(m.group(6)))
    rows = []
    for (L, cfg), ts in t.items():
        if len(ts) != 4: continue
        rows.append(dict(ds=ds, L=L, cfg=cfg, time=sum(ts) / 4, space=space(cfg),
                         strict=acc(L, cfg[3], "strict"), ties=acc(L, cfg[3], "acc")))
    return rows

rows = [r for ds in DS for r in load(ds)]
def pareto(rows, key):
    out = []
    for a in rows:
        if not any(b is not a and b["time"] <= a["time"] and b["space"] <= a["space"] and b[key] >= a[key]
                   and (b["time"] < a["time"] or b["space"] < a["space"] or b[key] > a[key]) for b in rows): out.append(a)
    return sorted(out, key=lambda r: (r["space"], r["time"]))
cn = lambda c: f"rz {c[1]} {c[2]}" if c[0] == "rz" else f"{c[0]} s={c[1]} {c[2]} {c[3]}"
name = lambda r: f"{r['ds']} {cn(r['cfg'])}"
for key in ("strict", "ties"):
    print(f"\nPareto set (accuracy = {key}, genus, mean of 4 regions; time = us/read; space = GB)")
    for r in pareto(rows, key):
        print(f"  L={r['L']:<3} {name(r):<25} space {r['space']:.2f}  time {r['time']:8.1f}  strict {r['strict']*100:.2f}  ties {r['ties']*100:.2f}")
print("\nall rows: L cfg space time strict ties")
for r in sorted(rows, key=lambda r: (r["ds"], r["L"], r["cfg"])):
    print(f"  {r['L']:<3} {name(r):<25} {r['space']:.2f} {r['time']:9.1f} {r['strict']*100:.2f} {r['ties']*100:.2f}")
