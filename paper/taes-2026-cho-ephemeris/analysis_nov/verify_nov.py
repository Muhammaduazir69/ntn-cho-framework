#!/usr/bin/env python3
"""Independent recheck of macros_nov.tex from the raw run outputs. Written separately from
analyze_nov.py: its own CSV parsing, its own Kaplan-Meier (product-limit over sorted unique
times), its own late-share and paired-difference code. Prints PASS/FAIL per macro checked
and writes verify_nov_report.txt."""
import csv, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
PAPER = os.path.dirname(SIM)
NOV = os.path.join(SIM, "runs_nov")
mac = {}
for line in open(os.path.join(PAPER, "macros_nov.tex")):
    m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s*$", line)
    if m:
        mac[m.group(1)] = m.group(2)
out, nfail = [], 0


def check(name, val, d, tol=None):
    global nfail
    if name not in mac:
        out.append(f"SKIP {name} (no macro)")
        return
    try:
        printed = float(mac[name].replace("\\ensuremath{-}", "-"))
    except ValueError:
        out.append(f"SKIP {name} = {mac[name]}")
        return
    tol = tol if tol is not None else 0.5 * 10 ** (-d) + 1e-9
    ok = abs(printed - val) <= tol
    nfail += not ok
    out.append(f"{'PASS' if ok else 'FAIL'} {name}: printed {mac[name]}, recomputed {val:.{d + 2}f}")


def km_restricted_mean(samples):
    """samples: list of (duration, observed)"""
    times = sorted({t for t, _ in samples})
    n_at = len(samples)
    S, area, prev = 1.0, 0.0, 0.0
    by = defaultdict(lambda: [0, 0])
    for t, o in samples:
        by[t][0] += 1
        by[t][1] += o
    for t in times:
        area += S * (t - prev)
        prev = t
        tot, ev = by[t]
        S *= (1 - ev / n_at)
        n_at -= tot
    return area


def run_rate(d):
    with open(os.path.join(d, "summary.csv")) as f:
        row = next(csv.DictReader(f))
    H, U = float(row["sim_time_s"]), int(row["num_ues"])
    hos = defaultdict(list)
    with open(os.path.join(d, "events.csv")) as f:
        for e in csv.DictReader(f):
            if e["kind"] == "ho" and e["success"] == "1":
                hos[int(e["ue"])].append(float(e["t_s"]))
    smp = []
    for ts in hos.values():
        ts.sort()
        smp += [(b - a, 1) for a, b in zip(ts, ts[1:])]
        smp.append((H - ts[-1], 0))
    return 60.0 / km_restricted_mean(smp), float(row["outage_frac"]) * 100, float(row["mean_sinr_db"])


def late(d):
    n = l = 0
    with open(os.path.join(d, "events.csv")) as f:
        for e in csv.DictReader(f):
            if e["kind"] == "ho":
                n += 1
                l += e["true_serv_exit_s"] not in ("", None) and float(e["true_serv_exit_s"]) <= 0
    return 100.0 * l / n


def seeds(p):
    return {int(s[4:]): os.path.join(p, s) for s in os.listdir(p) if s.startswith("seed")} if os.path.isdir(p) else {}


def mean(v):
    return sum(v) / len(v)


CH = {"los": "Los", "suburban": "Sub", "urban": "Urb", "denseurban": "Dense"}
PN = {"a3tf": "Atuned", "d2f": "Dtwof", "t1f": "Tonef", "zabf": "Rank", "mst": "Mst", "ut": "Ut",
      "tted2f": "Tted", "geof2": "Geo", "geof": "GeoTwoh", "geof10": "Near", "mst_k16": "MstK", "mst40": "Mstforty", "ut40": "Utforty"}
R = {}
for env in CH:
    for pol in PN:
        S = seeds(os.path.join(NOV, "conf", env, pol))
        if len(S) == 10:
            R[(env, pol)] = {s: run_rate(d) for s, d in S.items()}
            check(f"Nv{CH[env]}{PN[pol]}Eone", mean([v[0] for v in R[(env, pol)].values()]), 3)
            check(f"Nv{CH[env]}{PN[pol]}Out", mean([v[1] for v in R[(env, pol)].values()]),
                  3 if env in ("los", "suburban") else 2)
for env in ("los", "suburban"):
    for b in ("tted2f", "geof2", "zabf", "mst"):
        if (env, "ut") in R and (env, b) in R:
            d = [100 * (1 - R[(env, "ut")][s][0] / R[(env, b)][s][0]) for s in R[(env, "ut")]]
            check(f"Nv{CH[env]}UtVs{PN[b]}", mean(d), 1)
if ("los", "mst_k16") in R:
    d = [100 * (1 - R[("los", "ut")][s][0] / R[("los", "mst_k16")][s][0]) for s in R[("los", "ut")]]
    check("NvLosUtVsMstK", mean(d), 1)
for ch in ("los", "suburban"):
    for pol, nm in (("mst", "Mst"), ("tte", "Tte")):
        S = seeds(os.path.join(NOV, "oos", ch, pol))
        if len(S) == 10:
            check(f"NvOos{CH[ch]}{nm}Late", mean([late(d) for d in S.values()]), 2)
S1, S2 = seeds(os.path.join(NOV, "soft", "los", "mst")), seeds(os.path.join(NOV, "soft", "los", "tte"))
if len(S1) == 10 and len(S2) == 10:
    o = lambda d: float(next(csv.DictReader(open(os.path.join(d, "summary.csv"))))["outage_frac"]) * 100
    check("NvSoftMstOut", mean([o(d) for d in S1.values()]), 3)
    check("NvSoftSave", mean([o(S1[s]) - o(S2[s]) for s in S1]), 3)
# Amendment 4: per-pass vs median-matched fixed margin, late share recomputed from events
for k, n in (("r1", "One"), ("r3", "Three"), ("r5", "Five")):
    for arm, mname in (("tte", "LatePer"), ("fixedmed", "LateMed"), ("fixed", "LateFix")):
        S = seeds(os.path.join(NOV, "perpass", k, arm))
        if len(S) == 10:
            check(f"NvPp{n}{mname}", mean([late(d) for d in S.values()]), 2)
# Amendment 5: offset-model per-pass test, and the along-track drift of the propagated copy
for k, n in (("r1", "One"), ("r3", "Three"), ("r5", "Five")):
    for arm, mname in (("tte", "LatePer"), ("fixedmed", "LateMed")):
        S = seeds(os.path.join(NOV, "perpass_off", k, arm))
        if len(S) == 10:
            check(f"NvPoff{n}{mname}", mean([late(d) for d in S.values()]), 2)
    vals = []
    for d in seeds(os.path.join(NOV, "perpass", k, "tte")).values():
        vals += [float(r["dT_a_m"]) for r in csv.DictReader(open(os.path.join(d, "eph_offsets.csv")))]
    if vals:
        mu = sum(vals) / len(vals)
        check(f"NvDrift{n}Obs", math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals)) / 1e3, 1)
# Calibration of the exit-error lemma in the drift-free runs: sd of eps / sigma_tau
for k, n in (("r1", "One"), ("r3", "Three"), ("r5", "Five")):
    zz = []
    for d in seeds(os.path.join(NOV, "perpass_off", k, "tte")).values():
        for e in csv.DictReader(open(os.path.join(d, "events.csv"))):
            if e["kind"] == "ho" and e["trigger"] == "planned" and e["true_serv_exit_s"] not in ("", None):
                sg = (float(e["margin_s"]) - 0.2) / 2.0
                zz.append((float(e["pred_serv_tte_s"]) - float(e["true_serv_exit_s"])) / sg)
    if zz:
        mu = sum(zz) / len(zz)
        check(f"NvPoff{n}ZSd", math.sqrt(sum((v - mu) ** 2 for v in zz) / len(zz)), 2)
rep = "\n".join(out) + f"\n{len(out)} lines, {nfail} FAIL\n"
open(os.path.join(HERE, "verify_nov_report.txt"), "w").write(rep)
print(rep)
sys.exit(1 if nfail else 0)
