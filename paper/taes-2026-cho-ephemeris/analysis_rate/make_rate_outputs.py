#!/usr/bin/env python3
"""Macros (macros_rate.tex) and Fig. fig_rate.pdf for Section VII-A of the TAES manuscript:
handover-rate analysis (rate_model.py) against the exact-geometry Walker Monte Carlo
(mc_rate.csv, written by `rate_model.py --mc 40`) and against ns-3 (results_macros.tex).

    python3 make_rate_outputs.py
"""
import ast
import csv
import json
import math
import os
import re

import numpy as np

import rate_model as M
import walker_mc as W

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.abspath(os.path.join(HERE, "..", ".."))
DENSE_NBAR = 3.0  # "dense" shells for the error range of the exit-triggered analysis


def ns3_macros():
    m = {}
    for line in open(os.path.join(PAPER, "results_macros.tex")):
        g = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s*$", line)
        if g:
            m[g.group(1)] = g.group(2)
    return m


def ns3_geodesic():
    from scipy import stats
    root = os.path.join(HERE, "..", "runs_r2")
    paths = {792: "density/w792/los/geof", 1584: "main/los/geof", 3168: "density/w3168/los/geof"}
    out = {}
    for n, rel in paths.items():
        v = [float(next(csv.DictReader(open(f)))["ho_per_ue_min"])
             for f in sorted(__import__("glob").glob(os.path.join(root, rel, "seed*", "summary.csv")))]
        if len(v) >= 2:
            v = np.array(v)
            out[n] = (float(v.mean()), float(stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))))
    return out if len(out) == 3 else {}


def load_mc():
    # the Monte Carlo ran in two parts (resumed with --skip 3); merge the printed rows
    rows = []
    lines = []
    for fn in ("mc_rate_part1.log", "mc_rate_part2.log"):
        if os.path.exists(os.path.join(HERE, fn)):
            lines += open(os.path.join(HERE, fn)).readlines()
    for line in lines:
        if line.startswith("{"):
            d = ast.literal_eval(re.sub(r"np\.float64\(([^)]*)\)", r"\1", line.strip()))
            if "mc_ER" in d:
                rows.append({k: float(v) for k, v in d.items()})
    return rows


def input_check(N=1584, P=72, F=17, h=550, inc=53, eQ=25, epochs=6, T=6000):
    """Pass statistics of the real Walker geometry, for the inputs of Lemma 1."""
    cache = os.path.join(HERE, "input_check.json")
    if os.path.exists(cache):
        return json.load(open(cache))
    rng = np.random.default_rng(3)
    raan, u0, n, r, i = W.shell(N, P, F, h, inc)
    durs, nv, npass = [], [], 0
    for _ in range(epochs):
        u = u0 + n * rng.uniform(0, 86400)
        t = np.arange(0, T, 1.0)
        vis = W.elev(t, raan, u, n, r, i, np.radians(30), rng.uniform(0, 2 * np.pi)) >= eQ
        nv.append(vis.sum(0).mean())
        d = np.diff(vis.astype(int), axis=1)
        for s in range(N):
            st, en = np.where(d[s] == 1)[0], np.where(d[s] == -1)[0]
            npass += len(st)
            for a in st:
                b = en[en > a]
                if len(b):
                    durs.append(float(b[0] - a))
    out = dict(Nbar=float(np.mean(nv)), lam=npass / (epochs * T), ET=float(np.mean(durs)))
    json.dump(out, open(cache, "w"))
    return out


def main():
    ns = ns3_macros()
    mc = load_mc()
    key = lambda r: (int(r["N"]), int(r["eQ"]), int(r["h"]))
    an = {}
    for c in M.CASES:
        a = M.rates(c[0], c[3], c[4], c[5])
        a["Hc"] = 60 / M.poisson_chain(c[0], c[3], c[4], c[5], K=20000)
        a["Ha"] = 60 / M.asym(c[0], c[3], c[4], c[5])
        an[(c[0], c[5], c[3])] = a
    mcd = {key(r): r for r in mc}
    out = {}
    ref = an[(1584, 25, 550)]
    inp = input_check()
    out["RateAnNbar"] = f"{ref['Nvis']:.1f}"
    out["RateMcNbar"] = f"{inp['Nbar']:.1f}"
    out["RateAnLam"] = f"{60 * ref['lam']:.2f}"
    out["RateMcLam"] = f"{60 * inp['lam']:.2f}"
    out["RateAnET"] = f"{ref['ET']:.0f}"
    out["RateMcET"] = f"{inp['ET']:.0f}"
    out["RateNearAnRef"] = f"{ref['Hnear']:.3f}"
    out["RateExitAnRef"] = f"{ref['Hc']:.3f}"
    out["RateFortyAnRef"] = f"{an[(1584, 40, 550)]['Hc']:.3f}"
    out["RateFloorRef"] = f"{ref['Hinf']:.3f}"
    out["RateRatioRef"] = f"{ref['Hnear'] / ref['Hc']:.1f}"
    out["RateDenseNbar"] = f"{DENSE_NBAR:.0f}"
    geo, mst = float(ns["MainLosGeofHo"]), float(ns["MainLosMstHo"])
    out["RateNsRatio"] = f"{geo / mst:.1f}"
    out["RateNsNearErr"] = f"{100 * abs(ref['Hnear'] - geo) / geo:.0f}"
    # link margin at the coverage edge (LOS, deterministic part of eq. snr)
    Rk, h, f = W.R, 550.0, 2.185
    er = math.radians(25.0)
    psi = math.asin(Rk / (Rk + h) * math.cos(er))
    d = math.sqrt((Rk + h) ** 2 - (Rk * math.cos(er)) ** 2) - Rk * math.sin(er)
    snr = (34 + 30 + 10 * math.log10(20) + 12 * math.log10(math.cos(psi)) - (32.45 + 20 * math.log10(f)
           + 20 * math.log10(d * 1e3)) - 6.22 / f ** 1.5 - 3 - (-174 + 10 * math.log10(20e6) + 7))
    out["LinkMarginEmin"] = f"{snr + 8:.1f}"
    if mc:
        out["RateMcRuns"] = f"{40}"
        errs_near = [100 * abs(an[k]["Hnear"] - r["mc_Hnear"]) / r["mc_Hnear"] for k, r in mcd.items() if k in an]
        out["RateNearErrMax"] = f"{max(errs_near):.0f}"
        dn = [100 * abs(an[k]["Hnear"] - r["mc_Hnear"]) / r["mc_Hnear"] for k, r in mcd.items()
              if k in an and an[k]["Nvis"] >= DENSE_NBAR and k[1] == 25]
        out["RateNearErrDense"] = f"{max(dn):.0f}"
        if (396, 25, 550) in mcd:
            r = mcd[(396, 25, 550)]
            out["RateNearErrSparse"] = f"{100 * abs(an[(396, 25, 550)]['Hnear'] - r['mc_Hnear']) / r['mc_Hnear']:.0f}"
        dense = [100 * (an[k]["Hc"] - r["mc_Hexit"]) / r["mc_Hexit"] for k, r in mcd.items()
                 if k in an and an[k]["Nvis"] >= DENSE_NBAR]
        out["RateExitErrMin"] = f"{min(dense):.0f}"
        out["RateExitErrMax"] = f"{max(dense):.0f}"
        allerr = [100 * (an[k]["Hc"] - r["mc_Hexit"]) / r["mc_Hexit"] for k, r in mcd.items() if k in an]
        out["RateExitErrAll"] = f"{max(abs(x) for x in allerr):.0f}"
        if (1584, 25, 550) in mcd:
            out["RateExitMcRef"] = f"{mcd[(1584, 25, 550)]['mc_Hexit']:.3f}"
        if (1584, 40, 550) in mcd:
            out["RateFortyMcRef"] = f"{mcd[(1584, 40, 550)]['mc_Hexit']:.3f}"
        if (792, 25, 550) in mcd and (6336, 25, 550) in mcd:
            out["RateNearGrowth"] = f"{mcd[(6336, 25, 550)]['mc_Hnear'] / mcd[(792, 25, 550)]['mc_Hnear']:.1f}"
            out["RateExitDrop"] = f"{100 * (1 - mcd[(6336, 25, 550)]['mc_Hexit'] / mcd[(792, 25, 550)]['mc_Hexit']):.0f}"
    ph = os.path.join(HERE, "phasing.json")
    if os.path.exists(ph):
        P = json.load(open(ph))
        errs = [100 * (v["an"] - v["mean"]) / v["mean"] for k, v in P.items() if k.isdigit()]
        out["RatePhaseErrMin"] = f"{min(errs):.0f}"
        out["RatePhaseErrMax"] = f"{max(errs):.0f}"
        out["RatePhaseLoRef"] = f"{P['1584']['lo']:.2f}"
        out["RatePhaseHiRef"] = f"{P['1584']['hi']:.2f}"
        out["RatePhaseMeanRef"] = f"{P['1584']['mean']:.2f}"
        out["RatePhaseNF"] = f"{len(P['1584']['F'])}"
        ex = list(P["exit_3168"].values())
        out["RateExitPhaseSpread"] = f"{100 * (max(ex) - min(ex)) / np.mean(ex):.1f}"
    # ns-3 geodesic CHO at three shell sizes (runs_r2/main and density/*/los/geof, ten seeds)
    gn = ns3_geodesic()
    if gn:
        out["RateNsGeoSeven"] = f"{gn[792][0]:.3f}"
        out["RateNsGeoSevenCi"] = f"{gn[792][1]:.3f}"
        out["RateNsGeoThirty"] = f"{gn[3168][0]:.3f}"
        out["RateNsGeoThirtyCi"] = f"{gn[3168][1]:.3f}"
        mcn = {k[0]: r["mc_Hnear"] for k, r in mcd.items() if k[1] == 25 and k[2] == 550}
        ag = [100 * abs(gn[n][0] - mcn[n]) / mcn[n] for n in gn if n in mcn]
        out["RateNsMcGeoMaxErr"] = f"{max(ag):.0f}"
        er = [100 * abs(an[(n, 25, 550)]["Hnear"] - gn[n][0]) / gn[n][0] for n in gn]
        out["RateNsGeoErrMin"] = f"{min(er):.0f}"
        out["RateNsGeoErrMax"] = f"{max(er):.0f}"
    with open(os.path.join(HERE, "macros_rate.tex"), "w") as fo:
        fo.write("% written by sim/analysis_rate/make_rate_outputs.py; do not edit\n")
        for k in sorted(out):
            fo.write(f"\\newcommand{{\\{k}}}{{{out[k]}}}\n")
    json.dump({f"{k[0]}_{k[1]}_{k[2]}": {kk: (float(vv) if isinstance(vv, (float, np.floating)) else vv)
                                          for kk, vv in v.items()} for k, v in an.items()},
              open(os.path.join(HERE, "analysis_rate.json"), "w"), indent=1)
    fig(an, mcd, ns)
    print(json.dumps(out, indent=1))


def fig(an, mcd, ns):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use(os.path.abspath(os.path.join(PAPER, "../templates/ieee.mplstyle")))
    plt.rcParams.update({"font.family": "serif", "pdf.fonttype": 42, "ps.fonttype": 42, "legend.fontsize": 6.5,
                         "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})
    BLUE, ORANGE, GREEN, INK = "#2a78d6", "#eb6834", "#1baf7a", "#1f2933"
    fig, ax = plt.subplots(figsize=(3.45, 2.6))
    Ns = np.geomspace(396, 6336, 13)
    for eQ, col, lab in ((25, BLUE, r"exit rule, $e_Q=25^\circ$"), (40, ORANGE, r"exit rule, $e_Q=40^\circ$")):
        Hc = [60 / M.poisson_chain(int(n), 550, 53, eQ, K=4000) for n in Ns]
        Ha = [60 / M.asym(int(n), 550, 53, eQ) if M.rates(int(n), 550, 53, eQ)["Nvis"] >= 2.0 else np.nan for n in Ns]
        ax.plot(Ns, Hc, color=col, lw=1.1, label=lab)
        ax.plot(Ns, Ha, color=col, lw=0.9, ls="--")
        pts = sorted((k[0], r) for k, r in mcd.items() if k[1] == eQ and k[2] == 550)
        if pts:
            x = [p[0] for p in pts]
            y = [p[1]["mc_Hexit"] for p in pts]
            lo = [p[1]["mc_Hexit"] - p[1]["mc_Hexit_lo"] for p in pts]
            hi = [p[1]["mc_Hexit_hi"] - p[1]["mc_Hexit"] for p in pts]
            ax.errorbar(x, y, yerr=[lo, hi], fmt="o", ms=3.5, mfc="white", color=col, capsize=1.5, lw=0.8)
    Hn = [M.rates(int(n), 550, 53, 25)["Hnear"] for n in Ns]
    ax.plot(Ns, Hn, color=GREEN, lw=1.1, ls=":", label="nearest-satellite rule")
    pts = sorted((k[0], r) for k, r in mcd.items() if k[1] == 25 and k[2] == 550)
    if pts:
        ax.errorbar([p[0] for p in pts], [p[1]["mc_Hnear"] for p in pts],
                    yerr=[p[1]["mc_Hnear_ci"] for p in pts], fmt="s", ms=3.2, mfc="white", color=GREEN,
                    capsize=1.5, lw=0.8)
    # ns-3 (ten seeds)
    ax.plot([1584], [float(ns["MainLosMstHo"])], "o", color=BLUE, ms=4.5)
    ax.plot([792, 1584, 3168], [float(ns["DensHalfLosMstfortyHo"]), float(ns["MainLosMstfortyHo"]),
                                float(ns["DensDoubleLosMstfortyHo"])], "o", color=ORANGE, ms=4.5)
    gn = ns3_geodesic()
    if gn:
        ax.errorbar(sorted(gn), [gn[n][0] for n in sorted(gn)], yerr=[gn[n][1] for n in sorted(gn)],
                    fmt="s", color=GREEN, ms=4.5, capsize=1.5, lw=0.8)
    ax.plot([], [], color=INK, lw=0.9, ls="--", label=r"dense limit, $\bar{N}\geq 2$")
    ax.plot([], [], "o", mfc="white", color=INK, ms=3.5, label="Walker Monte Carlo (circle: exit, square: nearest)")
    ax.plot([], [], "o", color=INK, ms=4.5, label="ns-3 (MST, MST-40, geodesic)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([396, 792, 1584, 3168, 6336])
    ax.set_xticklabels(["396", "792", "1584", "3168", "6336"])
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    ax.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.0]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels(["0.2", "0.3", "0.4", "0.6", "1", "1.5", "2"])
    ax.set_ylim(0.18, 2.1)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("satellites in the shell, $N_s$")
    ax.set_ylabel("handovers per terminal-minute")
    ax.legend(loc="upper left", frameon=False, ncol=1, fontsize=6)
    fig.tight_layout(pad=0.3)
    os.makedirs(os.path.join(PAPER, "figures_v7"), exist_ok=True)
    fig.savefig(os.path.join(PAPER, "figures_v7", "fig_rate.pdf"))


if __name__ == "__main__":
    main()
