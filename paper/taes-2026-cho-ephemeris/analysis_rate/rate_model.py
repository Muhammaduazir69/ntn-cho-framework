#!/usr/bin/env python3
"""Constellation-level handover-rate analysis (Section IV-A of the TAES manuscript) and its
validation against an exact-geometry Monte Carlo of the Walker shell (walker_mc.py).

Analysis (Lemma 1, Propositions 1-2):
  rho_Q   ground range of the contour e_Q;  T(y) = (2R/v) arccos(cos(rho/R)/cos(y/R))
  sigma   density of sub-satellite points at latitude phi
  lambda  = 2 rho sigma v_rel, pass arrival rate;  Lambda(r) = lambda int_r^inf Fbar_T(s) ds
  exit rule, random shell (Poisson):  P(R <= r) = exp(-Lambda(r))
  exit rule, Walker shell (lattice):  P(R <= r) ~= max(0, 1 - Lambda(r))
  nearest-satellite rule:             H = 4 sqrt(sigma) v_rel / pi
Writes analysis_rate.csv, mc_rate.csv and macros_rate.tex here.

    python3 rate_model.py [--mc RUNS]
"""
import argparse
import os
import numpy as np
from scipy import integrate

import walker_mc as W

HERE = os.path.dirname(os.path.abspath(__file__))
R, MU, WE = W.R, W.MU, W.WE


def geom(N, h, inc, eQ, lat=30.0):
    r = R + h
    v = np.sqrt(MU / r)
    psi = np.arcsin(R / r * np.cos(np.radians(eQ)))
    rho = R * (np.pi / 2 - np.radians(eQ) - psi)
    i, p = np.radians(inc), np.radians(lat)
    vg = v * R / r
    cb = np.cos(i) / np.cos(p)
    sb = np.sqrt(1 - cb ** 2)
    vrel = np.hypot(vg * cb - WE * R * np.cos(p), vg * sb)
    sig = N / (2 * np.pi ** 2 * R ** 2) / np.sqrt(np.sin(i) ** 2 - np.sin(p) ** 2)
    return dict(rho=rho, vrel=vrel, sig=sig, lam=2 * rho * sig * vrel, Tmax=2 * rho / vrel)


def rates(N, h, inc, eQ, lat=30.0):
    g = geom(N, h, inc, eQ, lat)
    rho, vrel, Tm = g["rho"], g["vrel"], g["Tmax"]

    def Fb(s):
        if s >= Tm:
            return 0.0
        if s <= 0:
            return 1.0
        return R * np.arccos(min(1.0, np.cos(rho / R) / np.cos(s * vrel / (2 * R)))) / rho

    L = lambda r: g["lam"] * integrate.quad(Fb, r, Tm, limit=200)[0]
    ET = integrate.quad(Fb, 0, Tm, limit=200)[0]
    Ep = integrate.quad(lambda r: 1 - np.exp(-L(r)), 0, Tm, limit=200)[0]
    El = integrate.quad(lambda r: min(1.0, L(r)), 0, Tm, limit=200)[0]
    return dict(g, ET=ET, Nvis=g["lam"] * ET, ERp=Ep, ERl=El, Hp=60 / Ep, Hl=60 / El,
                Hnear=60 * 4 * np.sqrt(g["sig"]) * vrel / np.pi, Hinf=60 * vrel / (2 * rho))


def poisson_chain(N, h, inc, eQ, K=40000, seed=5):
    """Residence chain of Proposition 1 under Poisson passes, evaluated by sampling its kernel:
    given the previous residence w, the candidates are the passes that arrived in the last w
    seconds; if none is still inside the contour the serving cell is kept until the next pass."""
    g = geom(N, h, inc, eQ)
    rho, v, lam, Tm = g["rho"], g["vrel"], g["lam"], g["Tmax"]
    rng = np.random.default_rng(seed)
    T = lambda y: 2 * R / v * np.arccos(np.minimum(1.0, np.cos(rho / R) / np.cos(y / R)))
    w, res = Tm / 2, []
    for _ in range(K):
        n = rng.poisson(lam * w)
        rem = T(rng.uniform(0, rho, n)) - rng.uniform(0, w, n)
        if n and rem.max() > 0:
            r = rem.max()
        else:
            if res:
                res[-1] += rng.exponential(1 / lam)
            r = T(rng.uniform(0, rho))
        res.append(r)
        w = r
    return float(np.mean(res[200:]))


def asym(N, h, inc, eQ):
    """Dense-regime closed form, E[R] = Tmax (1 - Gamma(5/3) (8 sqrt2 Nbar / (3 pi))^(-2/3))."""
    from math import gamma
    a = rates(N, h, inc, eQ)
    return a["Tmax"] * (1 - gamma(5 / 3) * (8 * np.sqrt(2) * a["Nvis"] / (3 * np.pi)) ** (-2 / 3))


def mc_nearest(N, P, F, h, inc, eQ, H, rng, lat=30.0, dt=1.0):
    raan, u0, n, r, i = W.shell(N, P, F, h, inc)
    u0 = u0 + n * rng.uniform(0, 86400)
    t = np.arange(0, H, dt)
    E = W.elev(t, raan, u0, n, r, i, np.radians(lat), rng.uniform(0, 2 * np.pi))
    best = np.argmax(E, axis=0)
    return int(np.count_nonzero(np.diff(best)))


def mc(N, P, F, h, inc, eQ, runs, H=2400.0, lat=30.0):
    rng = np.random.default_rng(7)
    res, hn = [], []
    for _ in range(runs):
        _, rr = W.run(N, P, F, h, inc, eQ, H, rng, lat, dt=1.0)
        res += rr[1:]
        hn.append(mc_nearest(N, P, F, h, inc, eQ, H, rng, lat))
    res = np.array(res)
    m, se = res.mean(), res.std(ddof=1) / np.sqrt(len(res))
    hn = 60 * np.array(hn) / H
    return dict(ER=m, Hexit=60 / m, Hexit_lo=60 / (m + 1.96 * se), Hexit_hi=60 / (m - 1.96 * se),
                Hnear=hn.mean(), Hnear_ci=1.96 * hn.std(ddof=1) / np.sqrt(runs), nres=len(res))


CASES = [  # N, P, F, h, inc, eQ
    (396, 36, 17, 550, 53, 25), (792, 36, 17, 550, 53, 25), (1584, 72, 17, 550, 53, 25),
    (3168, 72, 17, 550, 53, 25), (6336, 144, 17, 550, 53, 25),
    (792, 36, 17, 550, 53, 40), (1584, 72, 17, 550, 53, 40), (3168, 72, 17, 550, 53, 40),
    (6336, 144, 17, 550, 53, 40), (648, 18, 1, 1200, 87.9, 25),
]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", type=int, default=0)
    ap.add_argument("--skip", type=int, default=0, help="skip the first SKIP cases (resume)")
    a = ap.parse_args()
    rows = []
    for c in CASES[a.skip:]:
        an = rates(c[0], c[3], c[4], c[5])
        row = dict(N=c[0], P=c[1], h=c[3], inc=c[4], eQ=c[5], **{k: an[k] for k in
                   ("rho", "vrel", "sig", "lam", "Tmax", "ET", "Nvis", "ERp", "Hp", "Hnear", "Hinf")})
        row["ERc"] = poisson_chain(c[0], c[3], c[4], c[5])
        row["Hc"] = 60 / row["ERc"]
        row["ERa"] = asym(c[0], c[3], c[4], c[5])
        row["Ha"] = 60 / row["ERa"]
        if a.mc:
            row.update({"mc_" + k: v for k, v in mc(*c, runs=a.mc).items()})
        rows.append(row)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    import csv
    with open(os.path.join(HERE, "mc_rate.csv" if a.mc else "analysis_rate.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
