#!/usr/bin/env python3
"""Margin analysis of the November 2026 campaign: the outage model of Eq. (20)
evaluated exactly as printed,

    J(z) = J0 + nu0 (1 + k sbar z) (c_h + ct q(z)),  ct = c_l + beta c_h,

with q(z) the share of exits whose prediction error exceeds z sigma. Exits are planned
decisions plus escapes taken after the true exit (the latter are exits the margin missed;
their error is at least the predicted time left at the escape). Every constant is listed:
J0 (outage not caused by exits, mostly blockage), c_h (measured mean interruption),
c_l (fitted), beta (fitted), nu0 and k (fitted). Functions take a list of runs, so the same
code fits on one seed set and predicts another."""
import csv, json, math, os, sys, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
ZS = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 6]
TX = 0.2


def zdir(ch, z, base="runs_r2"):
    if z == 0:
        return os.path.join(SIM, base, "zsweep", ch, "z0")
    if z == 2:
        return os.path.join(SIM, base, "realtle", ch, "tte")
    return os.path.join(SIM, base, "zsweep", ch, f"z{z:g}")


RE_KM, MU_KM, OM = 6378.137, 398600.4418, 7.2921159e-5


def eta_v(h_km=470.0, lat=30.0, inc=53.0):
    rs = RE_KM + h_km
    vg = math.sqrt(MU_KM / rs) * RE_KM / rs
    la = math.radians(lat)
    sa = min(1.0, math.cos(math.radians(inc)) / math.cos(la))
    eta = vg / math.hypot(vg * math.sqrt(1 - sa * sa), vg * sa - OM * RE_KM * math.cos(la))
    return eta, math.sqrt(MU_KM / rs) * 1e3


ETA, VSAT = eta_v()


def load_run(d, z):
    S = next(csv.DictReader(open(os.path.join(d, "summary.csv"))))
    U, T = float(S["num_ues"]), float(S["sim_time_s"])
    off = {int(o["sat"]): o for o in csv.DictReader(open(os.path.join(d, "eph_offsets.csv")))}
    ev = [e for e in csv.DictReader(open(os.path.join(d, "events.csv"))) if e["kind"] == "ho"]
    n_ho = len(ev)
    exits, late, eps, sig, intr = 0, 0, [], [], []
    for e in ev:
        tr = e["trigger"]
        true = e["true_serv_exit_s"]
        if true in ("", None):
            continue
        true = float(true)
        is_exit = tr == "planned" or (tr == "escape" and true <= 0.0)
        if e["interruption_ms"] not in ("", None):
            intr.append(float(e["interruption_ms"]) / 1e3)
        if not is_exit:
            continue
        exits += 1
        if true <= 0.0:
            late += 1
        # exit error from the element-set offsets via Lemma 2 (an input, not an outcome):
        # along-track offset of the serving satellite at the decision, interpolated over the run
        o = off.get(int(e["from_sat"]))
        if o is not None:
            t = float(e["t_s"])
            da = float(o["d0_a_m"]) + (float(o["dT_a_m"]) - float(o["d0_a_m"])) * t / T
            eps.append(-ETA * da / VSAT)
        if z > 0 and tr == "planned":
            sig.append((float(e["margin_s"]) - TX) / z)
    return dict(z=z, U=U, T=T, n_ho=n_ho, exits=exits, late=late, eps=eps, sig=sig,
                intr=intr, out=float(S["outage_frac"]) * 100.0,
                out_s=float(S["outage_frac"]) * T * U, seed=int(S["run"]))


def load(ch, seeds, zs=ZS, base="runs_r2"):
    R = []
    for z in zs:
        for s in seeds:
            d = os.path.join(zdir(ch, z, base), f"seed{s:02d}")
            if os.path.exists(os.path.join(d, "summary.csv")):
                R.append(load_run(d, z))
    return R


def fit(R, sbar=None):
    c_h = float(np.mean([x for r in R for x in r["intr"]]))
    if sbar is None:
        sbar = float(np.median([x for r in R for x in r["sig"]]))
    # outage seconds = J0 U T + c_h exits + ct late; J0 holds every outage that is not an
    # exit (blockage escapes and their interruptions), ct the full cost of a late exit
    # within-seed fit (seed fixed effects): the seed fixes the drop, the blockage process and the
    # element-set realization, so ct is estimated from how outage moves with late exits across z
    Y = np.array([r["out_s"] - c_h * r["exits"] for r in R])
    Lt = np.array([r["late"] for r in R], float)
    sd = np.array([r["seed"] for r in R])
    Yd, Ld = Y.copy(), Lt.copy()
    for s_ in set(sd):
        m = sd == s_
        Yd[m] -= Y[m].mean()
        Ld[m] -= Lt[m].mean()
    ct = float(np.dot(Ld, Yd) / np.dot(Ld, Ld))
    J0 = float(np.mean([(r["out_s"] - c_h * r["exits"] - ct * r["late"]) / (r["U"] * r["T"]) for r in R]))
    # exit rate per terminal-second = nu0 (1 + k sbar z)
    zz = np.array([r["z"] for r in R], float)
    er = np.array([r["exits"] / (r["U"] * r["T"]) for r in R])
    b1, b0 = np.polyfit(zz * sbar, er, 1)
    k_fit = float(b1 / b0)
    # nu0: exit rate at z = 0; k is not fitted but derived: a margin m shortens every residence
    # by m - T_x, so nu(m) = 1 / (E[R] - (m - T_x)) ~ nu0 (1 + nu0 (m - T_x)), i.e. k = nu0
    z0 = [r for r in R if r["z"] == 0]
    nu0 = float(np.mean([r["exits"] / (r["U"] * r["T"]) for r in (z0 or R)]))
    k = nu0
    # extra handovers per late exit: n_ho - exits = esc0 U T + beta late
    X2 = np.array([[r["U"] * r["T"], r["late"]] for r in R])
    Y2 = np.array([r["n_ho"] - r["exits"] for r in R], float)
    (esc0, beta), *_ = np.linalg.lstsq(X2, Y2, rcond=None)
    eps = np.array([e for r in R for e in r["eps"]])
    return dict(c_h=c_h, J0=float(J0), c_l=float(ct - beta * c_h), nu0=nu0, k=k, k_fit=k_fit, beta=float(beta),
                ct=float(ct), sbar=sbar, eps=eps)


TS = 0.1  # decision period: a planned exit is decided at the first epoch with tau_hat <= m,
          # on average T_s/2 below m, and is late when the exit error exceeds that (A4)


def q(P, z, eps=None):
    e = P["eps"] if eps is None else eps
    return float(np.mean(e > TX + z * P["sbar"] - TS / 2)) if len(e) else float("nan")


def J(P, z, eps=None):
    """outage share in %"""
    return 100.0 * (P["J0"] + P["nu0"] * (1 + P["k"] * P["sbar"] * z) * (P["c_h"] + P["ct"] * q(P, z, eps)))


def late_model(P, z, eps=None):
    return 100.0 * q(P, z, eps)


def zstar_numeric(P, eps=None, zmax=8.0):
    zf = np.linspace(0, zmax, 801)
    v = [J(P, z, eps) for z in zf]
    return float(zf[int(np.argmin(v))])


def zstar_gauss(P):
    k, s, ch, ct = P["k"], P["sbar"], P["c_h"], P["ct"]
    if k <= 0:
        return float("nan")
    R = ct / (math.sqrt(2 * math.pi) * k * s * ch)
    c0 = (TX - TS / 2) / s
    return math.sqrt(2 * math.log(R)) - c0 if R > 1 else float("nan")


def zstar_gauss_exact(P):
    """root of g(z) = kt (c_h + ct Q(z + c0)) - ct phi(z + c0)(1 + kt z), kt = nu0 sigma"""
    from scipy import stats as S
    kt, ch, ct = P["nu0"] * P["sbar"], P["c_h"], P["ct"]
    c0 = (TX - TS / 2) / P["sbar"]
    g = lambda z: kt * (ch + ct * S.norm.sf(z + c0)) - ct * S.norm.pdf(z + c0) * (1 + kt * z)
    if g(0) >= 0:
        return 0.0
    lo, hi = 0.0, 20.0
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


def zmeans(R):
    out = {}
    for z in sorted({r["z"] for r in R}):
        v = [r["out"] for r in R if r["z"] == z]
        L = [100.0 * r["late"] / max(1, r["n_ho"]) for r in R if r["z"] == z]
        out[z] = (float(np.mean(v)), float(np.std(v, ddof=1) / math.sqrt(len(v)) * 2.262 if len(v) > 1 else 0), len(v),
                  float(np.mean(L)))
    return out


def cv(ch, folds=((1, 2, 3, 4, 5), (6, 7, 8, 9, 10))):
    """two-fold seed cross-validation: model vs constant on held-out per-z means"""
    res = []
    for a, b in [(folds[0], folds[1]), (folds[1], folds[0])]:
        Ra, Rb = load(ch, a), load(ch, b)
        P = fit(Ra)
        const = float(np.mean([r["out"] for r in Ra]))
        ma = zmeans(Ra)
        step = {z: (ma[0][0] if z == 0 else float(np.mean([ma[x][0] for x in ma if x > 0]))) for x_ in [0] for z in ma}
        mb = zmeans(Rb)
        m0, p0 = mb[0][0], J(P, 0)
        for z, (m, _, _, _) in mb.items():
            res.append((z, m, J(P, z), const, m - m0, J(P, z) - p0, step[z] - step[0]))
    e_model = float(np.mean([(r[1] - r[2]) ** 2 for r in res]))
    e_const = float(np.mean([(r[1] - r[3]) ** 2 for r in res]))
    # profile: change from z = 0 within the held-out fold; the constant model predicts no change
    p_model = float(np.mean([(r[4] - r[5]) ** 2 for r in res if r[0] != 0]))
    p_const = float(np.mean([r[4] ** 2 for r in res if r[0] != 0]))
    p_step = float(np.mean([(r[4] - r[6]) ** 2 for r in res if r[0] != 0]))
    return dict(rmse_model=math.sqrt(e_model), rmse_const=math.sqrt(e_const),
                prof_model=math.sqrt(p_model), prof_const=math.sqrt(p_const), prof_step=math.sqrt(p_step), rows=res)


def boot_k(ch, B=2000, seed=7):
    R = load(ch, range(1, 11))
    rng = random.Random(seed)
    ks, zs = [], []
    for _ in range(B):
        ss = [rng.randint(1, 10) for _ in range(10)]
        Rb = [r for s in ss for r in R if r["seed"] == s]
        P = fit(Rb)
        ks.append(P["k_fit"])
        zs.append(zstar_numeric(P))  # minimizer of J with the offset tail
    ks = np.array(ks)
    return dict(k_lo=float(np.percentile(ks, 2.5)), k_hi=float(np.percentile(ks, 97.5)),
                p_nonpos=float(np.mean(ks <= 0)), z_lo=float(np.percentile(zs, 2.5)),
                z_hi=float(np.percentile(zs, 97.5)))


if __name__ == "__main__":
    out = {}
    for ch in ["los", "suburban"]:
        R = load(ch, range(1, 11))
        P = fit(R)
        zm = zmeans(R)
        rows = {z: dict(sim=zm[z][0], ci=zm[z][1], late_sim=zm[z][3], model=J(P, z), late_model=late_model(P, z) *
                        (sum(r["exits"] for r in R) / sum(r["n_ho"] for r in R))) for z in zm}
        c = cv(ch)
        b = boot_k(ch)
        out[ch] = dict(params={k: v for k, v in P.items() if k != "eps"}, zstar_gauss=zstar_gauss(P),
                       zstar_num=zstar_numeric(P), rows=rows, cv=dict(rmse_model=c["rmse_model"],
                       rmse_const=c["rmse_const"], prof_model=c["prof_model"], prof_const=c["prof_const"], prof_step=c["prof_step"]), boot=b, n_exits=int(sum(r["exits"] for r in R)),
                       q6=q(P, 6), share_pos=float(np.mean(P["eps"] > 0)))
        o = out[ch]
        print(ch, {k: round(v, 5) for k, v in o["params"].items()}, "z*g", o["zstar_gauss"], "z*n", o["zstar_num"])
        for z, r in o["rows"].items():
            print(f"  z={z}: sim {r['sim']:.4f}+-{r['ci']:.4f} model {r['model']:.4f}  late {r['late_sim']:.2f} vs {r['late_model']:.2f}")
        print("  cv", o["cv"], "boot", o["boot"], "q6", o["q6"], "pos", o["share_pos"])
    json.dump(out, open(os.path.join(HERE, "margin_model.json"), "w"), indent=1, default=float)
