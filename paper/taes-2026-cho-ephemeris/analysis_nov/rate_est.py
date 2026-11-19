"""Rate estimators shared by ns-3 and the exact-geometry Monte Carlo (protocol E1, E2).
E1: the intervals between consecutive completed handovers of a terminal after its first,
the last one right-censored at the horizon; Kaplan-Meier restricted mean pooled over the
terminals of a run; rate = 60 / mean (per terminal-minute).
E2: completed handovers per terminal-minute over [0, H] from attach."""
import csv, os, sys
import numpy as np
SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SIM, "analysis_rate"))


def km_mean(x, d):
    """restricted mean of a right-censored sample (x durations, d 1 = event)"""
    x, d = np.asarray(x, float), np.asarray(d, int)
    if len(x) == 0:
        return float("nan")
    o = np.lexsort((-d, x))
    x, d = x[o], d[o]
    n = len(x)
    S, t0, area = 1.0, 0.0, 0.0
    i = 0
    while i < n:
        t = x[i]
        area += S * (t - t0)
        t0 = t
        j = i
        ev = 0
        while j < n and x[j] == t:
            ev += d[j]
            j += 1
        S *= 1.0 - ev / (n - i)
        i = j
    return area


def intervals(times_per_ue, H):
    x, d = [], []
    for ts in times_per_ue:
        ts = sorted(ts)
        if not ts:
            continue
        for a, b in zip(ts[:-1], ts[1:]):
            x.append(b - a)
            d.append(1)
        x.append(H - ts[-1])
        d.append(0)
    return x, d


def ns3_run(d):
    S = next(csv.DictReader(open(os.path.join(d, "summary.csv"))))
    H, U = float(S["sim_time_s"]), int(S["num_ues"])
    per = {u: [] for u in range(U)}
    for e in csv.DictReader(open(os.path.join(d, "events.csv"))):
        if e["kind"] == "ho" and e["success"] == "1":
            per[int(e["ue"])].append(float(e["t_s"]))
    x, dd = intervals(per.values(), H)
    m = km_mean(x, dd)
    e2 = 60.0 * sum(len(v) for v in per.values()) / (U * H)
    return dict(e1=60.0 / m, e2=e2, out=100 * float(S["outage_frac"]), snr=float(S["mean_sinr_db"]),
                ho=float(S["ho_per_ue_min"]), seed=int(S["run"]), H=H)


def mc_exit(N, P, F, h, inc, eQ, H, epochs, seed=11):
    import walker_mc as W
    rng = np.random.default_rng(seed)
    e1, e2 = [], []
    for _ in range(epochs):
        ho, res = W.run(N, P, F, h, inc, eQ, H, rng, 30.0, dt=0.5)
        ts = list(np.cumsum(res))
        x, d = intervals([ts], H)
        e1.append((x, d))
        e2.append(60.0 * ho / H)
    # pool in groups of 10 epochs (100-terminal ns-3 runs pool over terminals; MC epochs are terminals)
    g = [60.0 / km_mean(sum((e1[i][0] for i in range(k, k + 10)), []), sum((e1[i][1] for i in range(k, k + 10)), []))
         for k in range(0, epochs, 10)]
    return dict(e1=float(np.mean(g)), e1_sd=float(np.std(g, ddof=1)), e2=float(np.mean(e2)),
                e2_ci=float(1.96 * np.std(e2, ddof=1) / np.sqrt(epochs)))


def mc_near(N, P, F, h, inc, H, epochs, seed=13, dt=1.0):
    import walker_mc as W
    rng = np.random.default_rng(seed)
    X, D, e2 = [], [], []
    for _ in range(epochs):
        raan, u0, n, r, i = W.shell(N, P, F, h, inc)
        u0 = u0 + n * rng.uniform(0, 86400)
        t = np.arange(0, H, dt)
        E = W.elev(t, raan, u0, n, r, i, np.radians(30.0), rng.uniform(0, 2 * np.pi))
        best = np.argmax(E, axis=0)
        ch = t[1:][np.diff(best) != 0]
        x, d = intervals([list(ch)], H)
        X.append(x); D.append(d)
        e2.append(60.0 * len(ch) / H)
    g = [60.0 / km_mean(sum(X[k:k + 10], []), sum(D[k:k + 10], [])) for k in range(0, epochs, 10)]
    return dict(e1=float(np.mean(g)), e2=float(np.mean(e2)), e2_ci=float(1.96 * np.std(e2, ddof=1) / np.sqrt(epochs)))


if __name__ == "__main__":
    import json
    out = {}
    for H in (600.0, 1200.0, 2400.0):
        for eQ in (25, 40):
            out[f"exit_e{eQ}_H{int(H)}"] = mc_exit(1584, 72, 17, 550, 53, eQ, H, 120)
            print(f"exit_e{eQ}_H{int(H)}", out[f"exit_e{eQ}_H{int(H)}"], flush=True)
        out[f"near_H{int(H)}"] = mc_near(1584, 72, 17, 550, 53, H, 120)
        print(f"near_H{int(H)}", out[f"near_H{int(H)}"], flush=True)
    json.dump(out, open(os.path.join(SIM, "analysis_nov", "mc_rate_like.json"), "w"), indent=1)
