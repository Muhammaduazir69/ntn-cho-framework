"""Latitude sweep: exact-geometry Monte Carlo of the 1584/72/17 shell at
550 km, 53 deg, for the exit rule (e_Q = 25 deg) and the nearest-satellite rule, read with
estimator E1 at 1200 s, against the analysis (window-free closed form and eq. near)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis_rate"))
import walker_mc as W
import rate_model as RM
import rate_est as E
out = {}
for lat in (0.0, 15.0, 30.0, 45.0, 50.0):
    rng = np.random.default_rng(31)
    X, D, Xn, Dn = [], [], [], []
    g = []
    for k in range(60):
        ho, res = W.run(1584, 72, 17, 550, 53, 25, 1200.0, rng, lat, dt=0.5)
        x, d = E.intervals([list(np.cumsum(res))], 1200.0); X.append(x); D.append(d)
    ex = [60.0 / E.km_mean(sum(X[i:i+10], []), sum(D[i:i+10], [])) for i in range(0, 60, 10)]
    rng = np.random.default_rng(37)
    for k in range(60):
        raan, u0, n, r, i = W.shell(1584, 72, 17, 550, 53)
        u0 = u0 + n * rng.uniform(0, 86400)
        t = np.arange(0, 1200.0, 1.0)
        Ev = W.elev(t, raan, u0, n, r, i, np.radians(lat), rng.uniform(0, 2 * np.pi))
        best = np.argmax(Ev, axis=0)
        ch = t[1:][np.diff(best) != 0]
        x, d = E.intervals([list(ch)], 1200.0); Xn.append(x); Dn.append(d)
    ne = [60.0 / E.km_mean(sum(Xn[i:i+10], []), sum(Dn[i:i+10], [])) for i in range(0, 60, 10)]
    a = RM.rates(1584, 550, 53, 25, lat)
    out[str(lat)] = dict(mc_exit=float(np.mean(ex)), mc_exit_sd=float(np.std(ex, ddof=1)), mc_near=float(np.mean(ne)),
                         mc_near_sd=float(np.std(ne, ddof=1)), an_exit=float(a["Hp"]), an_near=float(a["Hnear"]),
                         nbar=float(a["Nvis"]))
    print(lat, {k: round(v, 3) for k, v in out[str(lat)].items()}, flush=True)
json.dump(out, open(os.path.join(HERE, "mc_lat.json"), "w"), indent=1)
