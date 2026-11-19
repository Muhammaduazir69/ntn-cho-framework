"""Element-set error statistics of the epoch and held-out (+12 h) windows, and the margin
model's late share with the held-out exit-error sign reversed."""
import csv, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import margin_model as MM
def offs(d):
    return np.array([float(o["d0_a_m"]) for o in csv.DictReader(open(os.path.join(d, "eph_offsets.csv")))])
out = {}
for nm, base, seeds in (("epoch", os.path.join(SIM, "runs_r2/realtle/los/mst"), range(1, 11)),
                        ("held", os.path.join(SIM, "runs_nov/oos/los/mst"), range(11, 21))):
    p68, mean, med, big, sd = [], [], [], [], []
    allv = []
    for s in seeds:
        a = offs(os.path.join(base, f"seed{s:02d}")) / 1e3
        a = a[np.abs(a) < 100]
        p68.append(np.percentile(np.abs(a), 68.27)); allv += list(a)
    allv = np.array(allv)
    z = 0 if nm == "epoch" else 0
    R = [MM.load_run(os.path.join(base, f"seed{s:02d}"), 0) for s in seeds]
    eps = np.array([e for r in R for e in r["eps"]])
    out[nm] = dict(p68_lo=min(p68), p68_hi=max(p68), mean=float(allv.mean()), median=float(np.median(allv)),
                   big=float(np.mean(np.abs(allv) > 20)), lateside=float(np.mean(eps > 0)),
                   eps68=float(np.percentile(np.abs(eps), 68.27)), n_exits=len(eps))
# model with epoch constants, held-out offsets, sign as flown and reversed
P = MM.fit(MM.load("los", range(1, 11)))
Ro = []
for pol, z in (("mst", 0), ("tte", 2)):
    for s in range(11, 21):
        Ro.append(MM.load_run(os.path.join(SIM, "runs_nov/oos/los", pol, f"seed{s:02d}"), z))
e = np.array([x for r in Ro for x in r["eps"]])
em = np.array([x for r in Ro if r["z"] == 2 for x in r["eps"]])
for k, v in (("flown", e), ("mirrored", -e)):
    out[k] = {f"z{z}": 100 * MM.q(P, z, v) for z in (0, 2, 6)}
    zz = np.linspace(0, 20, 401)
    out[k]["z_for_0p1"] = float(next((z for z in zz if 100 * MM.q(P, z, v) <= 0.1), float("nan")))
out["margin_arm_only"] = {f"z{z}": 100 * MM.q(P, z, em) for z in (0, 2)}
json.dump(out, open(os.path.join(HERE, "window_stats.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
