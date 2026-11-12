"""M4: late shares on the real shell with and without the heavy-tail satellites
(|along-track offset at the start| > 20 km, about 6 sigma), with a cluster bootstrap over
(seed, satellite) pairs (gate M8)."""
import csv, os, sys, random, json
import numpy as np
SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def runs(tag, base="runs_r2"):
    d = os.path.join(SIM, base, tag)
    return sorted(os.path.join(d, x) for x in os.listdir(d) if x.startswith("seed")) if os.path.isdir(d) else []
def analyse(tag, base="runs_r2", thr=20e3, B=2000):
    ev = []  # (seed, sat, late, big)
    for d in runs(tag, base):
        off = {int(o["sat"]): abs(float(o["d0_a_m"])) for o in csv.DictReader(open(os.path.join(d, "eph_offsets.csv")))}
        s = os.path.basename(d)
        for e in csv.DictReader(open(os.path.join(d, "events.csv"))):
            if e["kind"] != "ho" or e["true_serv_exit_s"] in ("", None): continue
            sat = int(e["from_sat"])
            ev.append((s, sat, float(e["true_serv_exit_s"]) <= 0, off.get(sat, 0) > thr))
    if not ev: return None
    late = [x for x in ev if x[2]]
    cl = {}
    for x in ev: cl.setdefault((x[0], x[1]), []).append(x[2])
    keys = list(cl); rng = random.Random(3); bs = []
    for _ in range(B):
        ks = [keys[rng.randrange(len(keys))] for _ in keys]
        n = sum(len(cl[k]) for k in ks); l = sum(sum(cl[k]) for k in ks); bs.append(100 * l / n)
    lp = {}
    for x in late: lp[(x[0], x[1])] = lp.get((x[0], x[1]), 0) + 1
    top = sorted(lp.values(), reverse=True)
    small = [x for x in ev if not x[3]]
    return dict(n=len(ev), late=100 * len(late) / len(ev), ci=(float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))),
                pairs=len(lp), sats=len({k[1] for k in lp}), top5=sum(top[:5]) / max(1, len(late)),
                late_big=sum(1 for x in late if x[3]), late_wo=100 * sum(1 for x in small if x[2]) / max(1, len(small)),
                share_big_ev=sum(1 for x in ev if x[3]) / len(ev))
out = {}
for tag in sys.argv[1:]:
    r = analyse(tag)
    # post hoc threshold sensitivity (10 and 30 km)
    r["wo10"] = analyse(tag, thr=10e3, B=10)["late_wo"]
    r["wo30"] = analyse(tag, thr=30e3, B=10)["late_wo"]
    out[tag] = r
    print(tag, json.dumps(r))
json.dump(out, open(os.path.join(SIM, "analysis_nov", "late_clusters.json"), "w"), indent=1)
