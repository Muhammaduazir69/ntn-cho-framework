#!/usr/bin/env python3
"""Independent verification of the Paper 2 round-2 results.

Reads the raw run outputs under runs_r2/ with its own code (it does not import
analyze.py) and checks
  A. run integrity: every run returned 0, has its outputs, seed counts as
     designed, and only the configurations run_campaign.py defines exist
     (retired interim policies would be flagged, not counted);
  B. per-run invariants: counters consistent, rate formula, outage covers the
     booked handover interruptions, ping-pong <= handovers, trigger column
     consistent with the summary's escape share;
  C. identities the design implies (common random numbers and the switches of
     each mechanism), e.g. MST + margin == MST at zero ephemeris error;
  D. headline macros in ../results_macros.tex against an independent
     recomputation (tolerance = rounding of the macro);
  E. paired statistics: exact Wilcoxon p (scipy) and the Holm adjustment of the
     headline families, against the printed macros.
Exit code 0 only if every check passes. Writes verify_report.txt.
"""
import csv
import glob
import json
import math
import os
import re
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, os.environ.get("RUNS_DIR", "runs_r2"))
MACROS = os.path.join(HERE, "..", "results_macros.tex")
out = []
fails = 0


def check(ok, msg):
    global fails
    out.append(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails += 1


def summ(tag):
    rows = {}
    for d in sorted(glob.glob(os.path.join(RUNS, tag, "seed*"))):
        f = os.path.join(d, "summary.csv")
        if os.path.exists(f):
            r = next(csv.DictReader(open(f)))
            r["_dir"] = d
            rows[int(r["run"])] = r
    return rows


def col(tag, metric):
    """Independent metric definitions (not analyze.py's)."""
    o = {}
    for k, r in summ(tag).items():
        ue_h = float(r["num_ues"]) * float(r["sim_time_s"]) / 3600.0
        ok = float(r["ho_success"])
        v = {"ho": ok / (float(r["num_ues"]) * float(r["sim_time_s"]) / 60.0),
             "out": 100.0 * float(r["outage_frac"]),
             "hof": float(r["hof"]) / ue_h,
             "rlf": float(r["rlf"]) / ue_h,
             "ppone": 100.0 * float(r["pp1"]) / max(1.0, ok),
             "snr": float(r["mean_sinr_db"])}[metric]
        o[k] = v
    return o


def mean(tag, metric):
    v = list(col(tag, metric).values())
    return float(np.mean(v)) if v else float("nan")


# ---------------------------------------------------------------- A. integrity
sys.path.insert(0, HERE)
os.environ.setdefault("RUNS_DIR", os.path.basename(RUNS))
import run_campaign as RC  # noqa: E402  (only for the list of designed configurations)

designed = set()
seeds_per = {"pareto": 5, "trace": 1}
for name, f in RC.EXPS.items():
    for tag, s, _ in f(list(range(1, seeds_per.get(name, 10) + 1))):
        designed.add(tag)
bad_rc, missing = [], []
for cj in glob.glob(os.path.join(RUNS, "**", "cmd.json"), recursive=True):
    j = json.load(open(cj))
    d = os.path.dirname(cj)
    if j.get("rc", 1) != 0:
        bad_rc.append(d)
    acc = any("--mode=tte-accuracy" in a for a in j["cmd"])
    if not os.path.exists(os.path.join(d, "tte_accuracy.csv" if acc else "summary.csv")):
        missing.append(d)
check(not bad_rc, f"A1 all runs returned 0 ({len(bad_rc)} non-zero)")
check(not missing, f"A2 every run has its output file ({len(missing)} missing)")
present = {os.path.relpath(os.path.dirname(d), RUNS) for d in glob.glob(os.path.join(RUNS, "**", "seed*"), recursive=True)}
extra = sorted(present - designed)
check(not extra, f"A3 no configuration outside run_campaign.py ({len(extra)} extra: {extra[:4]})")
for exp, n in [("main", 10), ("urban", 10), ("eph", 10), ("ephemp", 10), ("ephblk", 10), ("realtle", 10),
               ("acc", 10), ("lb", 10), ("sens", 10), ("density", 10), ("losmodel", 10), ("reest", 10),
               ("pareto", 5), ("interf", 10), ("ephage", 10)]:
    tags = sorted(t for t in designed if t.startswith(exp + "/"))
    short = [t for t in tags if len(glob.glob(os.path.join(RUNS, t, "seed*", "cmd.json"))) < n]
    check(not short, f"A4 {exp}: {len(tags)} configurations, each with >= {n} seeds ({len(short)} short)")

# ------------------------------------------------------------ B. invariants
nb = 0
viol = {"att": 0, "rate": 0, "pp": 0, "out": 0, "trig": 0}
for sf in glob.glob(os.path.join(RUNS, "**", "summary.csv"), recursive=True):
    r = next(csv.DictReader(open(sf)))
    d = os.path.dirname(sf)
    nb += 1
    att, ok, hof = float(r["ho_attempts"]), float(r["ho_success"]), float(r["hof"])
    viol["att"] += abs(att - ok - hof) > 0.5
    rate = ok / (float(r["num_ues"]) * float(r["sim_time_s"]) / 60.0)
    viol["rate"] += abs(rate - float(r["ho_per_ue_min"])) > 6e-4
    viol["pp"] += float(r["pp5"]) > ok + 0.5 or float(r["pp1"]) > float(r["pp5"]) + 0.5
    ev = os.path.join(d, "events.csv")
    us = os.path.join(d, "ue_summary.csv")
    if os.path.exists(ev) and os.path.exists(us):
        intr, ho, esc = {}, 0, 0
        for e in csv.DictReader(open(ev)):
            if e["kind"] == "ho":
                ho += 1
                esc += e["trigger"] in ("escape", "netesc")
                intr[e["ue"]] = intr.get(e["ue"], 0.0) + float(e["interruption_ms"]) / 1000.0
        for u in csv.DictReader(open(us)):
            if float(u["outage_s"]) + 1e-3 < intr.get(u["ue"], 0.0):
                viol["out"] += 1
                break
        share = esc / ho if ho else 0.0
        viol["trig"] += abs(share - float(r["escape_share"])) > 1e-3
check(viol["att"] == 0, f"B1 attempts = successes + failures in all {nb} runs ({viol['att']} violations)")
check(viol["rate"] == 0, f"B2 handover rate = successes / terminal-minutes ({viol['rate']} violations)")
check(viol["pp"] == 0, f"B3 pp1 <= pp5 <= handovers ({viol['pp']} violations)")
check(viol["out"] == 0, f"B4 per-terminal outage >= booked interruptions ({viol['out']} runs violate)")
check(viol["trig"] == 0, f"B5 escape share = escape triggers in events.csv ({viol['trig']} violations)")


# ------------------------------------------------------- C. design identities
def same(a, b, metric="ho"):
    A, B = col(a, metric), col(b, metric)
    return len(A) > 0 and set(A) == set(B) and all(abs(A[k] - B[k]) < 1e-12 for k in A)


for a, b, why in [("main/los/tte", "main/los/mst", "MST + margin = MST at zero error (LOS)"),
                  ("main/suburban/tte", "main/suburban/mst", "MST + margin = MST at zero error (suburban)"),
                  ("main/los/tte40", "main/los/mst40", "MST-40 + margin = MST-40 at zero error (LOS)"),
                  ("main/los/mst_ute", "main/los/mst", "LOS: U = tau, so U ranking changes nothing")]:
    check(same(a, b), f"C {why}: per seed")
for ch in ["suburban", "urban"]:
    exp = "main" if ch == "suburban" else "urban"
    ca = [[u["class"] for u in csv.DictReader(open(os.path.join(r["_dir"], "ue_summary.csv")))]
          for r in summ(f"{exp}/{ch}/d2f").values()]
    cb = [[u["class"] for u in csv.DictReader(open(os.path.join(r["_dir"], "ue_summary.csv")))]
          for r in summ(f"{exp}/{ch}/ut40").values()]
    check(len(ca) > 0 and ca == cb, f"C {ch}: terminal classes identical across policies, seed by seed")
    # M8: aircraft and ships in line of sight (env 0), trains rural (env 1)
    bad = 0
    for r in summ(f"{exp}/{ch}/ut40").values():
        for u in csv.DictReader(open(os.path.join(r["_dir"], "ue_summary.csv"))):
            want = {"aviation": "0", "maritime": "0", "hst": "1"}.get(u["class"])
            bad += want is not None and u["env"] != want
    check(bad == 0, f"C {ch}: per-class environment as stated ({bad} terminals differ)")

# ----------------------------------------------------------- D. macros
M = {}
for line in open(MACROS):
    m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s*$", line)
    if m:
        M[m.group(1)] = m.group(2)


def num_macro(k):
    v = M.get(k, "").replace("{,}", "").replace("$", "").replace("p=", "").replace("p<", "")
    try:
        return float(v)
    except ValueError:
        return float("nan")


def close(k, val):
    s = M.get(k, "")
    mv = num_macro(k)
    dec = len(s.split(".")[1]) if "." in s else 0
    check(abs(mv - val) <= 0.5 * 10 ** (-dec) + 1e-9, f"D {k} = {s}, independent {val:.4f}")


cut = lambda a, b: 100 * (1 - a / b)
close("MainLosUtHo", mean("main/los/ut", "ho"))
close("MainLosMstHo", mean("main/los/mst", "ho"))
close("MainLosToneHo", mean("main/los/t1", "ho"))
close("CutUtVsMstLos", cut(mean("main/los/ut", "ho"), mean("main/los/mst", "ho")))
close("CutUtVsToneLos", cut(mean("main/los/ut", "ho"), mean("main/los/t1", "ho")))
close("CutMstVsToneLos", cut(mean("main/los/mst", "ho"), mean("main/los/t1", "ho")))
close("MainSubUtfortyOut", mean("main/suburban/ut40", "out"))
for ch, C in [("urban", "Urb"), ("denseurban", "Dense")]:
    for p, P in [("ut40", "Utforty"), ("mst40", "Mstforty"), ("d2f", "Dtwof"), ("t1f40", "Toneforty")]:
        close(f"Urb{C}{P}Out", mean(f"urban/{ch}/{p}", "out"))
        close(f"Urb{C}{P}Hof", mean(f"urban/{ch}/{p}", "hof"))
    close(f"Urb{C}UtfortyPpone", mean(f"urban/{ch}/ut40", "ppone"))
# radial term of eq:sigma, finite difference vs closed form
RE_, h_ = 6378137.0, 550e3


def rho(h, e):
    rs = RE_ + h
    return RE_ * (math.pi / 2 - e - math.asin(RE_ / rs * math.cos(e)))


for edeg in (25.0, 40.0):
    e = math.radians(edeg)
    fd = (rho(h_ + 1.0, e) - rho(h_ - 1.0, e)) / 2.0
    rs = RE_ + h_
    cf = RE_ / rs * math.tan(math.asin(RE_ / rs * math.cos(e)))
    check(abs(fd - cf) < 1e-4, f"D radial derivative at e_Q={edeg:g}: finite difference {fd:.4f}, closed form {cf:.4f}")
# Earth-rotation factor at 30 N, 53 deg inclination, 550 km
mu, om = 3.986004418e14, 7.2921159e-5
vg = math.sqrt(mu / (RE_ + h_)) * RE_ / (RE_ + h_)
saz = math.cos(math.radians(53)) / math.cos(math.radians(30))
rel = math.hypot(vg * math.sqrt(1 - saz ** 2), vg * saz - om * RE_ * math.cos(math.radians(30)))
close("RotKappa", vg / rel)
# accuracy: rms over predicted exits > 0
for key, k in [("sa5000", "AccFiveRmsPos"), ("sa0", "AccZeroRmsPos")]:
    e = []
    for f in glob.glob(os.path.join(RUNS, "acc", key, "seed*", "tte_accuracy.csv")):
        for r in csv.DictReader(open(f)):
            if 0 < float(r["true_tte_s"]) < 890 and float(r["pred_tte_s"]) > 0:
                e.append(float(r["pred_tte_s"]) - float(r["true_tte_s"]))
    if e:
        close(k, float(np.sqrt(np.mean(np.square(e)))))
# late handovers at sigma_A = 5 km
for pol, k in [("mst", "EphMstFiveLate"), ("tte", "EphTteFiveLate"), ("d2", "EphDtwoFiveLate")]:
    lates = []
    for r in summ(f"eph/sa5000/{pol}").values():
        ev = [e for e in csv.DictReader(open(os.path.join(r["_dir"], "events.csv"))) if e["kind"] == "ho"]
        lates.append(100.0 * sum(1 for e in ev if float(e["true_serv_exit_s"]) <= 0) / max(1, len(ev)))
    if lates:
        close(k, float(np.mean(lates)))
# measured TLE error
rows = [r for r in csv.DictReader(open(os.path.join(HERE, "tle", "tle_error_samples.csv"))) if r["dt_min"] == "0"]
E = np.array([[float(r["err_r_m"]), float(r["err_a_m"]), float(r["err_c_m"])] for r in rows])
E = E[np.all(np.abs(E) <= 1e5, axis=1)]
close("TleSigA", float(np.percentile(np.abs(E[:, 1]), 68.27)) / 1000)
close("TleSigC", float(np.percentile(np.abs(E[:, 2]), 68.27)))


# ----------------------------------------------------------- E. statistics
def wil(a, b):
    """Exact two-sided signed-rank p by enumerating all sign flips, scipy mid-ranks."""
    import itertools
    k = sorted(set(a) & set(b))
    d = np.array([round(a[i] - b[i], 9) for i in k])
    d = d[d != 0]
    if len(d) == 0:
        return 1.0
    r = stats.rankdata(np.abs(d))
    T = r.sum()
    obs = min(r[d > 0].sum(), T - r[d > 0].sum())
    cnt = 0
    for s in itertools.product([0, 1], repeat=len(d)):
        v = float(np.dot(s, r))
        cnt += min(v, T - v) <= obs + 1e-9
    return cnt / 2 ** len(d)


def holm(ps):
    m = len(ps)
    order = np.argsort(ps)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * ps[i]))
        adj[i] = run
    return adj


# urban outage family of UT-CHO-40 (as declared in analyze.fam_urban)
for ch, C in [("urban", "Urb"), ("denseurban", "Dense")]:
    base = ["a3tf", "d2f", "t1f", "elevf", "tted2f", "a3tfh", "d2fh", "t1fh", "elevfh", "t1f40", "elevf40",
            "d2f4050", "mst40"]
    names = {"a3tf": "Atunedf", "d2f": "Dtwof", "t1f": "Tonef", "elevf": "Elevf", "tted2f": "Ttedtwof",
             "a3tfh": "Atunedfh", "d2fh": "Dtwofh", "t1fh": "Tonefh", "elevfh": "Elevfh", "t1f40": "Toneforty",
             "elevf40": "Elevforty", "d2f4050": "Dtwoforty", "mst40": "Mstforty"}
    u = col(f"urban/{ch}/ut40", "out")
    ps = [wil(u, col(f"urban/{ch}/{b}", "out")) for b in base]
    adj = holm(ps)
    for b, p, pa in zip(base, ps, adj):
        k = f"{C}UtfortyVs{names[b]}Out"
        sp, sh = M.get("P" + k, ""), M.get("H" + k, "")
        val = lambda s: 0.0005 if "<0.001" in s else num_macro_str(s)
        def num_macro_str(s):
            s2 = s.replace("$", "").replace("p=", "")
            try:
                return float(s2)
            except ValueError:
                return float("nan")
        okp = (p < 0.001 and "<0.001" in sp) or abs(p - num_macro_str(sp)) <= 0.0015
        okh = (pa < 0.001 and "<0.001" in sh) or abs(pa - num_macro_str(sh)) <= 0.0015
        check(okp and okh, f"E {k}: p {sp} vs scipy {p:.4f}; Holm {sh} vs {pa:.4f}")
        # the paired difference and its CI
        bb = col(f"urban/{ch}/{b}", "out")
        kk = sorted(set(u) & set(bb))
        d = np.array([u[i] - bb[i] for i in kk])
        if len(d) > 1:
            close("D" + k, float(d.mean()))
            close("D" + k + "Ci", float(stats.t.ppf(0.975, len(d) - 1) * d.std(ddof=1) / math.sqrt(len(d))))

# co-channel interference. Means and the urban UT-CHO-40 family per variant,
# recomputed here; and the regression identity: no interference option means the main runs.
IV = {"k3c176": "Kthree", "k1c176": "Kone", "k1c10": "Konehi"}
for v, V in IV.items():
    for ch, C, pols in [("los", "Los", ["mst", "ut", "t1", "a3tf"]), ("urban", "Urb", ["ut40", "d2f", "mst40", "d2f4050"])]:
        for pol, P in [(x, {"mst": "Mst", "ut": "Ut", "t1": "Tone", "a3tf": "Atunedf", "ut40": "Utforty",
                            "d2f": "Dtwof", "mst40": "Mstforty", "d2f4050": "Dtwoforty"}[x]) for x in pols]:
            tag = f"interf/{ch}/{v}/{pol}"
            if summ(tag):
                close(f"If{C}{V}{P}Ho", mean(tag, "ho"))
                close(f"If{C}{V}{P}Out", mean(tag, "out"))
                close(f"If{C}{V}{P}Snr", mean(tag, "snr"))
    bases = [b for b in ["d2f", "t1f", "elevf", "d2f4050", "mst40"] if summ(f"interf/urban/{v}/{b}")]
    names = {"d2f": "Dtwof", "t1f": "Tonef", "elevf": "Elevf", "d2f4050": "Dtwoforty", "mst40": "Mstforty"}
    u = col(f"interf/urban/{v}/ut40", "out")
    if u and bases:
        ps = [wil(u, col(f"interf/urban/{v}/{b}", "out")) for b in bases]
        adj = holm(ps)
        for b, pr, pa in zip(bases, ps, adj):
            k = f"IfUrb{V}UtfortyVs{names[b]}Out"
            sh = M.get("H" + k, "")
            s2 = sh.replace("$", "").replace("p=", "")
            okh = (pa < 0.001 and "<0.001" in sh) or (s2 and abs(pa - float(s2)) <= 0.0015)
            check(bool(okh), f"E {k}: Holm {sh} vs {pa:.4f}")


# ---------------------------------------------- F. corrections and control runs
def late_pct(tag):
    o = {}
    for k, r in summ(tag).items():
        ev = [e for e in csv.DictReader(open(os.path.join(r["_dir"], "events.csv"))) if e["kind"] == "ho"]
        o[k] = 100.0 * sum(1 for e in ev if e["true_serv_exit_s"] and float(e["true_serv_exit_s"]) <= 0) / max(1, len(ev))
    return o


def pdiff(A, B):
    k = sorted(set(A) & set(B))
    return np.array([A[i] - B[i] for i in k])


def ci95(d):
    return float(stats.t.ppf(0.975, len(d) - 1) * d.std(ddof=1) / math.sqrt(len(d)))


PN = {"mst": "Mst", "tte": "Tte", "ut": "Ut", "a3tf": "Atunedf", "d2f": "Dtwof", "t1f": "Tonef", "elevf": "Elevf",
      "mst40": "Mstforty", "tte40": "Tteforty", "ut40": "Utforty"}
# issue 4: the real-shell control (MST + margin at the coverage edge), seeds 1-10
for ch, C in [("los", "Los"), ("suburban", "Sub")]:
    for pol in ["mst", "tte", "ut"]:
        L_ = late_pct(f"realtle/{ch}/{pol}")
        check(len(L_) == 10, f"F realtle/{ch}/{pol}: ten seeds ({len(L_)})")
        close(f"Tle{C}{PN[pol]}Late", float(np.mean(list(L_.values()))))
        close(f"Tle{C}{PN[pol]}Out", mean(f"realtle/{ch}/{pol}", "out"))
        close(f"Tle{C}{PN[pol]}Ho", mean(f"realtle/{ch}/{pol}", "ho"))
    for a, b in [("tte", "mst"), ("ut", "tte")]:
        d = pdiff(late_pct(f"realtle/{ch}/{a}"), late_pct(f"realtle/{ch}/{b}"))
        close(f"Tle{C}{PN[a]}Vs{PN[b]}Late", float(d.mean()))
        close(f"Tle{C}{PN[a]}Vs{PN[b]}LateCi", ci95(d))
        do = pdiff(col(f"realtle/{ch}/{a}", "out"), col(f"realtle/{ch}/{b}", "out"))
        close(f"DTle{C}{PN[a]}Vs{PN[b]}Out", float(do.mean()))
        close(f"DTle{C}{PN[a]}Vs{PN[b]}OutCi", ci95(do))
        sp = M.get(f"Tle{C}{PN[a]}Vs{PN[b]}LateP", "")
        p = wil(late_pct(f"realtle/{ch}/{a}"), late_pct(f"realtle/{ch}/{b}"))
        s2 = sp.replace("$", "").replace("p=", "")
        check((p < 0.001 and "<0.001" in sp) or (s2 not in ("", "--") and abs(p - float(s2)) <= 0.0015),
              f"E Tle{C}{PN[a]}Vs{PN[b]}LateP: {sp} vs scipy {p:.4f}")
# the eight-test family of the control (outage and handover rate), Holm
fam = []
for ch, C in [("los", "Los"), ("suburban", "Sub")]:
    for a, b in [("tte", "mst"), ("ut", "tte")]:
        for mk, MK in [("out", "Out"), ("ho", "Ho")]:
            fam.append((f"Tle{C}{PN[a]}Vs{PN[b]}{MK}", wil(col(f"realtle/{ch}/{a}", mk), col(f"realtle/{ch}/{b}", mk))))
adj = holm([p for _, p in fam])
for (k, p), pa in zip(fam, adj):
    sh = M.get("H" + k, "")
    s2 = sh.replace("$", "").replace("p=", "")
    check((pa < 0.001 and "<0.001" in sh) or (s2 not in ("", "--") and abs(pa - float(s2)) <= 0.0015),
          f"E {k}: Holm {sh} vs {pa:.4f}")
# issue 11: UT-CHO at the coverage edge under blockage with measured 0-6 h error
for ch, C in [("suburban", "Sub"), ("urban", "Urb")]:
    L_ = late_pct(f"ephblk/{ch}/ut")
    check(len(L_) == 10, f"F ephblk/{ch}/ut: ten seeds ({len(L_)})")
    close(f"Blk{C}UtLate", float(np.mean(list(L_.values()))))
    close(f"Blk{C}UtOut", mean(f"ephblk/{ch}/ut", "out"))
    d = pdiff(col(f"ephblk/{ch}/ut", "out"), col(f"ephblk/{ch}/mst", "out"))
    close(f"Blk{C}UtVsMstOut", float(d.mean()))
    exp0 = "main" if ch == "suburban" else "urban"
    d = pdiff(col(f"ephblk/{ch}/ut", "out"), col(f"{exp0}/{ch}/ut", "out"))
    close(f"Blk{C}UtDOut", float(d.mean()))
# issue 1: under blockage, late shares away from the coverage edge equal the exact-ephemeris ones
for ch, C in [("suburban", "Sub"), ("urban", "Urb")]:
    exp0 = "main" if ch == "suburban" else "urban"
    for pol in ["a3tf", "d2f", "t1f", "elevf", "mst40", "tte40", "ut40"]:
        la = float(np.mean(list(late_pct(f"ephblk/{ch}/{pol}").values())))
        lb = float(np.mean(list(late_pct(f"{exp0}/{ch}/{'mst40' if pol == 'tte40' else pol}").values())))
        check(abs(la - lb) < 0.005, f"F issue 1: ephblk/{ch}/{pol} late {la:.4f} % = exact {lb:.4f} % to 0.01")
        close(f"Blk{C}{PN[pol]}LateExact", lb)
# ... and the standardized conditions were never late in the line-of-sight error sweeps
nl = 0
for exp in ["eph", "ephemp"]:
    for pol in ["d2", "t1", "elev"]:
        for f in glob.glob(os.path.join(RUNS, exp, "*", pol, "seed*", "events.csv")):
            for e in csv.DictReader(open(f)):
                nl += e["kind"] == "ho" and bool(e["true_serv_exit_s"]) and float(e["true_serv_exit_s"]) <= 0
check(nl == 0, f"F issue 1: T1, D2 and elevation never late in the LOS error sweeps ({nl} late)")
# real shell: only A3-tuned (which reads no ephemeris) and the coverage-edge rules are ever late
for ch in ["los", "suburban"]:
    for pol in ["d2f", "t1f", "elevf", "mst40", "tte40", "ut40"]:
        v = list(late_pct(f"realtle/{ch}/{pol}").values())
        check(len(v) == 10 and max(v) == 0.0, f"F issue 1: realtle/{ch}/{pol} never late")
# issue 2: ping-pong (urban max over the table's policies; dense urban at the 40 deg contour)
urbp = ["a3tf", "d2f", "t1f", "elevf", "tted2f", "a3tfh", "d2fh", "t1fh", "elevfh", "t1f40", "elevf40", "d2f4050",
        "mst", "mst40", "ut", "ut40", "ut40_norank", "ut40_noprep", "ut40_noreserve", "ut40_qesc", "ut40_althold",
        "ut40_oldesc"]
close("UrbPponeMax", max(mean(f"urban/urban/{p}", "ppone") for p in urbp))
dp = [mean(f"urban/denseurban/{p}", "ppone") for p in ["elevf40", "d2f4050", "mst40", "ut40"]]
close("DensePponeFortyMin", min(dp))
close("DensePponeFortyMax", max(dp))
check(mean("urban/denseurban/t1f40", "ppone") < 1.0, "F issue 2: dense-urban T1-40 below 1 %")
# issue 3 and 10: interference
for v, V in [("k3c176", "Kthree"), ("k1c176", "Kone")]:
    dr = []
    for pol in ["a3tf", "d2f", "t1f", "elevf", "t1", "d2f4050", "mst", "mst40", "ut", "ut40"]:
        ref = f"main/los/{pol}"
        dr.append(-float(pdiff(col(f"interf/los/{v}/{pol}", "snr"), col(ref, "snr")).mean()))
    check(min(dr) > 0, f"F issue 3: SINR falls for every policy at {v}")
    close(f"IfLosDrop{V}Min", min(dr))
    close(f"IfLosDrop{V}Max", max(dr))
ch_ = []
for pol in ["a3tf", "d2f", "t1f", "elevf", "d2f4050", "mst", "mst40", "ut", "ut40"]:
    ch_.append(float(pdiff(col(f"interf/urban/k1c176/{pol}", "out"), col(f"urban/urban/{pol}", "out")).mean()))
check(max(ch_) < 0, f"F issue 10: urban outage falls at reuse 1 for every policy with the fallback (max {max(ch_):.4f})")
close("IfUrbKoneOutChgMax", -min(ch_))
# issue 5: MST-40 against the standardized rows, urban outage, detected and lower
base = ["a3tf", "d2f", "t1f", "elevf", "tted2f", "a3tfh", "d2fh", "t1fh", "elevfh", "t1f40", "elevf40", "d2f4050"]
m40 = col("urban/urban/mst40", "out")
ps = [wil(m40, col(f"urban/urban/{b}", "out")) for b in base]
adj = holm(ps)
low = sum(1 for b, pa in zip(base, adj) if pa < 0.05 and pdiff(m40, col(f"urban/urban/{b}", "out")).mean() < 0)
check(str(low) == M.get("FamUrbMstOutSigLower"), f"F issue 5: MST-40 detectably lower than {low} of 12 ({M.get('FamUrbMstOutSigLower')})")
# issue 6: conditional reconfigurations of the standardized conditions, LOS
rr = {p: float(np.mean([float(r["rrc_per_ue_min"]) for r in summ(f"main/los/{p}").values()]))
      for p in ["a3", "a3t", "d2", "t1", "elev", "a3tf", "d2f", "t1f", "elevf"]}
close("StdRrcLosMin", rr["elev"])
close("StdRrcLosMax", rr["a3"])
check(min(rr, key=rr.get) in ("elev", "elevf") and max(rr, key=rr.get) == "a3", "F issue 6: extremes are elevation and A3")
# issue 7: dense-urban HOF of UT-CHO-40 against D2 + net
close("UrbDenseUtfortyHof", mean("urban/denseurban/ut40", "hof"))
close("UrbDenseAtunedfHof", mean("urban/denseurban/a3tf", "hof"))
# issue 13: projection of the inertial offset on the relative track
lat, inc = math.radians(30), math.radians(53)
nn = math.sqrt(mu / (RE_ + h_) ** 3)
fac = (1 - om / nn * math.cos(inc)) / (rel / vg) ** 2
close("RotKappaProj", fac)
ui = np.array([math.sqrt(1 - saz ** 2), saz])
ur = np.array([vg * math.sqrt(1 - saz ** 2), vg * saz - om * RE_ * math.cos(lat)])
ur = ur / np.linalg.norm(ur)
close("RotNormalPct", 100 * abs(ui[0] * ur[1] - ui[1] * ur[0]))
# issue 16: UT-CHO-40 and MST-40 under line of sight
d = np.abs(pdiff(col("main/los/ut40", "ho"), col("main/los/mst40", "ho")))
close("SeedMaxLosUtfortyVsMstfortyHo", float(d.max()))
check(str(int((d < 1e-12).sum())) == M.get("SeedSameLosUtfortyVsMstfortyHo"), "F issue 16: identical seeds")

# ---------------------------------------------- G. pre-specified age-aware margin
# recompute_agemargin.py (no analyze.py import) recomputes every metric and test of the
# pre-registration; here its numbers are compared with the printed macros, and the design
# of the experiment is checked on the raw outputs.
import subprocess as _sp  # noqa: E402
_sp.run([sys.executable, os.path.join(HERE, "recompute_agemargin.py")], check=True, capture_output=True)
AMJ = json.load(open(os.path.join(HERE, "recompute_agemargin.json")))
AMSET = {"AmLos": "realtle/los", "AmSub": "realtle/suburban", "AmGau": "ephage/los"}
AMPOL = {"mst": "Mst", "tte": "Fix", "ttefz": "Fz", "tteage": "Age", "ut": "Ut", "utage": "Utage"}
for S, pre in AMSET.items():
    for pol, P in AMPOL.items():
        mm = AMJ["means"].get(f"{pre}/{pol}")
        if not mm:
            continue
        check(mm["n"] == 10, f"G {pre}/{pol}: ten seeds ({mm['n']})")
        close(f"{S}{P}Late", mm["late"])
        close(f"{S}{P}Aband", mm["aband"])
        close(f"{S}{P}Ho", mm["ho"])
        close(f"{S}{P}Out", mm["out"])
def pclose(k, val):
    s = M.get(k, "")
    if "<" in s:
        check(val < num_macro(k), f"D {k} = {s}, independent {val:.5f}")
    else:
        v = s.replace("$", "").replace("p=", "")
        dec = len(v.split(".")[1]) if "." in v else 0
        check(abs(num_macro(k) - val) <= 0.5 * 10 ** (-dec) + 1e-9, f"D {k} = {s}, independent {val:.4f}")


FAMN = {"P": "AmFamP", "S1": "AmFamSone", "S2": "AmFamStwo", "S3": "AmFamSthree", "G": "AmFamG", "G2": "AmFamGtwo"}
NM = {"P": {"los_late_ni": "AmLosAgeVsFixLate", "los_aband": "AmLosAgeVsFixAband",
            "sub_late_ni": "AmSubAgeVsFixLate", "sub_aband": "AmSubAgeVsFixAband"},
      "G": {"mst_late_ni": "AmGauAgeVsFixLate", "mst_aband": "AmGauAgeVsFixAband",
            "ut_late_ni": "AmGauUtageVsUtLate", "ut_aband": "AmGauUtageVsUtAband"},
      "G2": {"late": "AmGauAgeVsFzLate", "aband": "AmGauAgeVsFzAband"}}
for ch, C in [("los", "AmLos"), ("suburban", "AmSub")]:
    NM.setdefault("S1", {}).update({f"{ch}_ut_late_ni": f"{C}UtageVsUtLate", f"{ch}_ut_aband": f"{C}UtageVsUtAband",
                                    f"{ch}_ut_ho": f"{C}UtageVsUtHo", f"{ch}_ut_out": f"{C}UtageVsUtOut"})
    NM.setdefault("S2", {}).update({f"{ch}_ho": f"{C}AgeVsFixHo", f"{ch}_out": f"{C}AgeVsFixOut",
                                    f"{ch}_hof": f"{C}AgeVsFixHof"})
    NM.setdefault("S3", {}).update({f"{ch}_age_vs_fz_late": f"{C}AgeVsFzLate", f"{ch}_age_vs_fz_aband": f"{C}AgeVsFzAband",
                                    f"{ch}_fz_vs_fixed_late": f"{C}FzVsFixLate", f"{ch}_fz_vs_fixed_aband": f"{C}FzVsFixAband"})
for fam, names in NM.items():
    R_ = AMJ[fam]
    for k, mname in names.items():
        if k not in R_:
            continue
        r = R_[k]
        close(mname, r["diff"])
        close(mname + "Ci", r["ci"])
        pclose(mname + "P", r["p_test"])
        pclose(mname + "H", r["p_holm"])
        check(str(r["lower"]) == M.get(mname + "Lower"), f"G {mname}: seeds lower {r['lower']} ({M.get(mname + 'Lower')})")
# the decision rule, recomputed
P_ = AMJ["P"]
ok = all(r["p_holm"] < 0.05 for r in P_.values()) and all(r["diff"] < 0 for k, r in P_.items() if r["kind"] == "sup")
check(("supported" if ok else "not supported") == M.get("AmVerdict"), f"G pre-registered verdict: {'supported' if ok else 'not supported'} ({M.get('AmVerdict')})")
# design: the age-aware runs budget by age (margin rises with the bin sigma in every run),
# the fixed runs do not (margin spread below 0.2 s at 5-95 %), and every age run logs its ages
for pre in ["realtle/los", "realtle/suburban", "ephage/los"]:
    rising = spread = logged = 0
    for r in summ(f"{pre}/tteage").values():
        a = {int(x["sat"]): float(x["sig_a_m"]) for x in csv.DictReader(open(os.path.join(r["_dir"], "sat_ages.csv")))}
        logged += 1
        ev = [e for e in csv.DictReader(open(os.path.join(r["_dir"], "events.csv"))) if e["kind"] == "ho" and e["trigger"] == "planned"]
        lo = [float(e["margin_s"]) for e in ev if a[int(e["from_sat"])] < 2700]
        hi = [float(e["margin_s"]) for e in ev if a[int(e["from_sat"])] > 8000]
        rising += bool(lo and hi and np.median(hi) > np.median(lo))
    check(logged == 10, f"G {pre}/tteage: sat_ages.csv in every run ({logged})")
    check(rising >= 9, f"G {pre}/tteage: margin larger for old than young sets in >= 9 of 10 runs ({rising})")
    ms = [float(e["margin_s"]) for r in summ(f"{pre}/ttefz").values()
          for e in csv.DictReader(open(os.path.join(r["_dir"], "events.csv"))) if e["kind"] == "ho"]
    q = np.percentile(ms, [5, 95]) if ms else [0, 99]
    check(q[1] - q[0] < 0.2, f"G {pre}/ttefz: fixed margin, 5-95 % spread {q[1] - q[0]:.3f} s")
# late share by the element-set age of the serving satellite (< 6 h and >= 12 h), pooled seeds
for ch, C in [("los", "AmLos"), ("suburban", "AmSub")]:
    for pol, P in [("tte", "Fix"), ("tteage", "Age")]:
        cnt = {"y": [0, 0], "o": [0, 0]}
        tot = 0
        for k, r in summ(f"realtle/{ch}/{pol}").items():
            af = os.path.join(RUNS, f"realtle/{ch}/tteage", os.path.basename(r["_dir"]), "sat_ages.csv")
            ag = {int(x["sat"]): float(x["age0_h"]) for x in csv.DictReader(open(af))}
            for e in csv.DictReader(open(os.path.join(r["_dir"], "events.csv"))):
                if e["kind"] != "ho" or not e["true_serv_exit_s"]:
                    continue
                tot += 1
                a_ = ag[int(e["from_sat"])]
                b_ = "y" if a_ < 6 else ("o" if a_ >= 12 else None)
                if b_:
                    cnt[b_][0] += 1
                    cnt[b_][1] += float(e["true_serv_exit_s"]) <= 0
        close(f"{C}{P}LateYoung", 100 * cnt["y"][1] / cnt["y"][0])
        close(f"{C}{P}LateOld", 100 * cnt["o"][1] / cnt["o"][0])
        close(f"{C}{P}ShareYoung", 100 * cnt["y"][0] / tot)
# regression of the round-4 binary (regress_r4.py)
rl = open(os.path.join(HERE, "regress_r4.log")).read()
check("REGRESSION PASS" in rl and rl.count("identical") >= 12, "G regression of the cc0e3e190 binary: byte-identical")
# the 3.45 km fixed budget of the real-shell runs, and the measured all-age sigma
cj = json.load(open(os.path.join(RUNS, "realtle/los/tte/seed01/cmd.json")))["cmd"]
budget = [float(a.split("=")[1]) for a in cj if a.startswith("--assumedSigmaAlongM=")][0]
close("TleBudgetA", budget / 1000)

# ------------------------------------------------ H. recent-work baselines and Proposition 1 (Oct 2026)
REC = [("tted2f", "Ttedtwof"), ("geof", "Geof"), ("wangf", "Wangf"), ("zabf", "Zabf")]
ENV = [("main", "los", "Los"), ("main", "suburban", "Sub"), ("urban", "urban", "Urb"), ("urban", "denseurban", "Dense")]
fam, names_h = [], []
for key, KN in REC:
    for exp, ch, C in ENV:
        A = col(f"{exp}/{ch}/{key}", "ho")
        check(len(A) == 10, f"H {exp}/{ch}/{key}: ten seeds ({len(A)})")
        close(f"Rec{C}{KN}Ho", float(np.mean(list(A.values()))))
        close(f"Rec{C}{KN}Out", mean(f"{exp}/{ch}/{key}", "out"))
        ref, mk = ("ut", "ho") if exp == "main" else ("ut40", "out")
        a, b = col(f"{exp}/{ch}/{ref}", mk), col(f"{exp}/{ch}/{key}", mk)
        ks = sorted(set(a) & set(b))
        d = [a[i] - b[i] for i in ks]
        close(f"Rec{C}{KN}Diff", float(np.mean(d)))
        if mk == "ho":
            close(f"Rec{C}{KN}Rel", 100.0 * float(np.mean(d)) / float(np.mean([b[i] for i in ks])))
        fam.append(wil(a, b))
        names_h.append(f"Rec{C}{KN}P")
for k, pa in zip(names_h, holm(fam)):
    sp = M.get(k, "").replace("$", "").replace("p=", "")
    ok = ("<0.001" in sp and pa < 0.001) or (sp not in ("", "--") and abs(float(sp) - pa) <= 0.0015)
    check(ok, f"H {k}: Holm {sp} vs independent {pa:.4f}")
for key, KN in REC:
    v = late_pct(f"realtle/los/{key}")
    close(f"TleLos{KN}Late", float(np.mean(list(v.values()))))
# first-order exit error against the measured one (MST, real shell, planned exits)
eta = float(M["RotKappa"])
xs, ys = [], []
for ch in ["los", "suburban"]:
    for d in sorted(glob.glob(os.path.join(RUNS, f"realtle/{ch}/mst/seed*"))):
        cmd = json.load(open(os.path.join(d, "cmd.json")))["cmd"]
        alt = float([a.split("=")[1] for a in cmd if a.startswith("--altitudeKm=")][0]) * 1e3
        v = math.sqrt(3.986004418e14 / (6378137.0 + alt))
        off = {int(o["sat"]): o for o in csv.DictReader(open(os.path.join(d, "eph_offsets.csv")))}
        T = float(next(csv.DictReader(open(os.path.join(d, "summary.csv"))))["sim_time_s"])
        for e in csv.DictReader(open(os.path.join(d, "events.csv"))):
            if e["kind"] != "ho" or e["trigger"] != "planned" or e["true_serv_exit_s"] == "":
                continue
            o = off[int(e["from_sat"])]
            t = float(e["t_s"])
            da = float(o["d0_a_m"]) + (float(o["dT_a_m"]) - float(o["d0_a_m"])) * t / T
            xs.append(-eta * da / v)
            ys.append(float(e["pred_serv_tte_s"]) - float(e["true_serv_exit_s"]))
sl = stats.linregress(xs, ys)
close("PropSlope", sl.slope)
close("PropCorr", sl.rvalue)
check(str(len(xs)) == M.get("PropN"), f"H first-order sample size {len(xs)} = PropN {M.get('PropN')}")

# z sweep of Proposition 1 (real shell, line of sight and suburban): late share, outage and
# handover rate per z, z = 0 identical to MST seed by seed, the paired test of the best z against
# z = 2, and the closed form (eq:zstar) recomputed from the printed cost parameters
ZT = {0: "z0", 0.5: "z0.5", 1: "z1", 1.5: "z1.5", 2: None, 3: "z3", 4: "z4", 5: "z5", 6: "z6"}
ZN = {0: "Zero", 0.5: "Half", 1: "One", 1.5: "OneHalf", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six"}
for env, EN in (("los", "Los"), ("suburban", "Sub")):
    ztag = lambda z, env=env: f"zsweep/{env}/{ZT[z]}" if ZT[z] else f"realtle/{env}/tte"
    for z in ZT:
        L_ = late_pct(ztag(z))
        check(len(L_) == 10, f"H {env} {ztag(z)}: ten seeds ({len(L_)})")
        close(f"Zs{EN}{ZN[z]}Late", float(np.mean(list(L_.values()))))
        close(f"Zs{EN}{ZN[z]}Out", mean(ztag(z), "out"))
        close(f"Zs{EN}{ZN[z]}Ho", mean(ztag(z), "ho"))
    m0, mm = col(ztag(0), "out"), col(f"realtle/{env}/mst", "out")
    h0, hm = col(ztag(0), "ho"), col(f"realtle/{env}/mst", "ho")
    check(m0 == mm and h0 == hm, f"H {env} z = 0 reproduces MST seed by seed (outage and handover rate)")
    zb = min(ZT, key=lambda z: mean(ztag(z), "out"))
    check(abs(zb - float(M[f"Prop{EN}Zsim"])) < 1e-9, f"H {env} lowest simulated outage at z = {zb} (Prop{EN}Zsim {M['Prop' + EN + 'Zsim']})")
    d = pdiff(col(ztag(zb), "out"), col(ztag(2), "out"))
    close(f"Zs{EN}OptVsTwoOut", float(d.mean()))
    close(f"Zs{EN}OptVsTwoOutCi", ci95(d))
    pz = wil(col(ztag(zb), "out"), col(ztag(2), "out"))
    sp = M.get(f"Zs{EN}OptVsTwoOutP", "").replace("$", "").replace("p=", "")
    check(sp not in ("", "--") and abs(float(sp) - pz) <= 0.0015, f"H {env} z sweep best vs z=2: p {sp} vs {pz:.4f}")
    ce, kk, sg, chh = float(M[f"Prop{EN}Ceff"]), float(M[f"Prop{EN}K"]) * 1e-3, float(M[f"Prop{EN}Sigma"]), float(M[f"Prop{EN}Ch"]) * 1e-3
    check(abs(ce - (float(M[f"Prop{EN}Cl"]) + float(M[f"Prop{EN}Beta"]) * chh)) <= 0.011, f"H {env} c_eff = c_l + beta c_h")
    check(ce > math.sqrt(2 * math.pi) * kk * sg * (chh + ce / 2), f"H {env} condition of Proposition 1 holds")
    zs = math.sqrt(2 * math.log(ce / (math.sqrt(2 * math.pi) * kk * sg * chh)))
    check(abs(zs - float(M[f"Prop{EN}Zstar"])) <= 0.03, f"H {env} z* from printed parameters {zs:.3f} vs Prop{EN}Zstar {M['Prop' + EN + 'Zstar']}")


# ---------------------------------------------- I. handover-rate analysis (Section IV-A, Fig. 2)
# Independent recomputation of the analytic rates (own geometry code, not rate_model.py) and of
# the Monte Carlo-based macros from the raw Monte Carlo logs.
import ast as _ast
RM = {}
for _l in open(os.path.join(HERE, "analysis_rate", "macros_rate.tex")):
    _g = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s*$", _l)
    if _g:
        RM[_g.group(1)] = _g.group(2)
_R, _MU, _WE = 6378.137, 398600.4418, 7.2921159e-5
def _geo(N, h, inc, eQ, lat=30.0):
    r = _R + h; v = math.sqrt(_MU / r); vg = v * _R / r
    psi = math.asin(_R / r * math.cos(math.radians(eQ)))
    rho = _R * (math.pi / 2 - math.radians(eQ) - psi)
    i, p = math.radians(inc), math.radians(lat)
    cb = math.cos(i) / math.cos(p)
    vrel = math.hypot(vg * cb - _WE * _R * math.cos(p), vg * math.sqrt(1 - cb * cb))
    sig = N / (2 * math.pi ** 2 * _R ** 2 * math.sqrt(math.sin(i) ** 2 - math.sin(p) ** 2))
    return rho, vrel, sig, vg / vrel
_rho, _vrel, _sig, _eta = _geo(1584, 550, 53, 25)
check(abs(_eta - float(M["RotKappa"])) < 0.002, f"I eta from eq. vrel {_eta:.4f} = RotKappa {M['RotKappa']}")
_hn = 60 * 4 * math.sqrt(_sig) * _vrel / math.pi
check(abs(_hn - float(RM["RateNearAnRef"])) < 0.0015, f"I nu_near {_hn:.4f} = RateNearAnRef {RM['RateNearAnRef']}")
check(abs(60 * _vrel / (2 * _rho) - float(RM["RateFloorRef"])) < 0.0015, "I floor v_rel/(2 rho_Q)")
check(abs(math.pi * _rho ** 2 * _sig - float(RM["RateAnNbar"])) < 0.15, f"I Nbar planar {math.pi*_rho**2*_sig:.2f} ~ RateAnNbar {RM['RateAnNbar']}")
_geo_ho, _mst_ho = float(M["MainLosGeofHo"]), float(M["MainLosMstHo"])
check(f"{_geo_ho / _mst_ho:.1f}" == RM["RateNsRatio"], "I ns-3 ratio geodesic/MST")
check(f"{100 * abs(_hn - _geo_ho) / _geo_ho:.0f}" == RM["RateNsNearErr"], "I ns-3 nearest-rule error")
_mc = []
for _fn in ("mc_rate_part1.log", "mc_rate_part2.log"):
    _p = os.path.join(HERE, "analysis_rate", _fn)
    if os.path.exists(_p):
        for _l in open(_p):
            if _l.startswith("{"):
                _mc.append(_ast.literal_eval(re.sub(r"np\.float64\(([^)]*)\)", r"\1", _l.strip())))
for _d in _mc:
    _r2, _v2, _s2, _ = _geo(_d["N"], _d["h"], _d["inc"], _d["eQ"])
    check(abs(60 * 4 * math.sqrt(_s2) * _v2 / math.pi - _d["Hnear"]) < 1e-3, f"I {_d['N']}/{_d['eQ']}: logged nu_near = eq. near")
_d25 = {(_d["N"], _d["eQ"]): _d for _d in _mc if _d["h"] == 550}
if (1584, 25) in _d25 and "RateExitMcRef" in RM:
    check(f"{_d25[(1584, 25)]['mc_Hexit']:.3f}" == RM["RateExitMcRef"], "I Monte Carlo MST rate macro")
_dn = [100 * abs(_d["Hnear"] - _d["mc_Hnear"]) / _d["mc_Hnear"] for _d in _mc if _d["Nvis"] >= 3 and _d["eQ"] == 25]
if _dn:
    check(f"{max(_dn):.0f}" == RM["RateNearErrDense"], f"I nearest-rule error, dense shells, max {max(_dn):.1f} %")
_er, _ra, _f = 6378.137, math.radians(25.0), 2.185
_ps = math.asin(_er / (_er + 550) * math.cos(_ra)); _dd = math.sqrt((_er + 550) ** 2 - (_er * math.cos(_ra)) ** 2) - _er * math.sin(_ra)
_snr = 34 + 30 + 10 * math.log10(20) + 12 * math.log10(math.cos(_ps)) - (32.45 + 20 * math.log10(_f) + 20 * math.log10(_dd * 1e3)) \
    - 6.22 / _f ** 1.5 - 3 - (-174 + 10 * math.log10(20e6) + 7)
check(f"{_snr + 8:.1f}" == RM["LinkMarginEmin"], f"I link margin at e_min {_snr + 8:.2f} dB")

# ns-3 geodesic CHO at 792 and 3168 satellites (density/*/los/geof), recomputed from the summaries
for _sh, _mk in (("w792", "RateNsGeoSeven"), ("w3168", "RateNsGeoThirty")):
    _v = [float(next(csv.DictReader(open(_f)))["ho_per_ue_min"])
          for _f in glob.glob(os.path.join(HERE, "runs_r2", "density", _sh, "los", "geof", "seed*", "summary.csv"))]
    check(len(_v) == 10 and f"{np.mean(_v):.3f}" == RM[_mk], f"I ns-3 geodesic {_sh}: {np.mean(_v):.3f} = {_mk}")

report = "\n".join(out) + f"\n\n{len(out)} checks, {fails} failed\n"
open(os.path.join(HERE, "verify_report.txt"), "w").write(report)
print(report)
sys.exit(1 if fails else 0)
