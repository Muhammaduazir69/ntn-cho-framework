#!/usr/bin/env python3
"""November 2026 campaign runs, specified in
analysis_protocol.md (read that first).

    python3 run_nov.py tune            # group T, seeds 1-10, 600 s (exploratory)
    python3 run_nov.py select          # apply the pre-specified rival rule -> runs_nov/rival_choice.json
    python3 run_nov.py conf            # group C, seeds 11-20, 1200 s (confirmatory)
    python3 run_nov.py oos             # group O, real shell +12 h, seeds 11-20
    python3 run_nov.py soft            # group S, soft coverage edge, seeds 1-10

One binary (bin/p2-campaign-nov) with the frozen libraries in bin/lib; its sha256 is
written into every cmd.json. At most 3 jobs at nice 19."""
import csv, hashlib, itertools, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_campaign as rc  # noqa: E402  (policy definitions only)

NS3 = rc.NS3
BIN = os.path.join(HERE, "bin", "p2-campaign-nov")
LIB = os.path.join(HERE, "bin", "lib")
RUNS = os.path.join(HERE, "runs_nov")
R2 = os.path.join(HERE, "runs_r2")
P = rc.POLICIES
JOBS = 3
SHA = hashlib.sha256(open(BIN, "rb").read()).hexdigest()
ENVS = ["los", "suburban", "urban", "denseurban"]


def tted2(a, b):
    return ["--policy=d2", f"--d2ServElevDeg={a}", f"--d2CandElevDeg={b}", "--tteRank=1"] + rc.NET


def geo(hm):
    return ["--policy=geo", f"--geoHysM={hm}"] + rc.NET


TUNE = {"tted2f_25_35": tted2(25, 35), "tted2f_30_40": tted2(30, 40),
        "geof_h50": geo(50000), "geof_h100": geo(100000), "geof_h200": geo(200000)}
# reference settings already run on seeds 1-10 (runs_r2)
REF = {"tted2f": "tted2f_35_45", "geof": "geof_h10"}
FAM = {"tted2f": ["tted2f_35_45", "tted2f_25_35", "tted2f_30_40"],
       "geof": ["geof_h10", "geof_h50", "geof_h100", "geof_h200"]}
ALLSET = dict(TUNE, tted2f_35_45=P["tted2f"], geof_h10=P["geof"])
# Amendment 1: the geodesic choice sat at the grid edge (200 km); extend the family
TUNE2 = {"geof_h300": geo(300000), "geof_h500": geo(500000), "geof_h1000": geo(1000000)}
ALLSET.update(TUNE2)
FAM2 = {"geof": FAM["geof"] + list(TUNE2)}


def jobs_tune():
    for env, (name, args), s in itertools.product(ENVS, TUNE.items(), range(1, 11)):
        yield f"tune/{env}/{name}", s, args + [f"--channel={env}"] + rc.RUN


def _summ(d):
    f = os.path.join(d, "summary.csv")
    if not os.path.exists(f):
        return None
    r = list(csv.DictReader(open(f)))
    return r[0] if r and "ho_per_ue_min" in r[0] else None


def _fam_dir(env, setting):
    base = "main" if env in ("los", "suburban") else "urban"
    if setting == "tted2f_35_45":
        return os.path.join(R2, base, env, "tted2f")
    if setting == "geof_h10":
        return os.path.join(R2, base, env, "geof")
    return os.path.join(RUNS, "tune", env, setting)


def jobs_tune2():
    for env, (name, args), s in itertools.product(ENVS, TUNE2.items(), range(1, 11)):
        yield f"tune/{env}/{name}", s, args + [f"--channel={env}"] + rc.RUN


def select(fams=None, outname="rival_choice.json"):
    out = {"rule": "lowest mean HO rate among settings with mean outage <= 1.05 x the family's lowest "
                   "mean outage in that environment; ties to the reference setting", "choice": {}, "table": {}}
    for fam, sets in (fams or FAM).items():
        for env in ENVS:
            rows = {}
            for sname in sets:
                d = _fam_dir(env, sname)
                S = [_summ(os.path.join(d, f"seed{s:02d}")) for s in range(1, 11)]
                if any(x is None for x in S):
                    sys.exit(f"missing tuning rows: {d}")
                rows[sname] = (st.mean(float(x["ho_per_ue_min"]) for x in S),
                               st.mean(100 * float(x["outage_frac"]) for x in S))
            omin = min(v[1] for v in rows.values())
            ok = [k for k, v in rows.items() if v[1] <= 1.05 * omin]
            best = min(ok, key=lambda k: (rows[k][0], 0 if k == sets[0] else 1))
            out["choice"][f"{fam}/{env}"] = best
            out["table"][f"{fam}/{env}"] = rows
    os.makedirs(RUNS, exist_ok=True)
    json.dump(out, open(os.path.join(RUNS, outname), "w"), indent=1)
    print(json.dumps(out["choice"], indent=1))


CONF_LS = ["a3tf", "d2f", "t1f", "zabf", "mst", "ut"]
CONF_URB = ["d2f", "mst40", "ut40"]
RUN_C = ["--numUes=100", "--simTime=1200"]


def jobs_conf(with_rivals):
    seeds = range(11, 21)
    for env, s in itertools.product(["los", "suburban"], seeds):
        for name in CONF_LS:
            yield f"conf/{env}/{name}", s, P[name] + [f"--channel={env}"] + RUN_C
    for s in seeds:
        yield "conf/los/mst_k16", s, P["mst"] + ["--maxCandidates=16", "--channel=los"] + RUN_C
    for env, name, s in itertools.product(["urban", "denseurban"], CONF_URB, seeds):
        yield f"conf/{env}/{name}", s, P[name] + [f"--channel={env}"] + RUN_C
    if with_rivals:
        ch = json.load(open(os.path.join(RUNS, "rival_choice.json")))["choice"]
        for env, fam, s in itertools.product(ENVS, ["tted2f", "geof"], seeds):
            yield f"conf/{env}/{fam}", s, ALLSET[ch[f"{fam}/{env}"]] + [f"--channel={env}"] + RUN_C


SA = 3450.0
TLE = [f"--tleTruthFile={rc.TLE_SUP}", f"--tleNetFile={rc.TLE_GP}", "--altitudeKm=470"]


def jobs_oos():
    for ch, s in itertools.product(["los", "suburban"], range(11, 21)):
        yield f"oos/{ch}/mst", s, P["mst"] + [f"--channel={ch}"] + rc.RUN + TLE + ["--tleStartOffsetH=12"]
        yield f"oos/{ch}/tte", s, P["tte"] + [f"--channel={ch}"] + rc.RUN + TLE + rc.assumed(SA) + [
            "--tleStartOffsetH=12"]


SOFT = ["--softEdgeDeg=10", "--softRollDbPerDeg=2"]


def jobs_soft():
    for s in range(1, 11):
        yield "soft/los/mst_hard", s, P["mst"] + ["--channel=los"] + rc.RUN + TLE
        yield "soft/los/mst", s, P["mst"] + ["--channel=los"] + rc.RUN + TLE + SOFT
        yield "soft/los/tte", s, P["tte"] + ["--channel=los"] + rc.RUN + TLE + rc.assumed(SA) + SOFT


def jobs_conf_geo2():
    ch = json.load(open(os.path.join(RUNS, "rival_choice2.json")))["choice"]
    for env, s in itertools.product(ENVS, range(11, 21)):
        yield f"conf/{env}/geof2", s, ALLSET[ch[f"geof/{env}"]] + [f"--channel={env}"] + RUN_C


def jobs_near10():
    # validation row for Proposition 2 (not a rival test): geodesic at 10 km is the
    # nearest-satellite rule; LOS, seeds 11-20, 1200 s
    for s in range(11, 21):
        yield "conf/los/geof10", s, P["geof"] + ["--channel=los"] + RUN_C


def jobs_soft2():
    # Amendment 2 (post hoc, sensitivity): roll-off sweep and a hard-edge MST + margin rerun on the
    # same libraries, so every soft/hard pair shares one binary and one library set
    for s in range(1, 11):
        yield "soft/los/tte_hard", s, P["tte"] + ["--channel=los"] + rc.RUN + TLE + rc.assumed(SA)
        for r in (5, 10):
            sf = ["--softEdgeDeg=10", f"--softRollDbPerDeg={r}"]
            yield f"soft/los/mst_r{r}", s, P["mst"] + ["--channel=los"] + rc.RUN + TLE + sf
            yield f"soft/los/tte_r{r}", s, P["tte"] + ["--channel=los"] + rc.RUN + TLE + rc.assumed(SA) + sf


# Amendment 3: 2 x 2 of execution instant x target choice on the
# fresh seeds, and fixed vs per-pass margin under inflated cross-track and radial error
def jobs_attr():
    d2 = ["--policy=d2", "--d2ServElevDeg=30", "--d2CandElevDeg=40"] + rc.NET
    t1 = ["--policy=t1", "--t1DurationS=1"] + rc.NET
    for env, s in itertools.product(["los", "suburban"], range(11, 21)):
        yield f"conf/{env}/d2f_30_40", s, d2 + [f"--channel={env}"] + RUN_C
        yield f"conf/{env}/t1f_w1", s, t1 + [f"--channel={env}"] + RUN_C


PP = {"r1": (2000.0, 2000.0, 2000.0), "r3": (2000.0, 6000.0, 6000.0), "r5": (2000.0, 10000.0, 10000.0)}


def _eph(a, c, r):
    return [f"--ephSigmaAlongM={a}", f"--ephSigmaCrossM={c}", f"--ephSigmaRadialM={r}"]


def jobs_perpass():
    for (k, (a, c, r)), s in itertools.product(PP.items(), range(11, 21)):
        yield f"perpass/{k}/tte", s, P["tte"] + ["--channel=los"] + rc.RUN + _eph(a, c, r)


def jobs_fixed():
    # the fixed arm uses, per error setting, the mean margin of the per-pass arm on the same seeds
    for (k, (a, c, r)), s in itertools.product(PP.items(), range(11, 21)):
        ms = []
        for s2 in range(11, 21):
            f = os.path.join(RUNS, "perpass", k, "tte", f"seed{s2:02d}", "events.csv")
            ms += [float(e["margin_s"]) for e in csv.DictReader(open(f)) if e["kind"] == "ho" and e["trigger"] == "planned"]
        m = st.mean(ms)
        yield f"perpass/{k}/fixed", s, P["mst"] + [f"--execMarginS={m:.4f}", "--channel=los"] + rc.RUN + _eph(a, c, r)


def jobs_fig2():
    # Amendment 3(c): confirmatory-horizon ns-3 points for Fig. 2 (1200 s, seeds 11-20, LOS)
    shells = {"w792": ["--totalSats=792", "--numPlanes=36", "--phasing=17"],
              "w3168": ["--totalSats=3168", "--numPlanes=72", "--phasing=17"]}
    for (sh, sa), pol, s in itertools.product(shells.items(), ["mst", "mst40", "geof"], range(11, 21)):
        yield f"fig2/{sh}/{pol}", s, P[pol] + sa + ["--channel=los"] + RUN_C
    for s in range(11, 21):
        yield "conf/los/mst40", s, P["mst40"] + ["--channel=los"] + RUN_C


def jobs_k3168():
    # Amendment 3(d): MST with 32 prepared candidates on the 3168-satellite shell, to test whether the
    # K = 8 prepared set explains MST's rate above the Monte Carlo there
    sa = ["--totalSats=3168", "--numPlanes=72", "--phasing=17"]
    for s in range(11, 21):
        yield "fig2/w3168/mst_k32", s, P["mst"] + sa + ["--maxCandidates=32", "--channel=los"] + RUN_C


def jobs_fixedmed():
    # Amendment 4 (post hoc): fixed margin equal to the per-pass arm's MEDIAN planned margin
    # (the mean was pulled up by a few grazing passes); same seeds
    for (k, (a, c, r)), s in itertools.product(PP.items(), range(11, 21)):
        ms = []
        for s2 in range(11, 21):
            f = os.path.join(RUNS, "perpass", k, "tte", f"seed{s2:02d}", "events.csv")
            ms += [float(e["margin_s"]) for e in csv.DictReader(open(f)) if e["kind"] == "ho" and e["trigger"] == "planned"]
        m = st.median(ms)
        yield f"perpass/{k}/fixedmed", s, P["mst"] + [f"--execMarginS={m:.4f}", "--channel=los"] + rc.RUN + _eph(a, c, r)


def jobs_off_tte():
    # Amendment 5 (post hoc): the same per-pass test with the copy error held as a fixed RTN offset
    # (--ephModel=offset), so that radial error does not drift into along-track error and the
    # static sigma of the exit-error lemma is the right budget
    for (k, (a, c, r)), s in itertools.product(PP.items(), range(11, 21)):
        yield f"perpass_off/{k}/tte", s, P["tte"] + ["--channel=los", "--ephModel=offset"] + rc.RUN + _eph(a, c, r)


def jobs_off_fixedmed():
    for (k, (a, c, r)), s in itertools.product(PP.items(), range(11, 21)):
        ms = []
        for s2 in range(11, 21):
            f = os.path.join(RUNS, "perpass_off", k, "tte", f"seed{s2:02d}", "events.csv")
            ms += [float(e["margin_s"]) for e in csv.DictReader(open(f)) if e["kind"] == "ho" and e["trigger"] == "planned"]
        m = st.median(ms)
        yield f"perpass_off/{k}/fixedmed", s, P["mst"] + [f"--execMarginS={m:.4f}", "--channel=los", "--ephModel=offset"] + rc.RUN + _eph(a, c, r)


def run_one(tag, seed, args):
    out = os.path.join(RUNS, tag, f"seed{seed:02d}")
    done = os.path.join(out, "summary.csv")
    if os.path.exists(done):
        return tag, seed, "skip", 0.0
    os.makedirs(out, exist_ok=True)
    cmd = [BIN] + args + [f"--rngRun={seed}", f"--outputDir={out}"]
    env = dict(os.environ, LD_LIBRARY_PATH=LIB)
    t0 = time.time()
    with open(os.path.join(out, "stdout.txt"), "w") as fo:
        r = subprocess.call(["nice", "-n", "19"] + cmd, cwd=NS3, stdout=fo, stderr=subprocess.STDOUT, env=env)
    with open(os.path.join(out, "cmd.json"), "w") as fj:
        json.dump({"cmd": cmd, "rc": r, "bin_sha256": SHA, "lib": LIB}, fj)
    return tag, seed, ("ok" if r == 0 else f"rc={r}"), time.time() - t0


def free_gb():
    s = os.statvfs(HERE)
    return s.f_bavail * s.f_frsize / 1e9


def run(jobs, label):
    jobs = list(jobs)
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, "PROVENANCE.txt"), "a") as fp:
        fp.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} group={label} runs={len(jobs)} bin={BIN} sha256={SHA} "
                 f"lib={LIB} lib_manifest_sha256="
                 f"{hashlib.sha256(open(os.path.join(HERE, 'bin', 'lib_manifest.sha256'), 'rb').read()).hexdigest()}\n")
    print(f"{label}: {len(jobs)} runs, {JOBS} in parallel, bin {SHA[:16]}", flush=True)
    fails = 0
    with ThreadPoolExecutor(JOBS) as ex:
        futs = []
        for j in jobs:
            futs.append(ex.submit(lambda j=j: (run_one(*j) if free_gb() > 2.0 else (j[0], j[1], "nodisk", 0.0))))
        for i, f in enumerate(as_completed(futs), 1):
            tag, seed, stt, dt = f.result()
            if stt not in ("ok", "skip"):
                fails += 1
            if stt != "skip":
                print(f"[{i}/{len(jobs)}] {tag} seed{seed} {stt} {dt:.0f}s", flush=True)
    print(f"{label} done, {fails} failed", flush=True)
    return fails


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "tune":
        run(jobs_tune(), "tune")
    elif what == "select":
        select()
    elif what == "conf_base":
        run(jobs_conf(False), "conf_base")
    elif what == "conf":
        run(jobs_conf(True), "conf")
    elif what == "oos":
        run(jobs_oos(), "oos")
    elif what == "soft":
        run(jobs_soft(), "soft")
    elif what == "tune2":
        run(jobs_tune2(), "tune2 (amendment 1)")
    elif what == "select2":
        select(FAM2, "rival_choice2.json")
    elif what == "conf_geo2":
        run(jobs_conf_geo2(), "conf_geo2 (amendment 1)")
    elif what == "near10":
        run(jobs_near10(), "near10 (validation of Prop. 2, added after Amendment 1)")
    elif what == "soft2":
        run(jobs_soft2(), "soft2 (amendment 2, post hoc sensitivity)")
    elif what == "attr":
        run(jobs_attr(), "attr (amendment 3)")
    elif what == "perpass":
        run(jobs_perpass(), "perpass (amendment 3)")
    elif what == "fixed":
        run(jobs_fixed(), "fixed (amendment 3)")
    elif what == "fig2":
        run(jobs_fig2(), "fig2 (amendment 3c)")
    elif what == "k3168":
        run(jobs_k3168(), "k3168 (amendment 3d)")
    elif what == "fixedmed":
        run(jobs_fixedmed(), "fixedmed (amendment 4)")
    elif what == "offset":
        run(jobs_off_tte(), "per-pass, offset error model (amendment 5)")
        run(jobs_off_fixedmed(), "fixed median, offset error model (amendment 5)")
    elif what == "count":
        for g in (jobs_tune, jobs_oos, jobs_soft):
            print(g.__name__, len(list(g())))
        print("conf_base", len(list(jobs_conf(False))))
