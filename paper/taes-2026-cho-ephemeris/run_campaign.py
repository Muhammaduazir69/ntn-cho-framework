#!/usr/bin/env python3
"""Run the Paper 2 simulation campaign with ntn-cho-tte-campaign.

Every result in the manuscript comes from a run this script launches. A run is
skipped when its summary.csv already exists, so the campaign can be resumed.

    python3 run_campaign.py --exp all --jobs 10
    python3 run_campaign.py --exp main --seeds 10

Experiments
    main     policies at their reference settings, LOS and suburban channels
    pareto   handover rate against SINR: each policy family over its own knob
    eph      ephemeris-error sweep (true and assumed sigmas)
    acc      predicted against realized exit time
    lb       static terminals, with the offline lower bound on handovers
    sens     minimum elevation and a 1200 km polar-ish shell
    ephemp   measured TLE ephemeris errors by element-set age (N4)
    urban    urban and dense-urban environments
    ephblk   measured 0-6 h errors under blockage, every policy (R2)
    realtle  real 53-degree Starlink shell, SupGP truth, GP network copy (R2;
             MST + margin control)
    density  Walker shells of 792 and 3168 satellites
    losmodel Markov LOS process and the all-in-clutter terminal model (R2)
    reest    1 s re-establishment (R2)
    interf   co-channel inter-satellite interference, reuse 1 and 3
    agemargin age-aware margin, pre-specified:
             real shell (los, suburban) and Gaussian error with age (los)
"""
import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
NS3 = os.path.abspath(os.path.join(HERE, "../../../ns-3-dev"))
BIN = os.path.join(NS3, "build/contrib/ntn-cho/examples/ns3.43-ntn-cho-tte-campaign-optimized")
RUNS = os.path.join(HERE, os.environ.get("RUNS_DIR", "runs_r2"))

# Reference settings of each policy. a3 is the common 3GPP-style
# setting; a3t the best A3 found in the pilot (larger offset and TTT), kept so
# the comparison is against a tuned baseline and not a strawman. T1 and D2
# require the candidate to pass the quality floor Q, which is the A4 leg
# (a4-Threshold = Q, TTT 0) that TS 38.331 Rel-17 requires on the same
# candidate. Every exit-triggered policy (MST, TTE-CHO, UT-CHO) escapes to the
# strongest prepared candidate that is better than the serving cell
# (--escapeByQuality), the same target rule as the standardized
# conditions' fallback ("+ net").
NET = ["--measFallback=1"]
ESC = ["--escapeByQuality=1"]
# UT-CHO's escape (chosen in the R2 pilot, pilot_r2.py): rank the candidates within
# 3 dB of the strongest prepared candidate by the usable-time score
UT = ["--policy=tte-cho", "--rankByUsable=1", "--prepareByTte=1", "--prepStrongest=2", "--escapeNearBestDb=3"]
POLICIES = {
    "a3": ["--policy=a3", "--a3OffsetDb=3", "--a3TttMs=160"],
    "a3t": ["--policy=a3", "--a3OffsetDb=6", "--a3TttMs=640"],
    "d2": ["--policy=d2", "--d2ServElevDeg=35", "--d2CandElevDeg=45"],
    "t1": ["--policy=t1"],
    "t1w1": ["--policy=t1", "--t1DurationS=1"],
    "elev": ["--policy=elev", "--elevFloorDeg=30"],
    # the standardized conditions with the measurement fallback ("+ net")
    "a3tf": ["--policy=a3", "--a3OffsetDb=6", "--a3TttMs=640"] + NET,
    "d2f": ["--policy=d2", "--d2ServElevDeg=35", "--d2CandElevDeg=45"] + NET,
    "t1f": ["--policy=t1"] + NET,
    "elevf": ["--policy=elev", "--elevFloorDeg=30"] + NET,
    # ... and the escape hold-off of UT-CHO ("+ net + hold")
    "a3tfh": ["--policy=a3", "--a3OffsetDb=6", "--a3TttMs=640", "--escapeHoldoffS=5"] + NET,
    "d2fh": ["--policy=d2", "--d2ServElevDeg=35", "--d2CandElevDeg=45", "--escapeHoldoffS=5"] + NET,
    "t1fh": ["--policy=t1", "--escapeHoldoffS=5"] + NET,
    "elevfh": ["--policy=elev", "--elevFloorDeg=30", "--escapeHoldoffS=5"] + NET,
    # matched exit contour (40 deg). T1 is armed from the exit of the
    # 40 deg contour and prepares only satellites predicted inside it for
    # tau_min. (A first run that also reserved the two strongest satellites let
    # T1 hand over to a reserved satellite below the contour, re-arm at once
    # and cycle; it was discarded as a broken baseline, see PROVENANCE.)
    "t1f40": ["--policy=t1", "--tteExitElevDeg=40", "--prepAboveContour=1"] + NET,
    "elevf40": ["--policy=elev", "--elevFloorDeg=40"] + NET,
    "d2f4050": ["--policy=d2", "--d2ServElevDeg=40", "--d2CandElevDeg=50"] + NET,
    # prediction-based baselines
    "eager": ["--policy=tte-eager"] + ESC,
    "tted2": ["--policy=d2", "--d2ServElevDeg=35", "--d2CandElevDeg=45", "--tteRank=1"],
    "tted2f": ["--policy=d2", "--d2ServElevDeg=35", "--d2CandElevDeg=45", "--tteRank=1"] + NET,
    # classical maximum service time (Papapetrou et al. 2004): execute at the
    # predicted exit with the decision allowance T_x only, longest predicted
    # exit; with the same safety net as every other policy
    "mst": ["--policy=mst"] + ESC,
    "mst40": ["--policy=mst", "--tteExitElevDeg=40"] + ESC,
    "mst40s": ["--policy=mst", "--tteExitElevDeg=40", "--tteFloorTttMs=640"] + ESC,
    # MST plus the per-pass ephemeris margin (TTE-CHO). Identical to MST at
    # zero ephemeris error (checked per seed), so it runs only where error is
    # injected
    "tte": ["--policy=tte-cho"] + ESC,
    "tte40": ["--policy=tte-cho", "--tteExitElevDeg=40"] + ESC,
    # single mechanisms added to MST
    "mst_prep": ["--policy=mst", "--prepareByTte=1"] + ESC,
    "mst_prepr": ["--policy=mst", "--prepareByTte=1", "--prepStrongest=2"] + ESC,
    "mst_ute": ["--policy=mst", "--rankByUsable=1"] + ESC,
    # UT-CHO: per-pass margin, usable-time ranking, preparation by the score
    # with two slots for the strongest satellites, and the fixed escape. The
    # escape hold-off is decided by the R2 pilot (pilot_r2.py); see HOLD below.
    # Recent-work baselines (Oct 2026), re-implemented from the published
    # descriptions; each with the same measurement safety net (+ net).
    # Lee and Kim, IEEE TVT 2026: geodesic distance to the sub-satellite points
    "geof": ["--policy=geo"] + NET,
    # Wang et al., IEEE TVT 2026: entropy-weighted SINR and remaining service time
    "wangf": ["--policy=a3", "--a3OffsetDb=6", "--a3TttMs=640", "--rankEntropy=1"] + NET,
    # Zaboub et al., NetSoft 2026 (and Been et al., CSCN 2025): candidates
    # prioritized by predicted residence time, executed on the tuned A3 event
    "zabf": ["--policy=a3", "--a3OffsetDb=6", "--a3TttMs=640", "--tteRank=1"] + NET,
}
HOLD = []  # the R2 pilot chose no escape hold-off (pilot_r2_eval.py)


def _ut_family():
    P = POLICIES
    P["ut"] = UT + HOLD
    P["ut40"] = UT + HOLD + ["--tteExitElevDeg=40"]
    # one-factor ablations of UT-CHO-40
    P["ut40_norank"] = ["--policy=tte-cho", "--prepareByTte=1", "--prepStrongest=2", "--escapeNearBestDb=3"] + HOLD + [
        "--tteExitElevDeg=40"]
    P["ut40_noprep"] = ["--policy=tte-cho", "--rankByUsable=1", "--escapeNearBestDb=3"] + HOLD + ["--tteExitElevDeg=40"]
    P["ut40_noreserve"] = ["--policy=tte-cho", "--rankByUsable=1", "--prepareByTte=1", "--escapeNearBestDb=3"] + HOLD + [
        "--tteExitElevDeg=40"]
    # the escape rule as the one factor: escape to the strongest candidate better than serving
    P["ut40_qesc"] = ["--policy=tte-cho", "--rankByUsable=1", "--prepareByTte=1", "--prepStrongest=2"] + ESC + HOLD + [
        "--tteExitElevDeg=40"]
    # the hold-off as the one factor
    P["ut40_althold"] = UT + ["--escapeHoldoffS=5", "--tteExitElevDeg=40"]
    # the v6 escape (ranked by the score over every candidate above Q, ties by quality, with hold-off)
    P["ut40_oldesc"] = ["--policy=tte-cho", "--rankByUsable=1", "--prepareByTte=1", "--prepStrongest=2",
                        "--escapeHoldoffS=5", "--tteExitElevDeg=40"]


_ut_family()

# N4: measured TLE-vs-precise errors of 9697 Starlink satellites (68.3rd
# percentile of |err| per axis); cross and radial as measured ratios of along
TLE_DIR = os.path.join(HERE, "tle")
TLE_FILE = os.path.join(TLE_DIR, "tle_error_samples.csv")
TLE_SUP = os.path.join(TLE_DIR, "starlink_supgp_2026-10-16.tle")
TLE_GP = os.path.join(TLE_DIR, "starlink_gp_2026-10-16.tle")
CROSS_PER_ALONG = 0.043
RADIAL_PER_ALONG = 0.070
EPH_AGE_BINS = [(0, 6, 2290.0), (6, 12, 4392.0), (12, 24, 8712.0), (24, 72, 17523.0)]
RUN = ["--numUes=100", "--simTime=600"]

MAIN_POL = ["a3", "a3t", "d2", "t1", "t1w1", "elev", "a3tf", "d2f", "t1f", "elevf", "a3tfh", "d2fh", "t1fh",
            "elevfh", "t1f40", "elevf40", "d2f4050", "eager", "tted2", "tted2f", "mst", "mst40", "mst40s",
            "mst_prep", "mst_prepr", "mst_ute", "ut", "ut40", "ut40_norank", "ut40_noprep", "ut40_noreserve",
            "ut40_qesc", "ut40_althold", "ut40_oldesc"]
URB_POL = ["a3tf", "d2f", "t1f", "elevf", "tted2f", "a3tfh", "d2fh", "t1fh", "elevfh", "t1f40", "elevf40",
           "d2f4050", "mst", "mst40", "mst40s", "mst_prep", "mst_prepr", "mst_ute", "ut", "ut40", "ut40_norank",
           "ut40_noprep", "ut40_noreserve", "ut40_qesc", "ut40_althold", "ut40_oldesc"]
RECENT = ["geof", "wangf", "zabf"]
MAIN_POL += RECENT
URB_POL += RECENT


def jobs_main(seeds):
    for ch, name, s in itertools.product(["los", "suburban"], MAIN_POL, seeds):
        yield f"main/{ch}/{name}", s, POLICIES[name] + [f"--channel={ch}"] + RUN
    # TTE-CHO = MST at zero error: per-seed identity check
    for ch, name, s in itertools.product(["los", "suburban"], ["tte", "tte40"], seeds):
        yield f"main/{ch}/{name}", s, POLICIES[name] + [f"--channel={ch}"] + RUN


def jobs_urban(seeds):
    for ch, name, s in itertools.product(["urban", "denseurban"], URB_POL, seeds):
        yield f"urban/{ch}/{name}", s, POLICIES[name] + [f"--channel={ch}"] + RUN


def jobs_pareto(seeds):
    fam = []
    for off, ttt in itertools.product([0, 2, 4, 6, 8, 10], [160, 640]):
        fam.append((f"a3_o{off}_t{ttt}", ["--policy=a3", f"--a3OffsetDb={off}", f"--a3TttMs={ttt}"]))
    # the geometric conditions carry the same measurement fallback as the rule,
    # as in the reference comparison, so outage is comparable
    for f in [30, 35, 40, 45, 50, 55]:
        fam.append((f"elev_f{f}", ["--policy=elev", f"--elevFloorDeg={f}"] + NET))
    for a, b in [(30, 40), (35, 45), (40, 50), (45, 55), (50, 60)]:
        fam.append((f"d2_{a}_{b}", ["--policy=d2", f"--d2ServElevDeg={a}", f"--d2CandElevDeg={b}"] + NET))
    for e in [25, 30, 35, 40, 45, 50]:
        fam.append((f"t1_e{e}", ["--policy=t1", f"--tteExitElevDeg={e}"] + NET))
        # T1 with the preparation restricted to the contour
        fam.append((f"t1c_e{e}", ["--policy=t1", f"--tteExitElevDeg={e}", "--prepAboveContour=1"] + NET))
        fam.append((f"mst_e{e}", POLICIES["mst"] + [f"--tteExitElevDeg={e}"]))
        fam.append((f"ut_e{e}", POLICIES["ut"] + [f"--tteExitElevDeg={e}"]))
    # T1 window length
    for w in [1, 2, 5, 20, 30]:
        fam.append((f"t1_w{w}", ["--policy=t1", f"--t1DurationS={w}"] + NET))
    for ch, (name, args), s in itertools.product(["los", "suburban"], fam, seeds):
        yield f"pareto/{ch}/{name}", s, args + [f"--channel={ch}"] + RUN


def eph_args(sa):
    # Along-track dominates public LEO ephemerides; cross-track and radial are
    # set to the ratios measured on 9697 Starlink element sets (N4).
    return [f"--ephSigmaAlongM={sa}", f"--ephSigmaCrossM={CROSS_PER_ALONG * sa}",
            f"--ephSigmaRadialM={RADIAL_PER_ALONG * sa}"]


def assumed(sa):
    return [f"--assumedSigmaAlongM={sa}", f"--assumedSigmaCrossM={CROSS_PER_ALONG * sa}",
            f"--assumedSigmaRadialM={RADIAL_PER_ALONG * sa}"]


EPH_POL = ["mst", "tte", "ut", "t1", "d2", "elev"]


def jobs_eph(seeds):
    for sa, pol, s in itertools.product([0, 1000, 2000, 5000, 10000, 20000], EPH_POL, seeds):
        yield f"eph/sa{sa}/{pol}", s, POLICIES[pol] + ["--channel=los"] + RUN + eph_args(sa)
    # network budgets for a quarter of the error it actually has
    for sa, s in itertools.product([5000, 10000, 20000], seeds):
        yield f"eph/sa{sa}/tte-under", s, POLICIES["tte"] + ["--channel=los"] + RUN + assumed(sa / 4) + eph_args(sa)


def emp_args(lo, hi, sa):
    return [f"--ephErrorFile={TLE_FILE}", f"--ephAgeMinH={lo}", f"--ephAgeMaxH={hi}"] + assumed(sa)


def jobs_ephemp(seeds):
    # N4: per-satellite offsets drawn from the measured samples of each age
    # bin; the network budgets for that bin's measured sigmas
    for (lo, hi, sa), pol, s in itertools.product(EPH_AGE_BINS, EPH_POL, seeds):
        yield f"ephemp/age{lo}-{hi}/{pol}", s, POLICIES[pol] + ["--channel=los"] + RUN + emp_args(lo, hi, sa)


def jobs_ephblk(seeds):
    # measured element-set errors (0-6 h) under blockage, every policy
    lo, hi, sa = EPH_AGE_BINS[0]
    # UT-CHO at the coverage edge ("ut") added
    pols = ["a3tf", "d2f", "t1f", "elevf", "mst", "mst40", "tte40", "ut40", "ut"]
    for ch, pol, s in itertools.product(["suburban", "urban"], pols, seeds):
        yield f"ephblk/{ch}/{pol}", s, POLICIES[pol] + [f"--channel={ch}"] + RUN + emp_args(lo, hi, sa)


def jobs_realtle(seeds):
    # real 53-degree Starlink shell; truth SupGP, network GP, both SGP4.
    # The network budgets the measured all-age sigmas.
    sa = 3450.0
    # MST + margin at the coverage edge ("tte"), the
    # single-factor control that separates the margin from UT-CHO's preparation
    pols = ["a3tf", "d2f", "t1f", "elevf", "mst", "mst40", "tte40", "ut", "ut40", "tte",
            "tted2f"] + RECENT
    tle = [f"--tleTruthFile={TLE_SUP}", f"--tleNetFile={TLE_GP}", "--altitudeKm=470"] + assumed(sa)
    for ch, pol, s in itertools.product(["los", "suburban"], pols, seeds):
        args = POLICIES[pol] + [f"--channel={ch}"] + RUN + tle
        if pol in ("mst", "mst40", "a3tf", "d2f", "t1f", "elevf", "tted2f") or pol in RECENT:
            args = [a for a in args if not a.startswith("--assumedSigma")]
        yield f"realtle/{ch}/{pol}", s, args


def jobs_zsweep(seeds):
    # Oct 2026: margin multiplier z on the real shell, to check the optimal margin
    # of the renewal-reward analysis (Proposition 2) against simulation. z = 2 is
    # realtle/{ch}/tte.
    sa = 3450.0
    tle = [f"--tleTruthFile={TLE_SUP}", f"--tleNetFile={TLE_GP}", "--altitudeKm=470"] + assumed(sa)
    for ch, z, s in itertools.product(["los", "suburban"], [0, 0.5, 1, 1.5, 3, 4, 5, 6], seeds):
        yield f"zsweep/{ch}/z{z}", s, POLICIES["tte"] + [f"--marginZ={z}", f"--channel={ch}"] + RUN + tle


def jobs_acc(seeds):
    for sa, s in itertools.product([0, 1000, 5000, 20000], seeds):
        yield f"acc/sa{sa}", s, ["--mode=tte-accuracy", "--ueProfile=static", "--numUes=60",
                                 "--accuracySamples=400"] + eph_args(sa)
    # the pre-R2 re-anchored model, for the record
    for sa, s in itertools.product([5000], seeds):
        yield f"acc/off{sa}", s, ["--mode=tte-accuracy", "--ueProfile=static", "--numUes=60",
                                  "--accuracySamples=400", "--ephModel=offset"] + eph_args(sa)


def jobs_lb(seeds):
    for name in ["a3t", "d2", "t1", "elev", "mst", "mst_prepr", "ut"]:
        for s in seeds:
            yield f"lb/{name}", s, POLICIES[name] + ["--channel=los", "--ueProfile=static",
                                                     "--numUes=40", "--simTime=600"]
    for s in seeds:
        yield "lb/mst_k16", s, POLICIES["mst"] + ["--maxCandidates=16", "--channel=los",
                                                  "--ueProfile=static", "--numUes=40", "--simTime=600"]


def jobs_trace(seeds):
    # per-terminal elevation traces (seed 1, terminal 0)
    for name in ["a3t", "t1", "mst", "mst40", "ut", "ut40"]:
        yield f"trace/{name}", 1, POLICIES[name] + ["--channel=los"] + RUN + ["--traceUe=0"]


SENS_POL = ["a3t", "d2", "t1", "elev", "mst", "ut", "ut40"]


def jobs_sens(seeds):
    # Geometric thresholds must sit above the coverage edge or the condition can
    # never fire; at e_min = 40 deg they are raised by 15 deg (D2 50/60, floor 45).
    shifted = {"d2": ["--policy=d2", "--d2ServElevDeg=50", "--d2CandElevDeg=60"],
               "elev": ["--policy=elev", "--elevFloorDeg=45"]}
    for me, name, s in itertools.product([10, 40], SENS_POL, seeds):
        args = shifted[name] if (me == 40 and name in shifted) else POLICIES[name]
        yield f"sens/minel{me}/{name}", s, args + [f"--minElevDeg={me}", "--channel=los"] + RUN
    # OneWeb-like shell: 1200 km, 87.9 deg, 648 = 18 x 36
    shell = ["--totalSats=648", "--numPlanes=18", "--phasing=1", "--altitudeKm=1200",
             "--inclinationDeg=87.9"]
    for name, s in itertools.product(SENS_POL, seeds):
        yield f"sens/leo1200/{name}", s, POLICIES[name] + shell + ["--channel=los"] + RUN


def jobs_density(seeds):
    shells = {"w792": ["--totalSats=792", "--numPlanes=36", "--phasing=17"],
              "w3168": ["--totalSats=3168", "--numPlanes=72", "--phasing=17"]}
    for (sh, sargs), ch, name, s in itertools.product(shells.items(), ["los", "suburban"],
                                                      ["d2f", "t1f", "mst40", "ut40"], seeds):
        yield f"density/{sh}/{ch}/{name}", s, POLICIES[name] + sargs + [f"--channel={ch}"] + RUN


def jobs_densgeo(seeds):
    # 29 Oct 2026: nearest-satellite rule (geodesic CHO) at three shell sizes, to test the
    # Voronoi-crossing rate of Proposition 2 in ns-3 beyond the reference shell
    shells = {"w792": ["--totalSats=792", "--numPlanes=36", "--phasing=17"],
              "w3168": ["--totalSats=3168", "--numPlanes=72", "--phasing=17"]}
    for (sh, sargs), s in itertools.product(shells.items(), seeds):
        yield f"density/{sh}/los/geof", s, POLICIES["geof"] + sargs + ["--channel=los"] + RUN


def jobs_losmodel(seeds):
    # LOS-state process and terminal environment, urban
    variants = {"mk50": ["--losModel=markov", "--losCorrM=50", "--losCorrDeg=10"],
                "mk200": ["--losModel=markov", "--losCorrM=200", "--losCorrDeg=10"],
                "mk50a3": ["--losModel=markov", "--losCorrM=50", "--losCorrDeg=3"],
                "noclass": ["--classEnv=0"]}
    pols = ["a3tf", "d2f", "t1f", "mst40", "ut40"]
    for (v, vargs), name, s in itertools.product(variants.items(), pols, seeds):
        yield f"losmodel/{v}/{name}", s, POLICIES[name] + vargs + ["--channel=urban"] + RUN


def jobs_reest(seeds):
    # a longer re-establishment (cell search, system information and
    # RRC exchange) for the policies whose failures it prices
    pols = ["d2", "t1", "d2f", "t1f", "mst40", "ut40"]
    for ch, name, s in itertools.product(["suburban", "urban"], pols, seeds):
        yield f"reest/{ch}/{name}", s, POLICIES[name] + ["--reestBaseMs=1000", f"--channel={ch}"] + RUN


# co-channel interference. Satellite s on carrier s mod K; every
# other covering satellite on that carrier reaches the terminal at its
# main-lobe power less the coupling (17.6 dB: first sidelobe of the TR 38.811
# aperture pattern; 10 dB: a harsher overlap case), full load. K = 0 is the
# main / urban campaign itself (regression byte-identical, ns-3 f9a107d66).
INTERF = {"k1c176": ["--interfReuse=1", "--interfCouplingDb=17.6"],
          "k3c176": ["--interfReuse=3", "--interfCouplingDb=17.6"],
          "k1c10": ["--interfReuse=1", "--interfCouplingDb=10"]}
INTERF_POL = ["a3tf", "d2f", "t1f", "elevf", "t1", "d2f4050", "mst", "mst40", "ut", "ut40"]


# the harsher 10 dB coupling only for the policies of the headline comparisons
INTERF_HI_POL = {"los": ["a3tf", "t1", "mst", "ut"], "urban": ["d2f", "d2f4050", "mst40", "ut40"]}


def jobs_interf(seeds):
    for ch, (v, vargs), name, s in itertools.product(["los", "urban"], INTERF.items(), INTERF_POL, seeds):
        if v == "k1c10" and name not in INTERF_HI_POL[ch]:
            continue
        yield f"interf/{ch}/{v}/{name}", s, POLICIES[name] + vargs + [f"--channel={ch}"] + RUN


# the margin
# budgets each serving satellite by the measured sigmas of its element set's
# age bin (tle/age_sigma_table.csv, tle_age_table.py), z = 2.326 (1 % one-sided)
# fixed a priori. ttefz: the same z with the measured all-age sigmas, the
# control that separates age-awareness from z and from the 3450 m budget.
AGE_TABLE = os.path.join(TLE_DIR, "age_sigma_table.csv")
Z_AGE = 2.326
ALLAGE = {"r": 266.4, "a": 3797.3, "c": 165.1}  # measured all-age sigmas (tle_age_table.py)
AGE = [f"--ageSigmaTable={AGE_TABLE}", "--marginByAge=1", f"--marginZ={Z_AGE}"]
ALLAGE_ARGS = [f"--assumedSigmaAlongM={ALLAGE['a']}", f"--assumedSigmaCrossM={ALLAGE['c']}",
               f"--assumedSigmaRadialM={ALLAGE['r']}"]


def jobs_agemargin(seeds):
    sa = 3450.0
    tle = [f"--tleTruthFile={TLE_SUP}", f"--tleNetFile={TLE_GP}", "--altitudeKm=470"]
    real = {"tteage": POLICIES["tte"] + AGE, "utage": POLICIES["ut"] + AGE,
            "ttefz": POLICIES["tte"] + ALLAGE_ARGS + [f"--marginZ={Z_AGE}"]}
    for ch, (name, args), s in itertools.product(["los", "suburban"], real.items(), seeds):
        yield f"realtle/{ch}/{name}", s, args + [f"--channel={ch}"] + RUN + tle
    gauss = {"mst": POLICIES["mst"], "tte": POLICIES["tte"] + ALLAGE_ARGS,
             "ttefz": POLICIES["tte"] + ALLAGE_ARGS + [f"--marginZ={Z_AGE}"],
             "tteage": POLICIES["tte"] + AGE, "ut": POLICIES["ut"] + ALLAGE_ARGS,
             "utage": POLICIES["ut"] + AGE}
    for (name, args), s in itertools.product(gauss.items(), seeds):
        # the error model reads the table; only the age-aware rules budget by it
        extra = [] if "--marginByAge=1" in args else [f"--ageSigmaTable={AGE_TABLE}"]
        yield f"ephage/los/{name}", s, args + extra + ["--ephGaussByAge=1", "--channel=los"] + RUN


EXPS = {"main": jobs_main, "urban": jobs_urban, "pareto": jobs_pareto, "eph": jobs_eph, "ephemp": jobs_ephemp,
        "ephblk": jobs_ephblk, "realtle": jobs_realtle, "acc": jobs_acc, "lb": jobs_lb, "trace": jobs_trace,
        "sens": jobs_sens, "density": jobs_density, "losmodel": jobs_losmodel, "reest": jobs_reest,
        "interf": jobs_interf, "agemargin": jobs_agemargin, "zsweep": jobs_zsweep,
        "densgeo": jobs_densgeo}


def run_one(tag, seed, args):
    out = os.path.join(RUNS, tag, f"seed{seed:02d}")
    done = os.path.join(out, "tte_accuracy.csv" if "--mode=tte-accuracy" in args else "summary.csv")
    if os.path.exists(done):
        return tag, seed, "skip", 0.0
    os.makedirs(out, exist_ok=True)
    if os.environ.get("PLACEHOLDER") == "1":
        # a second runner on the same tree: claim the run so that a runner
        # started earlier, which checks for the output file, skips it; the
        # simulator overwrites the placeholder when it finishes
        with open(done, "w") as fh:
            fh.write("placeholder\n")
    cmd = [BIN] + args + [f"--rngRun={seed}", f"--outputDir={out}"]
    t0 = time.time()
    with open(os.path.join(out, "stdout.txt"), "w") as fo:
        rc = subprocess.call(cmd, cwd=NS3, stdout=fo, stderr=subprocess.STDOUT)
    with open(os.path.join(out, "cmd.json"), "w") as fj:
        json.dump({"cmd": cmd, "rc": rc}, fj)
    return tag, seed, ("ok" if rc == 0 else f"rc={rc}"), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", default="all")
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=10)
    a = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"binary missing: {BIN}")
    # Provenance: which code produced these numbers.
    sha = subprocess.run(["git", "-C", NS3, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", NS3, "status", "--porcelain", "contrib/ntn-cho",
                            "contrib/ntn-constellation"], capture_output=True, text=True).stdout.strip()
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, "PROVENANCE.txt"), "a") as fp:
        fp.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} exp={a.exp} ns-3 {sha} "
                 f"{'DIRTY' if dirty else 'clean'} binary mtime "
                 f"{time.ctime(os.path.getmtime(BIN))}\n")
    seeds = list(range(1, a.seeds + 1))
    SEEDS_PER = {"pareto": 5, "trace": 1}
    names = list(EXPS) if a.exp == "all" else a.exp.split(",")
    jobs = [j for n in names for j in EXPS[n](list(range(1, min(a.seeds, SEEDS_PER.get(n, a.seeds)) + 1)))]
    if os.environ.get("SKIP_UT") == "1":
        # everything that does not depend on the pilot's hold-off choice
        jobs = [j for j in jobs if "/ut" not in j[0]]
    if os.environ.get("ONLY_TAG"):
        jobs = [j for j in jobs if os.environ["ONLY_TAG"] in j[0]]
    if os.environ.get("ONLY_UT") == "1":
        jobs = [j for j in jobs if "/ut" in j[0]]
    print(f"{len(jobs)} runs, {a.jobs} in parallel, ns-3 {sha[:9]}{' DIRTY' if dirty else ''}", flush=True)
    fails = 0
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(run_one, *j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            tag, seed, st, dt = f.result()
            if st not in ("ok", "skip"):
                fails += 1
            if st != "skip":
                print(f"[{i}/{len(jobs)}] {tag} seed{seed} {st} {dt:.0f}s", flush=True)
    print(f"done, {fails} failed", flush=True)


if __name__ == "__main__":
    main()
