> Reviewer copy. Every span marked [redacted] removes internal wording; all other text is
> identical to the dated original whose SHA-256 values are logged in the run log, except the
> title, whose word "Pre-registration" is replaced by the article's term, pre-specified.

# Pre-specified analysis protocol: confirmatory runs [redacted]

**Written:** 10 Nov 2026, before any run of the groups below and before any result of
them was seen. Its SHA-256 is recorded in `sim/runs_nov/PROVENANCE.txt` before the first
run. Nothing above the "Amendments" section is edited after the first run; any change is
appended there with the time and the reason.

**Why.** [redacted]

Everything run on seeds 1-10 before today, and everything in group T below, is
**exploratory**. Groups C, O and S are **confirmatory**. The text labels both.

## 0. Binary, libraries, seeds

- Binary `sim/bin/p2-campaign-nov`, sha256 `ea8718be5399d46ea8bd520d13f87f28982d8955a1b9e1fdc14e73110a3a7af7`,
  built by `sim/ns3_patch/build.sh` from `sim/ns3_patch/orig/ntn-cho-tte-campaign.cc` (the exact
  source of the 29 Oct binary, sha256 `3b0e7840...`) plus two default-off options
  (`--softEdgeDeg/--softRollDbPerDeg`, `--tleStartOffsetH`). Frozen libraries in `sim/bin/lib`
  (`sim/bin/lib_manifest.sha256`). Every run records the binary sha256 in its `cmd.json`.
- Regression before any new run (`sim/regress_nov.py`): Walker rows must be byte-identical to
  `runs_r2`. Real-shell rows reproduce the event sequence; the libraries now carry the 30 Oct
  SGP4 time-precision fix, which moves positions by under a millimeter (logged, not hidden).
- Seeds: tuning T uses 1-10 (already exploratory). Confirmatory groups C and O use the fresh
  seeds **11-20**, never run before for any configuration. S uses 1-10 (it is a paired
  sensitivity against stored hard-edge rows, not a test of the rule).
- At most 3 jobs, nice 19, pid file `sim/run_nov.pid`.

## 1. Group T (exploratory tuning of the rivals, [redacted])

Seeds 1-10, 600 s, 100 terminals, same settings as `runs_r2` otherwise.
- TTE-D2 + net at D2 thresholds 25/35 and 30/40 deg (35/45 exists: `main|urban/*/tted2f`).
- Geodesic CHO + net at hysteresis 50, 100 and 200 km (10 km exists: `*/geof`).
- Environments: LOS, suburban, urban, dense urban.

**Selection rule, fixed now.** Per family and environment, choose the setting with the
lowest mean handover rate among the settings whose mean outage is at most 1.05 times the
lowest mean outage in that family and environment. Ties go to the reference setting. The
rule is applied by code (`run_nov.py select`) and its output is written to
`sim/runs_nov/rival_choice.json` before any confirmatory rival run starts.

## 2. Group C (confirmatory: stationary rates and fresh seeds, [redacted])

Seeds 11-20, `--simTime=1200`, 100 terminals.
- LOS and suburban: A3-tuned + net, D2 + net, T1 + net, residence priority + net (the A3-ranked
  mechanism; entropy-weighted is not rerun, since at seeds 1-10 it equals this row and A3-tuned),
  MST, UT-CHO, TTE-D2 + net (chosen setting), geodesic + net (chosen setting).
- LOS only: MST with K = 16 prepared candidates (`--maxCandidates=16`).
- Urban and dense urban: D2 + net, MST-40, UT-CHO-40, TTE-D2 + net (chosen), geodesic + net (chosen).

**Endpoints.**
- E1, stationary handover rate (primary for rates): for each terminal, the intervals between
  consecutive completed handovers after its first handover; the last interval of each terminal
  is right-censored at the end of the run. Per run, the mean interval is the Kaplan-Meier
  restricted mean pooled over terminals, and the rate is its inverse, per terminal-minute.
- E2, count rate: completed handovers per terminal-minute over [0, 1200] s, labeled as such.
- E3, outage share (unchanged definition).
- The exact-geometry Walker Monte Carlo is evaluated with the same estimators (E1 and E2) at
  the same horizon, so that ns-3 and the Monte Carlo are compared like with like.

**Families (seed-paired, exact two-sided Wilcoxon signed-rank, Holm within family).**
- F1 (headline rates, E1): UT-CHO against TTE-D2 + net, geodesic + net and residence priority
  + net, under LOS and suburban (6 tests). The rate cuts are reported as seed-paired mean
  relative differences with 95 % t intervals, whatever their sign.
- F2 (credit, E1): UT-CHO against MST (LOS, suburban) and against MST K = 16 (LOS) (3 tests).
- F3 (urban outage, E3): UT-CHO-40 against D2 + net, MST-40, TTE-D2 + net and geodesic + net,
  urban and dense urban (8 tests).
Whatever these show is what the abstract reports.

## 3. Group O (confirmatory: out-of-sample margin, [redacted])

No second element-set day is available locally (`sim/tle` holds only 16 Oct 2026), and no file
is downloaded for [redacted]. The hold-out is therefore in time and in seeds: the same
element-set pair, flown from **12 h after** the median truth epoch (`--tleStartOffsetH=12`), so
other satellites cross the box and the public sets are 12 h older, on seeds 11-20.
- MST and MST + margin (z = 2), the budget fixed at the all-age value sized at the epoch
  (sigma_A = 3.45 km, the same file), LOS and suburban, 600 s. 40 runs.
- Endpoints: late share and outage, MST + margin against MST, seed-paired (family F4, 4 tests,
  Holm).
- Prediction test: the outage model of Eq. (20), with every constant (J0, c_h, c_l, beta, k)
  fitted on the seeds 1-10 epoch sweep only, predicts the paired outage difference
  J(0) - J(2) of the hold-out window. The late tail Q(z) is taken (i) from the training
  window (pure transfer) and (ii) from the hold-out window's own element-set offsets (an input,
  not an outcome). Reported: predicted against observed difference with the observed 95 % CI,
  and the constant model (difference 0) as the null.
- In-sample, on the existing sweep (exploratory): two-fold cross-validation over seeds (fit on
  1-5, predict 6-10, and the reverse), held-out squared error of the model against the
  constant model, and a seed bootstrap interval for k.

## 4. Group S (sensitivity: soft coverage edge, [redacted])

Real shell at the epoch window, LOS, seeds 1-10, 600 s, `--softEdgeDeg=10
--softRollDbPerDeg=2`: the serving link survives below 25 deg with 2 dB extra loss per degree,
down to 15 deg. MST and MST + margin, plus MST at the hard edge rerun on this binary (to check
that the stored rows stand). 30 runs. Endpoints: outage of MST and the margin's outage saving,
soft against hard edge (paired by seed). Whatever the result, the late-decision claims are
scoped to the edge model in the text.

## 5. What will not be done

No change to UT-CHO, its escape, its reserved slots or z after this file. No seed above 20.
No further tuning on seeds 11-20.

## Amendments

### Amendment 1 (2026-11-10 02:56, before any confirmatory result was seen)
The rule of Section 1 picked the largest geodesic hysteresis of the grid (200 km) in all four
environments, with handover rate and LOS outage still falling at the grid's edge. A choice at
the boundary does not give the rival its best setting. Added: geodesic hysteresis 300, 500 and
1000 km on seeds 1-10 (exploratory, group T2), the same rule applied over the extended family
(10 to 1000 km) by code (`run_nov.py select2` -> `runs_nov/rival_choice2.json`), and
confirmatory runs of geodesic CHO at that choice (`conf/<env>/geof2`, seeds 11-20, 1200 s).
The extended choice is the geodesic row of the confirmatory table and of family F1/F3; the
200 km rows (`conf/<env>/geof`, as first registered) are reported beside it. TTE-D2's choice
(30/40 LOS and suburban, 25/35 urban and dense urban) is interior and stands. Nothing else changes.

### Amendment 2 (2026-11-10 04:37, after the soft-edge result of Section 4 was seen; post hoc)
[redacted] One roll-off value (2 dB/deg) cannot support a general statement about soft edges, and that the hard-edge MST + margin rows came from the earlier library build. Added, labeled post hoc: roll-off 5 and 10 dB/deg (MST and MST + margin, seeds 1-10) and MST + margin at the hard edge rerun on the frozen libraries (`soft/los/tte_hard`). Exploratory sensitivity only; no test family.

### Amendment 3 (2026-11-10 08:59, after [redacted]; exploratory additions)
(a) 2 x 2 of execution instant x target choice on the fresh seeds 11-20, 1200 s, LOS and suburban: D2 + net at 30/40 (condition, strongest) and T1 + net with a 1 s window (exit, strongest), beside the existing TTE-D2 + net 30/40 (condition, longest) and MST (exit, longest). Endpoint E1 rate; reported as cuts against D2 + net 30/40, paired, no family claim. 40 runs.
(b) Fixed against per-pass margin under inflated cross-track and radial error: Walker shell, LOS, Gaussian RTN error (sigma_A, sigma_C, sigma_R) = (2,2,2), (2,6,6), (2,10,10) km, seeds 11-20, 600 s. Per-pass arm: MST + margin, z = 2, budget = true sigmas. Fixed arm: MST with execMarginS equal to the per-pass arm's mean planned margin for that setting (computed from the per-pass runs by code before the fixed runs). Endpoints: late share and outage, paired. 60 runs.
(c) (09:01) Fig. 2 at the confirmatory horizon: MST, MST-40 and geodesic CHO at 10 km (nearest-satellite rule) on the 792- and 3168-satellite shells, and MST-40 on the reference shell, LOS, seeds 11-20, 1200 s (70 runs); read with E1 against the Monte Carlo at 1200 s. Exploratory figure data, no test family.
(d) (11:11, after the Fig. 2 points were seen; post hoc) MST with 32 prepared candidates on the 3168-satellite shell, seeds 11-20, 1200 s, LOS (10 runs): MST with K = 8 handed over 28 % above the Monte Carlo there; this checks whether the prepared-set size explains it.

### Amendment 4 (2026-11-10 11:46, after [redacted]; post hoc)
The fixed arm of Amendment 3(b) used the per-pass arm's mean margin, which a few grazing passes inflate (largest 273 s). Added: a fixed margin equal to the per-pass arm's median planned margin, same seeds 11-20 and settings (30 runs). A cap on sigma_tau inside the decision cannot be run: the decision margin is computed in the ntn-cho library, which [redacted] does not rebuild; the cap is evaluated offline only (how many passes it would change).

### Amendment 5 (2026-11-10 12:25; post hoc)
In the propagated copy model a radial error changes the semi-major axis and so drifts into along-track error at -1.5 n delta_R per second; over 600 s the along-track spread of the Amendment 3(b) runs grew from 2.0 km to 2.8, 6.2 and 10.0 km while the margin was budgeted at 2 km. Added: the same per-pass arm and a fixed arm at its median planned margin with the copy error held as a fixed RTN offset (--ephModel=offset, no drift), seeds 11-20, 600 s, (sigma_A, sigma_C, sigma_R) = (2,2,2), (2,6,6), (2,10,10) km (60 runs). Endpoint: late share, paired, Holm over the three settings.
