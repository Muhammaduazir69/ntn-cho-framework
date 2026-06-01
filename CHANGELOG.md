# Changelog

All notable changes to this module are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/) and this
project adheres to Semantic Versioning.

## [1.0.0] — 2026-04-23

### Added

- Initial public release of `ntn-cho`.

- **`NtnTteEstimator`** — Time-to-Exit estimator. For each candidate
  beam it propagates the serving/candidate satellite forward in time
  (via the SNS3 `satellite` module's SGP4 mobility) and uses a coarse
  forward step plus a **binary search** to locate the instant the beam
  geometry drops below the D1 distance threshold. The result is cached
  per candidate.

- **`NtnChoAlgorithm`** — 3GPP Release-17 conditional-handover state
  machine (`IDLE → PREPARE → MONITOR → EXEC`) with a `TriggerType`
  selector covering `TRIGGER_EVENT_A3`, `TRIGGER_LOCATION_D1`
  (condEventD1), `TRIGGER_TIME_BASED`, the proposed
  `TRIGGER_TTE_AWARE`, and `TRIGGER_THZ_BEAM_QUALITY`. Includes
  `SelectBaselineA3()` and a location-only selector for comparison
  against the TTE-aware path.

- **`NtnOrbitPredictor`** — thin wrapper around the `satellite`
  module's `SatSGP4MobilityModel` and antenna-gain-pattern container;
  predicts satellite sub-point / elevation / beam coverage over time.

- **`NtnMeasurementModel`** — per-UE RSRP / SINR / elevation / Doppler
  measurement generation used to drive the CHO triggers.

- **`NtnAiInterface`** — optional ns3-ai bridge exposing an
  observation struct (serving + candidate geometry, recent-HO history,
  reward channel) and an action struct, for learning-based HO policies.
  Requires the `ns3-ai` module; the C++ CHO baselines run without it.

- **`NtnChoHelper`** — convenience helper that wires the measurement
  model, TTE estimator, and CHO algorithm onto a node container.

- **`NtnRealisticMobilityHelper`** — per-class ground-UE motion
  generator for the seven 3GPP TR 38.811 §6.1.1.1 UE classes
  (`HANDHELD_STATIC`, `HANDHELD_PEDESTRIAN`, `VEHICULAR`,
  `HIGH_SPEED_TRAIN`, `MARITIME_VESSEL`, `AVIATION_COMMERCIAL`,
  `IOT_FIXED`).

- **`NtnRealisticTrafficHelper`** — per-class UDP traffic profiles so
  the data plane reflects the UE class (handheld bursts, HST video,
  IoT trickle, …).

- **Three example scripts** under `examples/`:
  - `ntn-cho-leo-basic` — minimal single-satellite CHO walkthrough.
  - `ntn-cho-full-constellation` — multi-satellite Walker constellation
    (built with the `satellite` module) running any of the
    `a3 / location / time / tte-aware` algorithms, writing per-event
    CSV/GeoJSON and a KPI summary.
  - `ntn-realistic-mobility-demo` — exercises the seven-class mobility
    helper.

- **One unit-test suite** (`test/ntn-cho-test-suite.cc`, suite name
  `ntn-cho`) with three `QUICK` cases: `NtnChoAlgorithmTestCase`,
  `NtnChoStateMachineTestCase`, and `NtnMeasurementModelTestCase`.

### Dependencies

- Requires the SNS3 [`satellite`](https://github.com/sns3/sns3-satellite)
  module (SGP4 mobility + antenna-gain patterns).
- `NtnAiInterface` additionally requires the `ns3-ai` module; it is only
  needed for the learning-based path.

### Notes

- SGP4 propagation and the Walker constellation geometry are provided
  by the `satellite` module — this module **consumes** them, it does not
  re-implement them.
- Aggregate KPIs quoted in the README (handover counts, ping-pong rate,
  success rate, confidence intervals) come from the associated
  manuscript (under review) and are **not** reproduced by the shipped
  unit-test suite, which validates the algorithm / state-machine logic.
