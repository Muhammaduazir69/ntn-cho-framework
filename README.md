<h1 align="center">NTN-CHO Framework</h1>

<p align="center"><strong>TTE-Aware 3GPP Release-17 Conditional Handover for 6G LEO Satellite Networks</strong></p>

<p align="center">
  <a href="https://www.nsnam.org"><img src="https://img.shields.io/badge/ns--3-3.43-blue.svg"/></a>
  <a href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"><img src="https://img.shields.io/badge/license-GPL--2.0-green.svg"/></a>
  <img src="https://img.shields.io/badge/3GPP-Rel--17-orange.svg"/>
  <img src="https://img.shields.io/badge/UE_classes-7%20(TR%2038.811)-purple.svg"/>
</p>

<p align="center">
  <img src="docs/architecture.png" alt="ntn-cho architecture" width="900"/>
</p>

---

## Why this module

LEO satellite handover is fundamentally an **orbital-dynamics** problem: a candidate cell that satisfies the 3GPP Release-17 conditional-handover (CHO) trigger at preparation can leave the user horizon before execution completes. This module integrates an **orbit-informed Time-to-Exit (TTE)** estimator into the standard Rel-17 CHO state machine, with the aim of suppressing ping-pong handovers and reducing the total handover count without breaking compliance with the existing RRC message set.

## Reported results

> **Provenance:** the figures below are reported in the associated manuscript
> (*under review* — see [Cite this work](#cite-this-work)). They come from a
> multi-seed sweep of the `ntn-cho-full-constellation` example over the four
> `--algorithm` settings; reproduce them by running that example across
> `--rngRun` seeds (see [INSTALL.md](INSTALL.md)). They are **not** asserted by
> the shipped `ntn-cho` unit-test suite, which validates the algorithm and
> state-machine logic.

| Metric (10-seed × 600-s, 66-sat Walker constellation) | TTE-aware (this work) | Location-only | A3 event |
|---|---|---|---|
| Handovers / seed | **135 ± 12** | 200 ± 40 | 463 ± 48 |
| Ping-pong rate | **0.0 %** | 57.1 % | 50.2 % |
| Success rate | **83.1 ± 4.6 %** | 98.7 ± 0.5 % | 70.3 % |
| TTE prediction cost | **O(log n) per candidate** | n/a | n/a |

## What it does

- 3GPP Rel-17 CHO state machine (`IDLE → PREPARE → MONITOR → EXEC`), `NtnChoAlgorithm`, with a `TriggerType` selector: **A3 event**, **location (condEventD1)**, **time-based**, the novel **TTE-aware** trigger, and a **THz-beam-quality** trigger
- `NtnTteEstimator` — Time-to-Exit via a coarse forward step plus **binary search** for the beam-exit instant (cached per candidate). Satellite positions come from the `satellite` module's SGP4 mobility — this module consumes that propagation rather than re-implementing it
- `NtnOrbitPredictor` — wraps the `satellite` module's `SatSGP4MobilityModel` + antenna-gain patterns to predict elevation / beam coverage
- `NtnRealisticMobilityHelper` — per-class UE motion for all **7 TR 38.811 §6.1.1.1 classes** (handheld-static / handheld-pedestrian / vehicular / HST / maritime / aviation / IoT-fixed), with a matching `NtnRealisticTrafficHelper`
- `NtnAiInterface` — *optional* ns3-ai bridge (observation + action structs covering serving/candidate geometry, recent-HO history, and a reward channel) for learning-based HO policies; requires the `ns3-ai` module. The C++ baselines run without it
- `ntn-cho-full-constellation` reference example over a Walker constellation (built with the `satellite` module), runnable with any of the `a3 / location / time / tte-aware` algorithms

## Install & run

See [**INSTALL.md**](INSTALL.md) for the full step-by-step guide (ns-3.43, SNS3 satellite module, build flags).

`ntn-cho` requires the SNS3 [`satellite`](https://github.com/sns3/sns3-satellite)
module (SGP4 mobility + antenna-gain patterns); the optional `NtnAiInterface`
additionally needs the [`ns3-ai`](https://github.com/hust-diangroup/ns3-ai)
module. Quick taste:

```bash
# from an ns-3.43 tree that already has contrib/satellite (and, for the
# AI path, contrib/ns3-ai):
git clone https://github.com/Muhammaduazir69/ntn-cho-framework.git contrib/ntn-cho
./ns3 configure --enable-examples --enable-tests
./ns3 build
./ns3 run "ntn-cho-full-constellation --algorithm=tte-aware --simTime=600 --rngRun=1"
```

## Documentation

- [INSTALL.md](INSTALL.md) — detailed setup + dependency notes
- [docs/architecture.png](docs/architecture.png) — module architecture
- Reference manuscript: *Time-to-Exit Conditional Handover for 6G LEO Satellite Networks* (under review)

## Cite this work

```bibtex
@misc{uzair2026ntncho,
  author = {Uzair, Muhammad},
  title  = {NTN-CHO Framework: TTE-Aware Conditional Handover for 6G LEO Satellite Networks},
  year   = {2026},
  url    = {https://github.com/Muhammaduazir69/ntn-cho-framework}
}
```

## Part of the ns3-ntn-toolkit

`ntn-cho` is one of the custom modules bundled in [**ns3-ntn-toolkit**](https://github.com/Muhammaduazir69/ns3-ntn-toolkit) — a pre-integrated ns-3.43 distribution for 6G NTN research. It also works as a standalone contrib module (see [INSTALL.md](INSTALL.md)). A few related modules (others are listed in the toolkit README):

| Module | Repo |
|---|---|
| Toolkit (umbrella) | [ns3-ntn-toolkit](https://github.com/Muhammaduazir69/ns3-ntn-toolkit) |
| **ntn-cho** | this repo |
| oran-ntn | [oran-ntn](https://github.com/Muhammaduazir69/oran-ntn) |
| thz-ntn | [ns3-thz-ntn](https://github.com/Muhammaduazir69/ns3-thz-ntn) |

## License

GPL-2.0-only — see [LICENSE](LICENSE).

## Acknowledgements

ns-3 core team · SNS3 maintainers · 3GPP RAN2 Rel-17 specifications.
