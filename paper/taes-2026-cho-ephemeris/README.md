# Code and data: CHO in LEO networks under ephemeris error (IEEE TAES, 2026)

Companion material for M. Uzair, "Conditional Handover in LEO Satellite Networks Under Ephemeris Error: Handover-Rate Analysis and Exit-Error Propagation," manuscript submitted to IEEE Transactions on Aerospace and Electronic Systems, 2026.

The handover rules (MST, exit-triggered prepared-candidate CHO and the tuned rivals) are in the `ntn-cho` module at the root of this repository. This folder holds everything that turns the module into the paper's numbers.

| Path | Contents |
|---|---|
| `ns3_patch/` | Campaign source of record (`ntn-cho-tte-campaign.cc`, also copied to `../../examples/`), the patch against the unpatched example, and `build.sh` |
| `run_campaign.py`, `run_nov.py` | Campaign drivers (seeds, shells, policies, ephemeris-error windows) |
| `run_summaries.tar` | Per-run summaries of every simulation in the paper (`cmd.json` plus `summary.csv` per run) |
| `analyze.py`, `analysis_nov/`, `analysis_rate/` | Statistics, the handover-rate model and its Walker Monte Carlo, the margin model and figure data |
| `analysis_protocol.md`, `rival_choice*.json` | Pre-specified analysis plan and the tuned settings of the rival rules |
| `verify_numbers.py`, `regress_nov.py` | Re-derive every number quoted in the paper from the summaries |
| `tle/` | Public Starlink element sets (CelesTrak GP and SupGP, 16 Oct 2026) and the element-set error tables |
| `MANIFEST.sha256`, `binary.sha256`, `lib_manifest.sha256`, `SOURCE_NOTE.txt` | Hashes of every file, of the campaign binary and of the linked libraries |

## Reproduce
1. Build ns-3.43 with this module in `contrib/`, then `ns3_patch/build.sh`.
2. `python3 run_campaign.py` and `python3 run_nov.py` regenerate the runs; or unpack `run_summaries.tar` to skip them.
3. `python3 analyze.py`, then the scripts in `analysis_nov/` and `analysis_rate/`, then `python3 verify_numbers.py`.

Raw per-event traces (about 2 GB) are not stored here; the scripts regenerate them, and they are available from the author on request.

License: as the repository (GPL-2.0-only for code). The element sets are public CelesTrak data.
Contact: Muhammad Uzair, muhammaduzairr69@gmail.com, ORCID 0009-0002-4104-2680.
