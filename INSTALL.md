# Install & run — ntn-cho

This guide walks you through installing the `ntn-cho` module on top of the
[ns3-ntn-toolkit](https://github.com/Muhammaduazir69/ns3-ntn-toolkit)
(or any vanilla ns-3.43 tree with the SNS3 `satellite` module).

---

## 1. System requirements

| Component | Version |
|---|---|
| OS | Linux (Ubuntu 22.04+ / Fedora 39+ recommended) |
| C++ compiler | gcc ≥ 11 or clang ≥ 14 |
| CMake | ≥ 3.24 |
| Python | ≥ 3.10 (3.13 supported) |
| ns-3 | **3.43** (other versions not supported on this branch) |
| Disk | ~6 GB after build (incl. SNS3 TLE data) |

---

## 2. Prerequisites

### 2a. Pull ns-3.43

If you don't already have it, the easiest path is to clone the umbrella
toolkit which bundles the patched ns-3.43 + mmwave + lte:

```bash
git clone https://github.com/Muhammaduazir69/ns3-ntn-toolkit.git
cd ns3-ntn-toolkit
```

Or pull a vanilla ns-3.43 from upstream — but you will then need to
manually apply the dual-connectivity LTE patches (see toolkit `INSTALL.md`).

### 2b. Add SNS3 satellite (REQUIRED)

`ntn-cho` depends on the SNS3 `satellite` module's `SatSGP4MobilityModel`
(SGP4 orbit propagation) and antenna-gain-pattern container:

```bash
cd contrib/
git clone https://github.com/sns3/sns3-satellite.git satellite
cd ..
```

> Size note: SNS3 + bundled TLE data is ~3.7 GB.

### 2c. Add ns3-ai (OPTIONAL — only for the learning-based path)

`NtnAiInterface` bridges to a Python agent via the `ns3-ai` module. The
C++ CHO baselines (A3 / location / time / tte-aware) build and run
**without** it; only add it if you want the AI interface:

```bash
cd contrib/
git clone https://github.com/hust-diangroup/ns3-ai.git
cd ..
```

---

## 3. Install the ntn-cho module

```bash
cd contrib/
git clone https://github.com/Muhammaduazir69/ntn-cho-framework.git ntn-cho
cd ..
```

---

## 4. Configure & build

```bash
./ns3 configure --enable-examples --enable-tests
./ns3 build ntn-cho
```

Verify the module is registered:

```bash
./ns3 show profile | grep ntn-cho
```

Expected output:

```
Modules to build: ... ntn-cho ...
```

---

## 5. Run an example

### 5a. Single-pass realistic mobility demo (~10 s wall clock)

```bash
./ns3 run "ntn-realistic-mobility-demo --simTime=600 --dt=1 --outputDir=/tmp/mob_demo"
```

This advances UEs spanning the seven TR 38.811 §6.1.1.1 classes through
600 s of motion and writes a per-step trajectory CSV
(`/tmp/mob_demo/mobility_trace.csv`).

### 5b. Full Walker-constellation handover (~30 s wall clock)

The constellation geometry (`--numPlanes`, `--satsPerPlane`) is built with
the `satellite` module; defaults give a multi-plane LEO Walker layout.

```bash
./ns3 run "ntn-cho-full-constellation \
  --algorithm=tte-aware \
  --simTime=600 \
  --numUes=30 \
  --rngRun=1 \
  --outputDir=results/seed1/"
```

Outputs land under `results/seed1/`:
`handover_events.csv` · `measurements.csv` · `tte_computations.csv` ·
`satellite_tracks.csv` · `ue_tracks.csv` · `kpi_timeseries.csv` ·
`kpi_summary.txt`, plus the GeoJSON layers `satellite_positions.geojson`,
`ue_positions.geojson`, `beam_footprints.geojson`, and
`handover_events.geojson`.

### 5c. Run the unit-test suite

```bash
./ns3 run "test-runner --suite=ntn-cho --verbose"
```

Expected: **3 / 3 passing**.

---

## 6. Reproduce the multi-seed comparison

The aggregate KPIs reported in the manuscript come from sweeping the
`ntn-cho-full-constellation` example over several seeds and the four
algorithms. There is no bundled harness — drive it with a shell loop:

```bash
for algo in a3 location time tte-aware; do
  for seed in $(seq 1 10); do
    ./ns3 run "ntn-cho-full-constellation \
      --algorithm=$algo --simTime=600 --numUes=30 \
      --rngRun=$seed --outputDir=results/$algo/seed$seed/"
  done
done

# Aggregate the per-run outputs into a comparison table + figures:
python3 tools/analyze_results.py --datadir results/ --output results/figures/
```

`tools/generate_cho_visualizations.py` turns the per-run CSV/GeoJSON into
figures if you want them.

---

## 7. Common issues

**`ntn-cho not found` after configure**
You forgot the SNS3 satellite module. ntn-cho will not register if
`contrib/satellite/` is missing.

**`Ns3AiMsgInterface` / `ns3-ai` errors when building**
Only the optional `NtnAiInterface` needs `ns3-ai` (step 2c). If you don't
want the AI path, you can build the rest of the module without that
module present.

**Build cache filtering modules**
If `ntn-cho` is silently excluded, run a clean configure:
`./ns3 configure --enable-modules='' --enable-tests --enable-examples`
(empty `--enable-modules=''` re-enables all autoloaded modules).

---

## 8. Uninstall

```bash
rm -rf contrib/ntn-cho
./ns3 configure --enable-examples
./ns3 build
```

---

## 9. Citing

See the [README](README.md#cite-this-work) for the BibTeX entry.
