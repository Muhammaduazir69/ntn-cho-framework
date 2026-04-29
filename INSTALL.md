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

`ntn-cho` depends on `SatMobilityModel` for SGP4 orbit propagation:

```bash
cd contrib/
git clone https://github.com/sns3/sns3-satellite.git satellite
cd ..
```

> Size note: SNS3 + bundled TLE data is ~3.7 GB.

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
./ns3 run "ntn-realistic-mobility-demo --simTime=600 --outputDir=/tmp/mob_demo"
```

This spawns 14 UEs across all 7 TR 38.811 §6.1.1.1 classes, advances them
through 600 s of motion, and writes a per-step trajectory CSV.

### 5b. Full 66-satellite Walker-Star handover (~30 s wall clock)

```bash
./ns3 run "ntn-cho-full-constellation \
  --algorithm=tte-aware \
  --simTime=600 \
  --numUes=30 \
  --rngRun=1 \
  --outputDir=results/seed1/"
```

CSV outputs land under `results/seed1/`:
`ntn_handover_events.csv` · `tte_computations.csv` · `cho_state_log.csv` ·
`satellite_tracks.csv` · `ntn_measurements.csv` · `kpi_summary.txt`.

### 5c. Run the unit-test suite

```bash
./ns3 run "test-runner --suite=ntn-cho --verbose"
```

Expected: **3 / 3 passing**.

---

## 6. Reproduce the paper Monte-Carlo

A 10-seed × 600-s campaign with all four algorithms is included:

```bash
cd papers/sim_runs/
./run_mc_sweep.sh           # ~5 min total wall-clock on a modern CPU
python3 build_figures.py    # generates the publication figures
```

---

## 7. Common issues

**`ntn-cho not found` after configure**
You forgot the SNS3 satellite module. ntn-cho will not register if
`contrib/satellite/` is missing.

**`MmWaveAmc not found`**
You're missing the mmWave dependency. If you used the toolkit, this is
already there. Otherwise: `cd contrib && git clone
https://gitlab.com/cttc-lena/nr.git mmwave`.

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
