#!/usr/bin/env python3
"""Produce the Nov 2026 patched campaign example from orig/ (the exact source of the
29 Oct binary). Two default-off options; with both at their defaults the run is
unchanged (checked by regress_nov.py).
  --softEdgeDeg, --softRollDbPerDeg: the serving link survives below e_min
      down to e_min - softEdgeDeg, with an extra roll-off loss of softRollDbPerDeg per
      degree below e_min. Only the serving satellite; candidates, the predictor's
      contour and the late-decision test (true exit of e_min) are unchanged.
  --tleStartOffsetH: start the real-shell run this many hours after the
      median truth epoch (temporal hold-out on the same element-set pair)."""
import os
H = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H, "orig", "ntn-cho-tte-campaign.cc")).read()
E = [
 ("    double tleStartSpreadH = 1.0;     // start at the median truth epoch +- this\n",
  "    double tleStartSpreadH = 1.0;     // start at the median truth epoch +- this\n"
  "    double tleStartOffsetH = 0.0;     // Nov 2026 (M5): shift the start by this (h)\n"
  "    double softEdgeDeg = 0.0;         // Nov 2026 (M7): serving link kept this far below e_min (0 = hard edge)\n"
  "    double softRollDbPerDeg = 2.0;    //   extra roll-off loss per degree below e_min\n"),
 ("        if (!Covers(uePos, sp))\n        {\n            continue;\n        }\n        const double elev = ntngeo::GeodeticElevationDeg(uePos, sp);\n",
  "        bool soft = false;\n"
  "        if (!Covers(uePos, sp))\n        {\n"
  "            // M7: soft edge for the serving link only (off by default)\n"
  "            if (g_p.softEdgeDeg > 0.0 && static_cast<int32_t>(s) == u.serving &&\n"
  "                ntngeo::GeodeticElevationDeg(uePos, sp) >= g_p.minElevDeg - g_p.softEdgeDeg)\n"
  "            {\n                soft = true;\n            }\n"
  "            else\n            {\n                continue;\n            }\n        }\n"
  "        const double elev = ntngeo::GeodeticElevationDeg(uePos, sp);\n"),
 ("        ls.scanLossDb = ScanLossDb(uePos, sp);\n        // line-of-sight state",
  "        ls.scanLossDb = ScanLossDb(uePos, sp) -\n"
  "                        (soft ? g_p.softRollDbPerDeg * (g_p.minElevDeg - elev) : 0.0);\n"
  "        // line-of-sight state"),
 ("        tleStartUnix = epochs[epochs.size() / 2] +\n",
  "        tleStartUnix = epochs[epochs.size() / 2] + 3600.0 * p.tleStartOffsetH +\n"),
 ('    cmd.AddValue("tleStartSpreadH", "R2: start within +- this of the median truth epoch (h)", p.tleStartSpreadH);\n',
  '    cmd.AddValue("tleStartSpreadH", "R2: start within +- this of the median truth epoch (h)", p.tleStartSpreadH);\n'
  '    cmd.AddValue("tleStartOffsetH", "Nov 2026: shift the real-shell start by this (h)", p.tleStartOffsetH);\n'
  '    cmd.AddValue("softEdgeDeg", "Nov 2026: serving link kept this far below e_min (0 = hard edge)", p.softEdgeDeg);\n'
  '    cmd.AddValue("softRollDbPerDeg", "Nov 2026: roll-off below e_min (dB/deg)", p.softRollDbPerDeg);\n'),
]
for a, b in E:
    assert src.count(a) == 1, a[:60]
    src = src.replace(a, b)
open(os.path.join(H, "ntn-cho-tte-campaign.cc"), "w").write(src)
print("patched")
