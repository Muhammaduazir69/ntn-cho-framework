#!/usr/bin/env python3
"""Nov 2026 regression: rerun stored seed-1 rows with a given binary and the frozen
libraries in bin/lib, and compare events.csv, ue_summary.csv (and eph_offsets.csv)
byte for byte, and summary.csv except the wall_s column.
usage: regress_nov.py BINARY OUTDIR"""
import json, os, subprocess, sys, filecmp, csv, hashlib
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
NS3 = os.path.abspath(os.path.join(HERE, "../../../ns-3-dev"))
BIN = os.path.abspath(sys.argv[1]); OUT = os.path.abspath(sys.argv[2])
ROWS = ["main/los/ut", "realtle/los/mst", "urban/urban/d2f"] + sys.argv[3:]
env = dict(os.environ, LD_LIBRARY_PATH=os.path.join(HERE, "bin/lib"))
def one(tag):
    ref = os.path.join(HERE, "runs_r2", tag, "seed01")
    args = json.load(open(os.path.join(ref, "cmd.json")))["cmd"][1:]
    args = [a for a in args if not a.startswith("--outputDir")]
    out = os.path.join(OUT, tag, "seed01"); os.makedirs(out, exist_ok=True)
    cmd = ["nice", "-n", "19", BIN] + args + [f"--outputDir={out}"]
    with open(os.path.join(out, "stdout.txt"), "w") as fo:
        rc = subprocess.call(cmd, cwd=NS3, stdout=fo, stderr=subprocess.STDOUT, env=env)
    res = []
    for f in ["events.csv", "ue_summary.csv", "eph_offsets.csv"]:
        a, b = os.path.join(ref, f), os.path.join(out, f)
        if os.path.exists(a):
            res.append(f"{f}:{'SAME' if os.path.exists(b) and filecmp.cmp(a, b, shallow=False) else 'DIFF'}")
    def rows(p):
        r = list(csv.reader(open(p))); i = r[0].index("wall_s")
        return [x[:i] + x[i+1:] for x in r]
    sa, sb = os.path.join(ref, "summary.csv"), os.path.join(out, "summary.csv")
    res.append("summary-wall:" + ("SAME" if os.path.exists(sb) and rows(sa) == rows(sb) else "DIFF"))
    return f"{tag} rc={rc} " + " ".join(res)
sha = hashlib.sha256(open(BIN, "rb").read()).hexdigest()
print(f"binary {BIN} sha256 {sha}", flush=True)
with ThreadPoolExecutor(3) as ex:
    for r in ex.map(one, ROWS):
        print(r, flush=True)
