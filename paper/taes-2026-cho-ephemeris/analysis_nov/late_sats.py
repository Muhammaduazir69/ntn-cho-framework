"""M4: which real-shell satellites carry the large along-track offsets behind the late exits
that no margin up to z = 6 covers. Pairs SupGP (truth) and GP (network) sets as the campaign
does (inclination 52.5-53.5 deg, mean altitude 455-485 km), propagates both with SGP4 to the
median truth epoch and reports the RTN offset, the two epochs, the GP set's age, mean-motion
derivative and B*. Run with the ntn-constellation venv (sgp4 package)."""
import sys, math, numpy as np
from sgp4.api import Satrec, jday
MU = 398600.4418; RE = 6378.137
def read(p):
    L = [l.rstrip() for l in open(p) if l.strip()]
    out = {}
    i = 0
    while i < len(L):
        if L[i].startswith("1 ") and L[i + 1].startswith("2 "):
            name, l1, l2 = "", L[i], L[i + 1]; i += 2
        else:
            name, l1, l2 = L[i].strip(), L[i + 1], L[i + 2]; i += 3
        out[int(l1[2:7])] = (name, l1, l2)
    return out
tru, net = read(sys.argv[1]), read(sys.argv[2])
rows = []
for cat in sorted(tru):
    if cat not in net: continue
    st = Satrec.twoline2rv(*tru[cat][1:]); sn = Satrec.twoline2rv(*net[cat][1:])
    n = st.no_kozai / 60.0
    a = (MU / n ** 2) ** (1 / 3); alt = a - RE; inc = math.degrees(st.inclo)
    if not (455 <= alt <= 485 and 52.5 <= inc <= 53.5): continue
    rows.append((cat, st, sn))
ep = np.median([r[1].jdsatepoch + r[1].jdsatepochF for r in rows])
big = []
for cat, st, sn in rows:
    e1, r1, v1 = st.sgp4(ep, 0.0); e2, r2, v2 = sn.sgp4(ep, 0.0)
    if e1 or e2: continue
    r1, v1, r2 = map(np.array, (r1, v1, r2)); d = r2 - r1
    R = r1 / np.linalg.norm(r1); N = np.cross(r1, v1); N /= np.linalg.norm(N); T = np.cross(N, R)
    da = float(d @ T)
    age_h = (ep - (sn.jdsatepoch + sn.jdsatepochF)) * 24
    big.append((cat, tru[cat][0], da, float(d @ R), age_h, (ep - (st.jdsatepoch + st.jdsatepochF)) * 24, sn.ndot, sn.bstar,
                float(np.linalg.norm(d))))
A = np.array([b[2] for b in big if b[8] < 1e5])
print(f"pairs in shell {len(big)}, kept (<100 km) {len(A)}, median epoch JD {ep:.4f}")
print(f"|along| 68.3 pct {np.percentile(np.abs(A),68.27)/1e0:.2f} km; share |along|>20 km {np.mean(np.abs(A)>20):.4f}; >30 km {np.mean(np.abs(A)>30):.4f}")
print("cat name along_km radial_km GPage_h SupGPage_h ndot bstar")
for b in sorted(big, key=lambda b: -abs(b[2])):
    if 20 < abs(b[2]) and b[8] < 1e5:
        print(f"{b[0]} {b[1]} {b[2]:.1f} {b[3]:.2f} {b[4]:.1f} {b[5]:.1f} {b[6]:.3e} {b[7]:.3e}")
