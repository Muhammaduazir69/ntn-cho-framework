"""Monte Carlo of the Walker shells of Fig. 2 read with estimator E1 at the ns-3 horizons."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rate_est as E
out = {}
P = {396: 36, 792: 36, 1584: 72, 3168: 72, 6336: 144}
for H in (600.0, 1200.0):
    for N, Pn in P.items():
        for eQ in (25, 40):
            if eQ == 40 and N == 396: continue
            out[f"exit_{N}_{eQ}_H{int(H)}"] = E.mc_exit(N, Pn, 17, 550, 53, eQ, H, 60, seed=21)
        out[f"near_{N}_H{int(H)}"] = E.mc_near(N, Pn, 17, 550, 53, H, 60, seed=23)
        print(H, N, {k: round(v["e1"], 4) for k, v in out.items() if f"_{N}_" in k and f"H{int(H)}" in k}, flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mc_fig.json"), "w"), indent=1)
