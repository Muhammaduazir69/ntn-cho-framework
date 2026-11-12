"""Statistics of the per-pass arm (margins, error-sigma correlation, escapes), an offline count of
passes a sigma cap would change, and the along-track drift of the propagated copy model against
the drift term of the exit-error lemma, sqrt(sA^2 + (1.5 n t sR)^2)."""
import csv, json, os, sys
import numpy as np
SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = {}
for k in ("r1", "r3", "r5"):
    M, E, S, esc = [], [], [], {}
    for arm in ("tte", "fixed"):
        n = 0
        for s in range(11, 21):
            for e in csv.DictReader(open(os.path.join(SIM, "runs_nov/perpass", k, arm, f"seed{s:02d}", "events.csv"))):
                if e["kind"] != "ho":
                    continue
                if e["trigger"] == "escape":
                    n += 1
                if arm == "tte" and e["trigger"] == "planned" and e["true_serv_exit_s"] not in ("", None):
                    m = float(e["margin_s"]); M.append(m)
                    E.append(abs(float(e["pred_serv_tte_s"]) - float(e["true_serv_exit_s"])))
                    S.append((m - 0.2) / 2.0)
        esc[arm] = n
    M, E, S = map(np.array, (M, E, S))
    along = 2000.0 * 1.041 / 7585.0  # eta sigma_A / v at 550 km
    out[k] = dict(mean=float(M.mean()), median=float(np.median(M)), max=float(M.max()), n_over20=int((M > 20).sum()),
                  share_below_mean=float(np.mean(M < M.mean())), corr=float(np.corrcoef(E, S)[0, 1]), n=len(M),
                  esc_per=esc["tte"], esc_fix=esc["fixed"], along_sigma=along,
                  share_capped_3x=float(np.mean(S > 3 * along)))
    print(k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in out[k].items()})
MU, RE = 3.986004418e14, 6378137.0
n = (MU / (RE + 550e3) ** 3) ** 0.5
for k, (sa, sr) in (("r1", (2000.0, 2000.0)), ("r3", (2000.0, 6000.0)), ("r5", (2000.0, 10000.0))):
    a0, aT = [], []
    for s in range(11, 21):
        for e in csv.DictReader(open(os.path.join(SIM, "runs_nov/perpass", k, "tte", f"seed{s:02d}", "eph_offsets.csv"))):
            a0.append(float(e["d0_a_m"])); aT.append(float(e["dT_a_m"]))
    out[k]["drift_sd0_km"] = float(np.std(a0)) / 1e3
    out[k]["drift_sdT_km"] = float(np.std(aT)) / 1e3
    out[k]["drift_pred_km"] = float(np.hypot(sa, 1.5 * n * 600.0 * sr)) / 1e3
    print(k, "drift sd0 %.2f sdT %.2f pred %.2f km" % (out[k]["drift_sd0_km"], out[k]["drift_sdT_km"], out[k]["drift_pred_km"]))
for k in ("r1", "r3", "r5"):
    E, S2, Z = [], [], []
    for s in range(11, 21):
        f = os.path.join(SIM, "runs_nov/perpass_off", k, "tte", f"seed{s:02d}", "events.csv")
        if not os.path.exists(f):
            break
        for e in csv.DictReader(open(f)):
            if e["kind"] == "ho" and e["trigger"] == "planned" and e["true_serv_exit_s"] not in ("", None):
                eps = float(e["pred_serv_tte_s"]) - float(e["true_serv_exit_s"]); sg = (float(e["margin_s"]) - 0.2) / 2.0
                E.append(abs(eps)); S2.append(sg); Z.append(eps / sg)
    if len(E) > 100:
        out[k]["off_corr"] = float(np.corrcoef(E, S2)[0, 1]); out[k]["off_marginmax"] = float(2 * max(S2) + 0.2)
        out[k]["off_zsd"] = float(np.std(Z))
        print(k, "offset corr %.3f" % out[k]["off_corr"])
json.dump(out, open(os.path.join(SIM, "analysis_nov", "perpass_stats.json"), "w"), indent=1)
