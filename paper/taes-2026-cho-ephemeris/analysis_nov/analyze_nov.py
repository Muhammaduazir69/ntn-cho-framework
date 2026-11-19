#!/usr/bin/env python3
"""November 2026 campaign: every new number of the manuscript (macros_nov.tex), the
confirmatory table (table_conf.tex), the rival-tuning and margin tables of the supplement,
and Figs. fig_rate.pdf and fig_margin.pdf (figures_v8/). Pre-registration:
analysis_protocol.md. An independent recheck is verify_nov.py.

    python3 analysis_nov/analyze_nov.py
"""
import csv, json, math, os, sys
import statistics as st
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
PAPER = os.path.dirname(SIM)
sys.path.insert(0, HERE)
import rate_est as RE_  # noqa: E402
import margin_model as MM  # noqa: E402

NOV = os.path.join(SIM, "runs_nov")
R2 = os.path.join(SIM, "runs_r2")
FIG = os.path.join(PAPER, "figures_v8")
os.makedirs(FIG, exist_ok=True)
M = {}
REPORT = {}


def put(k, v):
    assert k.isalpha(), k
    M[k] = v


def fmt(x, d):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "--"
    s = f"{x:.{d}f}"
    if s.startswith("-") and float(s) == 0.0:
        s = s[1:]
    if s.startswith("-"):
        s = "\\ensuremath{-}" + s[1:]  # math minus in text and in tables
    return s


def tci(v):
    v = [x for x in v if x is not None and not math.isnan(x)]
    n = len(v)
    if n < 2:
        return (v[0] if v else float("nan")), 0.0, n
    return st.mean(v), stats.t.ppf(0.975, n - 1) * st.stdev(v) / math.sqrt(n), n


def wil(d):
    d = [x for x in d if abs(x) > 1e-12]
    if len(d) < 1:
        return 1.0
    return float(stats.wilcoxon(d, method="exact").pvalue) if len(d) > 0 else 1.0


def holm(ps):
    m = len(ps)
    o = sorted(range(m), key=lambda i: ps[i])
    adj, run = [0.0] * m, 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (m - r) * ps[i]))
        adj[i] = run
    return adj


_cache = {}


def runs(path):
    """per-seed metrics of a run directory (E1, E2, outage, SNR)"""
    if path in _cache:
        return _cache[path]
    out = {}
    if os.path.isdir(path):
        for s in sorted(os.listdir(path)):
            d = os.path.join(path, s)
            f = os.path.join(d, "summary.csv")
            if s.startswith("seed") and os.path.exists(f) and "ho_per_ue_min" in open(f).readline():
                r = RE_.ns3_run(d)
                r["late"] = late_share(d)
                out[r["seed"]] = r
    _cache[path] = out
    return out


def late_share(d):
    n = l = 0
    for e in csv.DictReader(open(os.path.join(d, "events.csv"))):
        if e["kind"] != "ho":
            continue
        n += 1
        if e["true_serv_exit_s"] not in ("", None) and float(e["true_serv_exit_s"]) <= 0:
            l += 1
    return 100.0 * l / max(1, n)


def mean_of(path, key):
    R = runs(path)
    return tci([r[key] for r in R.values()])


def paired(pa, pb, key, rel=False):
    a, b = runs(pa), runs(pb)
    k = sorted(set(a) & set(b))
    if not k:
        return None
    if rel:  # cut of a relative to b, % per seed
        d = [100.0 * (1 - a[x][key] / b[x][key]) for x in k]
    else:
        d = [a[x][key] - b[x][key] for x in k]
    m, h, n = tci(d)
    return dict(d=m, ci=h, n=n, p=wil(d), lower=sum(1 for x in d if x < 0), higher=sum(1 for x in d if x > 0))


C = lambda env, pol: os.path.join(NOV, "conf", env, pol)
ENVS4 = ["los", "suburban", "urban", "denseurban"]
CH = {"los": "Los", "suburban": "Sub", "urban": "Urb", "denseurban": "Dense"}
PN = {"a3tf": "Atuned", "d2f": "Dtwof", "t1f": "Tonef", "zabf": "Rank", "mst": "Mst", "ut": "Ut",
      "tted2f": "Tted", "d2f_30_40": "Dtwoft", "t1f_w1": "Tonefw", "geof2": "Geo", "geof": "GeoTwoh", "geof10": "Near", "mst_k16": "MstK", "mst40": "Mstforty", "ut40": "Utforty"}


# ---------------------------------------------------------------- rates (C2, M3, M10)
def rate_block():
    mc = json.load(open(os.path.join(HERE, "mc_rate_like.json")))
    put("NvMcExitEone", fmt(mc["exit_e25_H1200"]["e1"], 3))
    put("NvMcExitEtwo", fmt(mc["exit_e25_H1200"]["e2"], 3))
    put("NvMcExitEtwoSix", fmt(mc["exit_e25_H600"]["e2"], 3))
    put("NvMcExitEoneSix", fmt(mc["exit_e25_H600"]["e1"], 3))
    put("NvMcExitStat", fmt(mc["exit_e25_H2400"]["e1"], 3))
    put("NvMcFortyEone", fmt(mc["exit_e40_H1200"]["e1"], 3))
    put("NvMcNearEone", fmt(mc["near_H1200"]["e1"], 3))
    put("NvMcRatio", fmt(mc["near_H1200"]["e1"] / mc["exit_e25_H1200"]["e1"], 1))
    # analysis (stationary) at the reference shell
    an = json.load(open(os.path.join(SIM, "analysis_rate", "analysis_rate.json")))["1584_25_550"]
    put("NvRatioAsym", fmt(8 / math.pi ** 1.5 * math.sqrt(an["Nvis"]), 2))
    # exploratory 600 s counts (what the earlier version reported) against the same runs read with E1
    for pol in ["mst", "ut", "geof"]:
        m2, _, _ = mean_of(os.path.join(R2, "main", "los", pol), "e2")
        m1, _, _ = mean_of(os.path.join(R2, "main", "los", pol), "e1")
        put(f"NvSix{PN[pol]}Etwo", fmt(m2, 3))
        put(f"NvSix{PN[pol]}Eone", fmt(m1, 3))
    # confirmatory: every policy, every environment
    for env in CH:
        for pol in PN:
            p = C(env, pol)
            if not runs(p):
                continue
            for key, d in (("e1", 3), ("e2", 3), ("out", 3), ("snr", 1)):
                m, h, n = mean_of(p, key)
                put(f"Nv{CH[env]}{PN[pol]}{ {'e1': 'Eone', 'e2': 'Etwo', 'out': 'Out', 'snr': 'Snr'}[key]}",
                    fmt(m, d if key != "out" or env in ("los", "suburban") else 2))
                if key in ("e1", "out"):
                    put(f"Nv{CH[env]}{PN[pol]}{ {'e1': 'Eone', 'out': 'Out'}[key]}Ci", fmt(h, 3))
    if runs(C("los", "mst")) and runs(C("los", "geof10")) and runs(C("los", "geof2")):
        put("NvRatioGeoMst", fmt(mean_of(C("los", "geof10"), "e1")[0] / mean_of(C("los", "mst"), "e1")[0], 1))
        put("NvLosMstVsMcPct", fmt(100 * (mean_of(C("los", "mst"), "e1")[0] / mc["exit_e25_H1200"]["e1"] - 1), 1))
        put("NvLosGeoVsMcPct", fmt(100 * (mean_of(C("los", "geof10"), "e1")[0] / mc["near_H1200"]["e1"] - 1), 1))
        gaps = [mean_of(C("los", b), "snr")[0] - mean_of(C("los", "ut"), "snr")[0] for b in ("tted2f", "geof2", "zabf")]
        put("NvLosSnrGapMin", fmt(min(gaps), 1)); put("NvLosSnrGapMax", fmt(max(gaps), 1))
        put("NvLosGeoSnrGain", fmt(mean_of(C("los", "geof2"), "snr")[0] - mean_of(C("los", "ut"), "snr")[0], 1))
        put("NvLosUtVsMcPct", fmt(100 * (mean_of(C("los", "ut"), "e1")[0] / mc["exit_e25_H1200"]["e1"] - 1), 1))
    # families
    fams = {}
    F1 = [(env, riv) for env in ("los", "suburban") for riv in ("tted2f", "geof2", "zabf")]
    fams["FOne"] = [("ut", env, riv, "e1", True) for env, riv in F1]
    fams["FTwo"] = [("ut", "los", "mst", "e1", True), ("ut", "suburban", "mst", "e1", True),
                  ("ut", "los", "mst_k16", "e1", True)]
    fams["FThree"] = [("ut40", env, riv, "out", False) for env in ("urban", "denseurban")
                  for riv in ("d2f", "mst40", "tted2f", "geof2")]
    res = {}
    for fn, items in fams.items():
        rows = []
        for a, env, b, key, rel in items:
            r = paired(C(env, a), C(env, b), key, rel)
            if r:
                rows.append(((a, env, b), r))
        if len(rows) != len(items):
            continue
        adj = holm([r["p"] for _, r in rows])
        for ((a, env, b), r), h in zip(rows, adj):
            r["holm"] = h
            tag = f"Nv{CH[env]}{PN[a]}Vs{PN[b]}"
            put(tag, fmt(r["d"], 1 if rel else 3))
            put(tag + "Ci", fmt(r["ci"], 1 if rel else 3))
            put(tag + "P", ("$p_{\\mathrm{H}}<0.05$" if h < 0.05 else f"$p_{{\\mathrm{{H}}}}={h:.2f}$"))
            put(tag + "Raw", f"$p={r['p']:.3f}$")
            res[(a, env, b)] = r
        put(f"NvFam{fn}Sig", str(sum(1 for (_, r) in rows if r["holm"] < 0.05)))
        put(f"NvFam{fn}Size", str(len(rows)))
    for env in ("los", "suburban"):
        cuts = [res[("ut", env, b)]["d"] for b in ("tted2f", "geof2", "zabf") if ("ut", env, b) in res]
        if cuts:
            put(f"Nv{CH[env]}CutMin", fmt(min(cuts), 1))
            put(f"Nv{CH[env]}CutMax", fmt(max(cuts), 1))
    REPORT["families"] = {f"{a}|{e}|{b}": r for (a, e, b), r in res.items()}
    return res


def table_conf(res):
    """Table III: confirmatory comparison, seeds 11-20, 1200 s"""
    ch = json.load(open(os.path.join(NOV, "rival_choice.json")))["choice"]
    if os.path.exists(os.path.join(NOV, "rival_choice2.json")):
        ch.update({k: v for k, v in json.load(open(os.path.join(NOV, "rival_choice2.json")))["choice"].items()})
    lab = {"tted2f": lambda e: "TTE-D2 + net, " + ch[f"tted2f/{e}"].split("_", 1)[1].replace("_", "/") + "$^\\circ$",
           "geof": lambda e: "Geodesic + net, " + ch[f"geof/{e}"].split("_h")[1] + " km"}

    def cell(env, pol, key, d):
        R = runs(C(env, pol))
        if not R:
            return "--"
        m, h, _ = mean_of(C(env, pol), key)
        return fmt(m, d)

    def dcell(env, a, b, key, rel, d):
        r = res.get((a, env, b))
        if not r:
            return ""
        star = "$^{\\ast}$" if r["holm"] < 0.05 else ""
        if rel and r["ci"] > 0:
            return fmt(r["d"], d) + "\\,$\\pm$\\," + fmt(r["ci"], d) + star
        return fmt(r["d"], d) + star

    rows = [("A3-tuned + net", "a3tf", "a3tf"), ("D2 35/45$^\\circ$ + net", "d2f", "d2f"), ("T1 + net", "t1f", "t1f"),
            ("A3 + residence rank$^{a}$", "zabf", None), ("TTE-D2 + net$^{b}$", "tted2f", "tted2f"),
            ("Geodesic + net$^{c}$", "geof2", "geof2"), ("MST", "mst", None),
            ("MST-40", None, "mst40"), ("UT-CHO$^{d}$", "ut", None), ("UT-CHO-40", None, "ut40")]
    L = []
    L.append("\\begin{table*}[!t]")
    L.append("\\caption{Confirmatory Comparison on Fresh Seeds (11 to 20, 1200 s, Exact Ephemeris)}")
    L.append("\\label{tab:conf}\\centering\\footnotesize\\setlength{\\tabcolsep}{3.0pt}")
    L.append("\\begin{tabular}{@{}lrrrr@{\\hspace{8pt}}rrr@{\\hspace{8pt}}rrr@{\\hspace{8pt}}rrr@{}}\\toprule")
    L.append(" & \\multicolumn{4}{c}{Line of sight} & \\multicolumn{3}{c}{Suburban} & \\multicolumn{3}{c}{Urban} & "
             "\\multicolumn{3}{c}{Dense urban}\\\\")
    L.append("\\cmidrule(lr){2-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}\\cmidrule(l){12-14}")
    L.append("Policy & HO & $\\Delta$HO$_{\\mathrm{UT}}$\\,\\% & Out.\\,\\% & SNR & HO & $\\Delta$HO$_{\\mathrm{UT}}$\\,\\% & Out.\\,\\% & HO & Out.\\,\\% & "
             "$\\Delta$Out. & HO & Out.\\,\\% & $\\Delta$Out.\\\\\\midrule")
    for name, pl, pu in rows:
        c = [name]
        for env in ("los", "suburban"):
            if pl and runs(C(env, pl)):
                c.append(cell(env, pl, "e1", 3))
                c.append(dcell(env, "ut", pl, "e1", True, 1) if pl not in ("ut",) else "")
                c.append(cell(env, pl, "out", 3))
                if env == "los":
                    c.append(cell(env, pl, "snr", 2))
            else:
                c += ["--", "", "--"] + (["--"] if env == "los" else [])
        for env in ("urban", "denseurban"):
            if pu and runs(C(env, pu)):
                c.append(cell(env, pu, "e1", 2))
                c.append(cell(env, pu, "out", 2))
                c.append(dcell(env, "ut40", pu, "out", False, 2) if pu != "ut40" else "")
            else:
                c += ["--", "--", ""]
        if name.startswith("UT-CHO$"):
            L.append("\\midrule")
        L.append(" & ".join(c) + "\\\\")
    L.append("\\bottomrule\\end{tabular}\\\\[2pt]")
    L.append("\\parbox{\\textwidth}{\\scriptsize Ten fresh seeds, 100 terminals. HO: handovers per terminal-minute, inverse "
             "Kaplan-Meier mean residence; $\\Delta$HO$_{\\mathrm{UT}}$: rate cut of UT-CHO against the row (pre-specified family; under LOS UT-CHO equals MST with $K=16$), seed-paired mean $\\pm$ 95\\,\\% $t$ interval, \\%; $^{d}$under LOS the same handovers as MST with $K=16$ in every seed; "
             "$\\Delta$Out.: outage of UT-CHO-40 minus the row, points; $^{\\ast}$: Holm $p<0.05$ in its declared family. "
             "Rule-family settings (Table~S16), LOS/suburban/urban/dense urban: TTE-D2 "
             + "; ".join(ch[f"tted2f/{e}"].split("_", 1)[1].replace("_", "/") for e in ENVS4) + "$^\\circ$; geodesic "
             + "; ".join(ch[f"geof/{e}"].split("_h")[1] for e in ENVS4) + "~km. Families representative of $^{a}$"
             "\\cite{zaboub2026multi,been2025priority,wang2026mdp}, $^{b}$\\cite{noreen2026location}, "
             "$^{c}$\\cite{lee2026geodesic}. SNR in dB.}")
    L.append("\\end{table*}")
    open(os.path.join(PAPER, "table_conf.tex"), "w").write("\n".join(L) + "\n")


def table_tune():
    ch = json.load(open(os.path.join(NOV, "rival_choice.json")))
    f2 = os.path.join(NOV, "rival_choice2.json")
    if os.path.exists(f2):  # Amendment 1: extended geodesic family
        c2 = json.load(open(f2))
        for env in ENVS4:
            ch["table"][f"geof/{env}"] = c2["table"][f"geof/{env}"]
            ch["choice"][f"geof/{env}"] = c2["choice"][f"geof/{env}"]
    L = ["\\begin{table}[!t]", "\\caption{Tuning of the Rule Families (Exploratory Seeds 1 to 10, 600 s)}", "\\label{tab:tune}",
         "\\centering\\scriptsize\\setlength{\\tabcolsep}{2.2pt}",
         "\\begin{tabular}{@{}lrrrrrrrr@{}}\\toprule",
         " & \\multicolumn{2}{c}{LOS} & \\multicolumn{2}{c}{Suburban} & \\multicolumn{2}{c}{Urban} & "
         "\\multicolumn{2}{c}{Dense}\\\\",
         "Setting & HO & Out. & HO & Out. & HO & Out. & HO & Out.\\\\\\midrule"]
    for fam in ("tted2f", "geof"):
        sets = list(ch["table"][f"{fam}/los"].keys())
        for s in sets:
            nm = ("TTE-D2 " + s.split("_", 1)[1].replace("_", "/")) if fam == "tted2f" else ("Geodesic " + s.split("_h")[1] + " km")
            c = [nm]
            for env in ("los", "suburban", "urban", "denseurban"):
                ho, out = ch["table"][f"{fam}/{env}"][s]
                b = ch["choice"][f"{fam}/{env}"] == s
                w = (lambda x: "\\textbf{" + x + "}") if b else (lambda x: x)
                c += [w(fmt(ho, 3)), w(fmt(out, 3 if env in ("los", "suburban") else 2))]
            L.append(" & ".join(c) + "\\\\")
        L.append("\\midrule" if fam == "tted2f" else "")
    L.append("\\bottomrule\\end{tabular}\\\\[2pt]")
    L.append("\\parbox{\\columnwidth}{\\scriptsize Means over ten seeds; HO is the count per terminal-minute over 600~s "
             "(the selection statistic fixed in advance), outage in \\%. Bold: the setting chosen per environment by "
             "the rule fixed in advance (lowest HO among settings with outage within 5\\,\\% of the family's lowest).}")
    L.append("\\end{table}")
    open(os.path.join(PAPER, "table_tune.tex"), "w").write("\n".join(L) + "\n")
    for k, v in ch["choice"].items():
        fam, env = k.split("/")
        put(f"NvChoice{ {'tted2f': 'Tted', 'geof': 'Geo'}[fam]}{CH[env]}",
            v.split("_", 1)[1].replace("_", "/") if fam == "tted2f" else v.split("_h")[1])


# ---------------------------------------------------------------- margin (C1, M4, M5, M7)
def margin_block():
    out = {}
    for ch in ("los", "suburban"):
        R = MM.load(ch, range(1, 11))
        P = MM.fit(R)
        c = MM.cv(ch)
        b = MM.boot_k(ch)
        zm = MM.zmeans(R)
        Cn = CH[ch]
        put(f"NvProp{Cn}Ch", fmt(P["c_h"] * 1e3, 0))
        put(f"NvProp{Cn}Ct", fmt(P["ct"], 2))
        put(f"NvProp{Cn}Jzero", fmt(100 * P["J0"], 3))
        put(f"NvProp{Cn}Nu", fmt(60 * P["nu0"], 3))
        put(f"NvProp{Cn}Sigma", fmt(P["sbar"], 2))
        put(f"NvProp{Cn}Zstar", fmt(MM.zstar_gauss(P), 2))
        put(f"NvProp{Cn}Zemp", fmt(MM.zstar_numeric(P), 1))
        put(f"NvProp{Cn}Zexact", fmt(MM.zstar_gauss_exact(P), 2))
        put(f"NvProp{Cn}Czero", fmt((MM.TX - MM.TS / 2) / P["sbar"], 2))
        put(f"NvProp{Cn}KfitLo", fmt(b["k_lo"] * 1e3, 1))
        put(f"NvProp{Cn}KfitHi", fmt(b["k_hi"] * 1e3, 1))
        put(f"NvProp{Cn}KfitNonpos", fmt(100 * b["p_nonpos"], 0))
        put(f"NvCv{Cn}Model", fmt(c["prof_model"], 4))
        put(f"NvCv{Cn}Const", fmt(c["prof_const"], 4))
        put(f"NvCv{Cn}Ratio", fmt(c["prof_const"] / c["prof_model"], 1))
        put(f"NvCv{Cn}Step", fmt(c["prof_step"], 4))
        put(f"NvProp{Cn}ZempLo", fmt(b["z_lo"], 1))
        put(f"NvProp{Cn}ZempHi", fmt(b["z_hi"], 1))
        sg = [x for r in R for x in r["sig"]]
        put(f"NvProp{Cn}SigLo", fmt(float(np.percentile(sg, 5)), 2))
        put(f"NvProp{Cn}SigHi", fmt(float(np.percentile(sg, 95)), 2))
        put(f"NvProp{Cn}ER", fmt(1.0 / P["nu0"], 0))
        put(f"NvProp{Cn}KfitMid", fmt(P["k_fit"] * 1e3, 1))
        put(f"NvProp{Cn}NuSec", fmt(P["nu0"] * 1e3, 1))
        put(f"NvProp{Cn}Rise", fmt(100 * P["nu0"] * P["sbar"] * 6, 1))
        put(f"NvProp{Cn}LateSide", fmt(100 * float(np.mean(P["eps"] > 0)), 0))
        put(f"NvProp{Cn}QSix", fmt(100 * MM.q(P, 6), 1))
        # plateau: J(2) - min over z in [2, 6] in percentage points (model)
        Js = [MM.J(P, z) for z in np.linspace(2, 6, 41)]
        put(f"NvProp{Cn}Plateau", fmt(MM.J(P, 2) - min(Js), 4))
        out[ch] = dict(P=P, zm=zm, cv=c, boot=b)
        # simulated outage flat from z = 1.5: paired tests
        for za, zb, nm in ((1.5, 2, "OneHalfTwo"), (2, 5, "TwoFive"), (0, 2, "ZeroTwo")):
            pa, pb = MM.zdir(ch, za), MM.zdir(ch, zb)
            r = paired(pa, pb, "out")
            put(f"NvZs{Cn}{nm}", fmt(r["d"], 3))
            put(f"NvZs{Cn}{nm}Ci", fmt(r["ci"], 3))
            put(f"NvZs{Cn}{nm}P", f"$p={r['p']:.3f}$")
        # absolute saving of the margin per terminal and 10 min, z = 0 -> 2
        r = paired(MM.zdir(ch, 0), MM.zdir(ch, 2), "out")
        put(f"NvSave{Cn}Sec", fmt(r["d"] / 100 * 600, 2))
    # held-out window (+12 h, seeds 11-20)
    for ch in ("los", "suburban"):
        Cn = CH[ch]
        pm, pt = os.path.join(NOV, "oos", ch, "mst"), os.path.join(NOV, "oos", ch, "tte")
        if len(runs(pm)) < 10 or len(runs(pt)) < 10:
            continue
        put(f"NvOos{Cn}MstLate", fmt(mean_of(pm, "late")[0], 2))
        put(f"NvOos{Cn}TteLate", fmt(mean_of(pt, "late")[0], 2))
        rl = paired(pm, pt, "late")
        ro = paired(pm, pt, "out")
        put(f"NvOos{Cn}LateD", fmt(rl["d"], 2))
        put(f"NvOos{Cn}LateDCi", fmt(rl["ci"], 2))
        put(f"NvOos{Cn}LateLower", str(rl["higher"]))
        put(f"NvOos{Cn}OutD", fmt(ro["d"], 3))
        put(f"NvOos{Cn}OutDCi", fmt(ro["ci"], 3))
        put(f"NvOos{Cn}MstOut", fmt(mean_of(pm, "out")[0], 3))
        put(f"NvOos{Cn}TteOut", fmt(mean_of(pt, "out")[0], 3))
        out[ch]["oos_raw"] = (rl, ro)
        # prediction with constants from the exploratory epoch sweep
        P = out[ch]["P"]
        Ro = []
        for pol, z in (("mst", 0), ("tte", 2)):
            for s in range(11, 21):
                d = os.path.join(NOV, "oos", ch, pol, f"seed{s:02d}")
                if os.path.exists(os.path.join(d, "summary.csv")):
                    Ro.append(MM.load_run(d, z))
        # exit errors of the held-out window from its element-set offsets (Lemma 2), an input
        eps_new = np.array([e for r in Ro for e in r["eps"]])
        pred_tr = MM.J(P, 0) - MM.J(P, 2)
        pred_new = MM.J(P, 0, eps_new) - MM.J(P, 2, eps_new)
        put(f"NvOos{Cn}PredTr", fmt(pred_tr, 3))
        put(f"NvOos{Cn}PredNew", fmt(pred_new, 3))
        # late per exit, predicted from the offsets against observed
        for z, pol, nm in ((0, "mst", "Zero"), (2, "tte", "Two")):
            obs = [r for r in Ro if r["z"] == z]
            put(f"NvOos{Cn}LateExitObs{nm}", fmt(100 * sum(r["late"] for r in obs) / sum(r["exits"] for r in obs), 2))
            put(f"NvOos{Cn}LateExitPred{nm}", fmt(100 * MM.q(P, z, eps_new), 2))
        out[ch]["oos"] = dict(pred_tr=pred_tr, pred_new=pred_new, eps_new=eps_new)
        # late shares without the heavy-tail satellites
    # four families F4: late and outage, mst vs tte, LOS and suburban (Holm)
    rows = [(ch, k, out[ch]["oos_raw"][i]) for ch in ("los", "suburban") if "oos_raw" in out[ch]
            for i, k in enumerate(("late", "out"))]
    if len(rows) == 4:
        adj = holm([r["p"] for _, _, r in rows])
        for (ch, k, r), h in zip(rows, adj):
            put(f"NvOos{CH[ch]}{ {'late': 'Late', 'out': 'Out'}[k]}H",
                "$p_{\\mathrm{H}}<0.05$" if h < 0.05 else f"$p_{{\\mathrm{{H}}}}={h:.2f}$")
    # soft edge (M7), real shell LOS, seeds 1-10
    S = os.path.join(NOV, "soft", "los")
    if all(len(runs(os.path.join(S, p))) == 10 for p in ("mst", "tte", "mst_hard")):
        hard_t = os.path.join(R2, "realtle", "los", "tte")
        put("NvSoftMstOut", fmt(mean_of(os.path.join(S, "mst"), "out")[0], 3))
        put("NvSoftTteOut", fmt(mean_of(os.path.join(S, "tte"), "out")[0], 3))
        put("NvHardMstOut", fmt(mean_of(os.path.join(R2, "realtle", "los", "mst"), "out")[0], 3))
        put("NvHardTteOut", fmt(mean_of(hard_t, "out")[0], 3))
        rs = paired(os.path.join(S, "mst"), os.path.join(S, "tte"), "out")
        rh = paired(os.path.join(R2, "realtle", "los", "mst"), hard_t, "out")  # stored pair, same libraries
        put("NvSoftSave", fmt(rs["d"], 3))
        put("NvSoftSaveCi", fmt(rs["ci"], 3))
        put("NvSoftSaveP", f"$p={rs['p']:.3f}$")
        put("NvHardSave", fmt(rh["d"], 3))
        put("NvHardSaveCi", fmt(rh["ci"], 3))
        put("NvSoftMstLate", fmt(mean_of(os.path.join(S, "mst"), "late")[0], 2))
        put("NvSoftTteLate", fmt(mean_of(os.path.join(S, "tte"), "late")[0], 2))
        # stored hard-edge rows reproduce? event sequences
        same = 0
        for s in range(1, 11):
            a = os.path.join(S, "mst_hard", f"seed{s:02d}", "events.csv")
            b = os.path.join(R2, "realtle", "los", "mst", f"seed{s:02d}", "events.csv")
            ka = [l.split(",")[:6] for l in open(a)]
            kb = [l.split(",")[:6] for l in open(b)]
            same += ka == kb
        put("NvHardSameSeq", str(same))
        rmh = paired(os.path.join(S, "mst_hard"), os.path.join(R2, "realtle", "los", "mst"), "out")
        put("NvHardStoredDiff", fmt(abs(rmh["d"]), 4))
    # Amendment 2: same-library hard edge and the roll-off sweep
    th = os.path.join(S, "tte_hard")
    if len(runs(th)) == 10 and len(runs(os.path.join(S, "mst_hard"))) == 10:
        put("NvSoftHardMstOut", fmt(mean_of(os.path.join(S, "mst_hard"), "out")[0], 3))
        r = paired(os.path.join(S, "mst_hard"), th, "out")
        put("NvSoftHardSave", fmt(r["d"], 3))
        put("NvSoftHardSaveCi", fmt(r["ci"], 3))
        parts = []
        for ro in (5, 10):
            a_, b_ = os.path.join(S, f"mst_r{ro}"), os.path.join(S, f"tte_r{ro}")
            if len(runs(a_)) == 10 and len(runs(b_)) == 10:
                rr = paired(a_, b_, "out")
                put(f"NvSoftR{ {5: 'Five', 10: 'Ten'}[ro]}Save", fmt(rr["d"], 3))
                put(f"NvSoftR{ {5: 'Five', 10: 'Ten'}[ro]}SaveCi", fmt(rr["ci"], 3))
                put(f"NvSoftR{ {5: 'Five', 10: 'Ten'}[ro]}MstOut", fmt(mean_of(a_, "out")[0], 3))
                parts.append(f"{fmt(rr['d'], 3)}\\,$\\pm$\\,{fmt(rr['ci'], 3)} at {ro}")
        put("NvSoftRollSentence", ("; at steeper roll-offs it was " + " and ".join(parts) + "~dB per degree") if parts else "")
    # cluster analysis (M4)
    lc = json.load(open(os.path.join(HERE, "late_clusters.json")))
    for tag, nm in (("realtle/los/mst", "LosMst"), ("realtle/los/tte", "LosTte"), ("realtle/los/ut", "LosUt"),
                    ("realtle/suburban/mst", "SubMst"), ("realtle/suburban/tte", "SubTte"),
                    ("zsweep/los/z6", "LosZsix")):
        r = lc[tag]
        put(f"NvCl{nm}Lo", fmt(r["ci"][0], 1))
        put(f"NvCl{nm}Hi", fmt(r["ci"][1], 1))
        put(f"NvCl{nm}Wo", fmt(r["late_wo"], 2))
        put(f"NvCl{nm}WoTen", fmt(r["wo10"], 2))
        put(f"NvCl{nm}WoThirty", fmt(r["wo30"], 2))
        put(f"NvCl{nm}Pairs", str(r["pairs"]))
        put(f"NvCl{nm}Sats", str(r["sats"]))
        put(f"NvCl{nm}BigShare", fmt(100 * r["share_big_ev"], 1))
    return out


def table_margin(mo):
    L = ["\\begin{table}[!t]", "\\caption{Constants of the Outage Model \\eqref{eq:jz} (Exploratory Sweep, Seeds 1 to 10) and Its Transfer Test}",
         "\\label{tab:margin}", "\\centering\\scriptsize\\setlength{\\tabcolsep}{3pt}",
         "\\begin{tabular}{@{}llrr@{}}\\toprule", "Quantity & How obtained & LOS & Suburban\\\\\\midrule"]
    def row(name, how, f):
        L.append(f"{name} & {how} & {f(mo['los'])} & {f(mo['suburban'])}\\\\")
    row("$c_h$ (ms)", "mean interruption, from the runs", lambda m: fmt(1e3 * m["P"]["c_h"], 0))
    row("$\\nu_0$ (per min)", "exit rate at $z=0$, from the runs", lambda m: fmt(60 * m["P"]["nu0"], 3))
    row("$k$ (s$^{-1}$)", "$=\\nu_0$, derived", lambda m: fmt(m["P"]["nu0"], 4))
    row("$\\sigma_\\tau$ (s)", "median of the per-pass $\\sigma_\\tau$", lambda m: fmt(m["P"]["sbar"], 2))
    row("$q(z)$", "offsets via Lemma~\\ref{lem:sigma}, $>\\Delta+z\\sigma_\\tau$", lambda m: f"{100*MM.q(m['P'],0):.1f}/{100*MM.q(m['P'],2):.1f}/{100*MM.q(m['P'],6):.1f}\\,\\%")
    row("$\\tilde{c}_\\ell$ (s)", "within-seed regression, fitted", lambda m: fmt(m["P"]["ct"], 2))
    row("$\\beta$", "extra handovers per late exit, fitted", lambda m: fmt(m["P"]["beta"], 2))
    row("$J_0$ (\\%)", "remainder of outage, fitted", lambda m: fmt(100 * m["P"]["J0"], 3))
    row("$k$, free fit ($10^{-3}$s$^{-1}$)", "seed bootstrap 95\\,\\% interval",
        lambda m: f"{1e3*m['boot']['k_lo']:.1f} to {1e3*m['boot']['k_hi']:.1f}")
    row("$z^\\star$", "\\eqref{eq:zstar}, Gaussian tail", lambda m: fmt(MM.zstar_gauss(m["P"]), 2))
    row("$z^\\star$, exact root", "Section~\\ref{sec:proofs}", lambda m: fmt(MM.zstar_gauss_exact(m["P"]), 2))
    row("$\\arg\\min J$", "offset tail, $z\\le8$ (bootstrap 95\\,\\%)", lambda m: f"{MM.zstar_numeric(m['P']):.1f} ({m['boot']['z_lo']:.1f}--{m['boot']['z_hi']:.1f})")
    row("held-out error (pts)", "two-fold seed CV, model", lambda m: fmt(m["cv"]["prof_model"], 4))
    row("held-out error (pts)", "two-fold seed CV, two-level step", lambda m: fmt(m["cv"]["prof_step"], 4))
    row("held-out error (pts)", "two-fold seed CV, no change", lambda m: fmt(m["cv"]["prof_const"], 4))
    g = lambda k: M.get(k, "--").replace("\\ensuremath{-}", "$-$")
    L.append(f"saving, held-out window (pts) & predicted, training-window tail & {g('NvOosLosPredTr')} & {g('NvOosSubPredTr')}\\\\")
    L.append(f"saving, held-out window (pts) & predicted, window's own offsets & {g('NvOosLosPredNew')} & {g('NvOosSubPredNew')}\\\\")
    L.append(f"saving, held-out window (pts) & observed, seed-paired 95\\,\\% & {g('NvOosLosOutD')}\\,$\\pm$\\,{g('NvOosLosOutDCi')} & {g('NvOosSubOutD')}\\,$\\pm$\\,{g('NvOosSubOutDCi')}\\\\")
    L.append("\\bottomrule\\end{tabular}\\\\[2pt]")
    L.append("\\parbox{\\columnwidth}{\\scriptsize Exits are planned decisions plus escapes taken after the true exit. "
             "$q$ at $z=0/2/6$; $\\Delta=T_{\\mathrm{x}}-T_s/2=0.15$~s. Held-out error: root-mean-square error of the predicted change in outage from $z=0$, "
             "fitted on seeds 1 to 5 and tested on 6 to 10 and the reverse; the two-level step takes the training means at $z=0$ and $z>0$; both folds come from the same element-set window.}")
    L.append("\\end{table}")
    open(os.path.join(PAPER, "table_margin.tex"), "w").write("\n".join(L) + "\n")


def window_block():
    """Element-set error of the two real-shell windows"""
    import subprocess
    w = json.load(open(os.path.join(HERE, "window_stats.json")))
    e, h = w["epoch"], w["held"]
    put("NvWinEpochSigLo", fmt(e["p68_lo"], 1)); put("NvWinEpochSigHi", fmt(e["p68_hi"], 1))
    put("NvWinHeldSigLo", fmt(h["p68_lo"], 1)); put("NvWinHeldSigHi", fmt(h["p68_hi"], 1))
    put("NvWinHeldRatio", fmt(0.5 * (h["p68_lo"] + h["p68_hi"]) / 3.45, 1))
    put("NvWinEpochMean", fmt(e["mean"], 1)); put("NvWinHeldMean", fmt(h["mean"], 1))
    put("NvWinEpochLateSide", fmt(100 * e["lateside"], 0)); put("NvWinHeldLateSide", fmt(100 * h["lateside"], 0))
    put("NvWinEpochBig", fmt(100 * e["big"], 1)); put("NvWinHeldBig", fmt(100 * h["big"], 1))
    put("NvWinMirTwo", fmt(w["mirrored"]["z2"], 1)); put("NvWinMirSix", fmt(w["mirrored"]["z6"], 1))
    put("NvWinMirZ", fmt(w["mirrored"]["z_for_0p1"], 0))
    put("NvWinArmZero", fmt(w["margin_arm_only"]["z0"], 2))


def attr_block():
    """2 x 2 of execution instant x target choice (Amendment 3a), E1, fresh seeds"""
    rows = []
    for env in ("los", "suburban"):
        base = C(env, "d2f_30_40")
        if len(runs(base)) < 10 or len(runs(C(env, "t1f_w1"))) < 10:
            return
        for pol, nm in (("tted2f", "Tted"), ("t1f_w1", "Tonefw"), ("mst", "Mst")):
            r = paired(C(env, pol), base, "e1", True)
            put(f"NvAttr{CH[env]}{nm}", fmt(r["d"], 1)); put(f"NvAttr{CH[env]}{nm}Ci", fmt(r["ci"], 1))
        # MST (K = 8) against the three rule families: cut ranges
        rr = [paired(C(env, "mst"), C(env, b), "e1", True) for b in ("tted2f", "geof2", "zabf")]
        cuts = [r["d"] for r in rr]
        put(f"NvMst{CH[env]}CutMin", fmt(min(cuts), 1)); put(f"NvMst{CH[env]}CutMax", fmt(max(cuts), 1))
        put(f"NvMst{CH[env]}CutCiMax", fmt(max(r["ci"] for r in rr), 1))
    L = ["\\begin{table}[!t]", "\\caption{Where the Saving Comes From: Execution Instant and Target (Post Hoc; Fresh Seeds, 1200~s)}",
         "\\label{tab:attr}\\centering\\footnotesize\\setlength{\\tabcolsep}{3pt}",
         "\\begin{tabular}{@{}llcc@{}}\\toprule",
         "Instant & Target & LOS & Suburban\\\\\\midrule"]
    def c(env, pol):
        if pol == "d2f_30_40":
            return fmt(mean_of(C(env, pol), "e1")[0], 3) + " (ref.)"
        r = paired(C(env, pol), C(env, "d2f_30_40"), "e1", True)
        return fmt(mean_of(C(env, pol), "e1")[0], 3) + f" ({fmt(r['d'], 1)}\\,\\%)"
    for inst, tgt, pol in (("D2 30/40$^\\circ$", "strongest", "d2f_30_40"), ("D2 30/40$^\\circ$", "longest", "tted2f"),
                           ("predicted exit$^{a}$", "strongest", "t1f_w1"), ("predicted exit", "longest", "mst")):
        L.append(f"{inst} & {tgt} & {c('los', pol)} & {c('suburban', pol)}\\\\")
    L.append("\\bottomrule\\end{tabular}\\\\[2pt]")
    L.append("\\parbox{\\columnwidth}{\\scriptsize Handovers per terminal-minute (E1) and seed-paired cut against "
             "D2 30/40$^\\circ$ + net; all with the safety net. $^{a}$T1, 1~s window ending at the predicted exit.}")
    L.append("\\end{table}")
    open(os.path.join(PAPER, "table_attr.tex"), "w").write("\n".join(L) + "\n")


def perpass_block():
    """fixed against per-pass margin at equal mean margin (Amendment 3b), Walker LOS, seeds 11-20"""
    nm = {"r1": "One", "r3": "Three", "r5": "Five"}
    PPMED = []
    for k, n in nm.items():
        a, b = os.path.join(NOV, "perpass", k, "tte"), os.path.join(NOV, "perpass", k, "fixed")
        if len(runs(a)) < 10 or len(runs(b)) < 10:
            return
        put(f"NvPp{n}LatePer", fmt(mean_of(a, "late")[0], 2)); put(f"NvPp{n}LateFix", fmt(mean_of(b, "late")[0], 2))
        rl = paired(b, a, "late"); ro = paired(b, a, "out")
        put(f"NvPp{n}LateD", fmt(rl["d"], 2)); put(f"NvPp{n}LateDCi", fmt(rl["ci"], 2)); put(f"NvPp{n}LateP", f"$p={rl['p']:.3f}$")
        put(f"NvPp{n}OutD", fmt(ro["d"], 3)); put(f"NvPp{n}OutDCi", fmt(ro["ci"], 3)); put(f"NvPp{n}OutP", f"$p={ro['p']:.3f}$")
        ms = [float(e["margin_s"]) for sd in sorted(os.listdir(a)) for e in csv.DictReader(open(os.path.join(a, sd, "events.csv")))
              if e["kind"] == "ho" and e["trigger"] == "planned"]
        put(f"NvPp{n}MarginMean", fmt(st.mean(ms), 2))
        put(f"NvPp{n}MarginLo", fmt(float(np.percentile(ms, 5)), 2)); put(f"NvPp{n}MarginHi", fmt(float(np.percentile(ms, 95)), 2))
        c = os.path.join(NOV, "perpass", k, "fixedmed")
        if len(runs(c)) == 10:
            put(f"NvPp{n}LateMed", fmt(mean_of(c, "late")[0], 2))
            rm = paired(c, a, "late")
            put(f"NvPp{n}LateMedD", fmt(rm["d"], 2)); put(f"NvPp{n}LateMedDCi", fmt(rm["ci"], 2))
            put(f"NvPp{n}LateMedP", f"$p={rm['p']:.3f}$")
            put(f"NvPp{n}MarginMed", fmt(st.median(ms), 2))
            PPMED.append(rm["p"])
    if len(PPMED) == 3:  # Holm over the three error settings (drift case, per pass vs fixed median)
        q = sorted(PPMED)
        put("NvPpMedHolmMin", f"{min(1.0, max(q[0] * 3, 0)):.3f}")
    ps = os.path.join(HERE, "perpass_stats.json")
    if os.path.exists(ps):
        J = json.load(open(ps))
        for k, n in nm.items():
            put(f"NvPp{n}Corr", fmt(J[k]["corr"], 2)); put(f"NvPp{n}EscPer", str(J[k]["esc_per"]))
            put(f"NvPp{n}EscFix", str(J[k]["esc_fix"])); put(f"NvPp{n}MarginMax", fmt(J[k]["max"], 0))
            if "drift_sdT_km" in J[k]:
                put(f"NvDrift{n}Obs", fmt(J[k]["drift_sdT_km"], 1)); put(f"NvDrift{n}Pred", fmt(J[k]["drift_pred_km"], 1))
            if "off_corr" in J[k]:
                put(f"NvPoff{n}Corr", fmt(J[k]["off_corr"], 2))
                put(f"NvPoff{n}ZSd", fmt(J[k]["off_zsd"], 2))
    # Amendment 5: the same test with the copy error held as a fixed offset (no drift)
    lm = []
    for k, n in nm.items():
        a, c = os.path.join(NOV, "perpass_off", k, "tte"), os.path.join(NOV, "perpass_off", k, "fixedmed")
        if len(runs(a)) < 10 or len(runs(c)) < 10:
            return
        put(f"NvPoff{n}LatePer", fmt(mean_of(a, "late")[0], 2)); put(f"NvPoff{n}LateMed", fmt(mean_of(c, "late")[0], 2))
        rm = paired(c, a, "late")
        lm.append((n, rm))
        put(f"NvPoff{n}LateMedD", fmt(rm["d"], 2)); put(f"NvPoff{n}LateMedDCi", fmt(rm["ci"], 2))
        put(f"NvPoff{n}LateMedP", f"$p={rm['p']:.3f}$")
        ro = paired(c, a, "out")
        put(f"NvPoff{n}OutD", fmt(ro["d"], 3)); put(f"NvPoff{n}OutDCi", fmt(ro["ci"], 3))
    ps_ = sorted([r["p"] for _, r in lm])
    adj = [min(1.0, max(p * (3 - i) for i, p in enumerate(ps_[: j + 1]))) for j in range(3)]
    put("NvPoffHolmMin", f"{adj[0]:.3f}")


def fig_rate2():
    """Fig. 2: (a) rate against shell size, Monte Carlo and ns-3 both read with E1
    at the confirmatory 1200 s horizon (ns-3 seeds 11-20); (b) rate against latitude, Monte Carlo
    only. Grayscale-safe: line styles and marker shapes carry the series; one serif font."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(SIM, "analysis_rate"))
    import rate_model as RM
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "pdf.fonttype": 42, "font.size": 8, "axes.labelsize": 8,
                         "legend.fontsize": 7, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": 0.6})
    K = "#1f2933"; G = "#6b7785"
    mcf = json.load(open(os.path.join(HERE, "mc_fig.json")))
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.05))
    Ns = np.geomspace(396, 6336, 13)
    sizes = [396, 792, 1584, 3168, 6336]
    a = ax[0]
    sty = {25: ("-", "o"), 40: ("-.", "^")}
    for eQ in (25, 40):
        ls, mk = sty[eQ]
        Hc = [60 / RM.poisson_chain(int(n), 550, 53, eQ, K=4000) for n in Ns]
        a.plot(Ns, Hc, color=K, lw=1.0, ls=ls)
        a.plot([], [], color=K, lw=1.0, ls=ls, marker=mk, ms=4, mfc="white", label=rf"exit rule, $e_Q={eQ}^\circ$")
        xs = [n for n in sizes if f"exit_{n}_{eQ}_H1200" in mcf]
        a.plot(xs, [mcf[f"exit_{n}_{eQ}_H1200"]["e1"] for n in xs], mk, ms=4, mfc="white", color=K)
    Hn = [RM.rates(int(n), 550, 53, 25)["Hnear"] for n in Ns]
    a.plot(Ns, Hn, color=G, lw=1.0, ls=":")
    a.plot([], [], color=G, lw=1.0, ls=":", marker="s", ms=4, mfc="white", label="nearest-satellite rule")
    xs = [n for n in sizes if f"near_{n}_H1200" in mcf]
    a.plot(xs, [mcf[f"near_{n}_H1200"]["e1"] for n in xs], "s", ms=4, mfc="white", color=G)
    pts = {("o", K): {792: "fig2/w792/mst", 1584: "conf/los/mst", 3168: "fig2/w3168/mst"},
           ("^", K): {792: "fig2/w792/mst40", 1584: "conf/los/mst40", 3168: "fig2/w3168/mst40"},
           ("s", G): {792: "fig2/w792/geof", 1584: "conf/los/geof10", 3168: "fig2/w3168/geof"}}
    errs_exit, errs_near = [], []
    for (mk, col), d in pts.items():
        for n, rel in d.items():
            R = runs(os.path.join(NOV, rel))
            if len(R) < 10:
                continue
            m, h, _ = mean_of(os.path.join(NOV, rel), "e1")
            a.errorbar([n], [m], yerr=[h], fmt=mk, ms=4.5, color=col, capsize=1.5, lw=0.8)
            if mk == "s":
                errs_near.append(100 * abs(m / mcf[f"near_{n}_H1200"]["e1"] - 1))
            else:
                eQ = 25 if mk == "o" else 40
                errs_exit.append(100 * abs(m / mcf[f"exit_{n}_{eQ}_H1200"]["e1"] - 1))
    if len(runs(os.path.join(NOV, "fig2/w3168/mst_k32"))) == 10:
        m32, h32, _ = mean_of(os.path.join(NOV, "fig2/w3168/mst_k32"), "e1")
        a.errorbar([3168 * 1.07], [m32], yerr=[h32], fmt="D", ms=4.5, color=K, mfc="#b8c0c8", capsize=1.5, lw=0.8)
    a.plot([], [], "o", mfc="white", color=K, ms=4, label="Monte Carlo (open)")
    a.plot([], [], "o", color=K, ms=4.5, label="ns-3 (filled)")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xticks(sizes); a.set_xticklabels([str(v) for v in sizes])
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    a.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.0])); a.yaxis.set_minor_locator(NullLocator())
    a.set_yticklabels(["0.2", "0.3", "0.4", "0.6", "1", "1.5", "2"]); a.set_ylim(0.18, 4.2)
    a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r"satellites in the shell, $N_s$"); a.set_ylabel("handovers per terminal-minute")
    a.set_title(r"(a) shell size, 30$^\circ$N", fontsize=8)
    a.legend(loc="upper left", frameon=False, ncol=2, handlelength=2.2, fontsize=7, columnspacing=0.8)
    b = ax[1]
    lat = json.load(open(os.path.join(HERE, "mc_lat.json")))
    L = sorted(float(k) for k in lat)
    lf = np.linspace(0, 50, 26)
    b.plot(lf, [RM.rates(1584, 550, 53, 25, x)["Hp"] for x in lf], color=K, lw=1.0, ls="-", label="exit rule, analysis")
    b.plot(lf, [RM.rates(1584, 550, 53, 25, x)["Hnear"] for x in lf], color=G, lw=1.0, ls=":", label="nearest, analysis")
    b.plot(L, [lat[str(x)]["mc_exit"] for x in L], "o", ms=4, mfc="white", color=K, label="exit rule, Monte Carlo")
    b.plot(L, [lat[str(x)]["mc_near"] for x in L], "s", ms=4, mfc="white", color=G, label="nearest, Monte Carlo")
    b.set_yscale("log"); b.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 0.7, 1.0, 1.5])); b.yaxis.set_minor_locator(NullLocator())
    b.set_yticklabels(["0.2", "0.3", "0.5", "0.7", "1", "1.5"]); b.set_ylim(0.18, 4.2)
    b.set_xlabel(r"terminal latitude $\varphi$ (deg)"); b.set_title("(b) latitude, 1584/72/17 shell", fontsize=8)
    b.legend(loc="upper left", frameon=False, ncol=2, handlelength=2.2, fontsize=7, columnspacing=0.8)
    for x in ax:
        x.grid(True, alpha=0.25, lw=0.4)
    fig.tight_layout(pad=0.3, w_pad=1.2)
    fig.savefig(os.path.join(FIG, "fig_rate.pdf"))
    plt.close(fig)
    if errs_exit:
        put("NvFigTwoExitErrMax", fmt(max(errs_exit), 0))
        m8, _, _ = mean_of(os.path.join(NOV, "fig2/w3168/mst"), "e1")
        mc = mcf["exit_3168_25_H1200"]["e1"]
        put("NvFigTwoDenseMstErr", fmt(100 * (m8 / mc - 1), 0))
        rest = sorted(errs_exit)[:-1]
        put("NvFigTwoExitErrRest", fmt(max(rest), 0))
        if len(runs(os.path.join(NOV, "fig2/w3168/mst_k32"))) == 10:
            m32, _, _ = mean_of(os.path.join(NOV, "fig2/w3168/mst_k32"), "e1")
            put("NvFigTwoDenseKErr", fmt(100 * (m32 / mc - 1), 0))
            pass
    if errs_near:
        put("NvFigTwoNearErrMax", fmt(max(errs_near), 0))
    ex = [100 * (lat[str(x)]["an_exit"] / lat[str(x)]["mc_exit"] - 1) for x in L]
    ne = [100 * (lat[str(x)]["an_near"] / lat[str(x)]["mc_near"] - 1) for x in L]
    put("NvLatExitErrMin", fmt(min(ex), 0)); put("NvLatExitErrMax", fmt(max(ex), 0))
    put("NvLatNearErrMin", fmt(min(ne), 0)); put("NvLatNearErrMax", fmt(max(ne), 0))
    rat = [lat[str(x)]["mc_near"] / lat[str(x)]["mc_exit"] for x in L]
    put("NvLatRatioMin", fmt(min(rat), 1)); put("NvLatRatioMax", fmt(max(rat), 1))


def eta(h_km, lat=30.0, inc=53.0):
    RE, MU, OM = 6378.137, 398600.4418, 7.2921159e-5
    rs = RE + h_km
    vg = math.sqrt(MU / rs) * RE / rs
    la = math.radians(lat)
    s = min(1.0, math.cos(math.radians(inc)) / math.cos(la))
    return vg / math.hypot(vg * math.sqrt(1 - s * s), vg * s - OM * RE * math.cos(la))


def fo_rows(ch, pol="mst", h_km=470.0):
    """planned exits of the real-shell runs: observed exit error against -eta delta_A / v (Lemma 2),
    eta and v at the shell's 470 km"""
    et = eta(h_km)
    v = math.sqrt(398600.4418e9 / ((6378.137 + h_km) * 1e3))
    out = []
    base = os.path.join(R2, "realtle", ch, pol)
    for s in sorted(os.listdir(base)):
        d = os.path.join(base, s)
        off = {int(o["sat"]): o for o in csv.DictReader(open(os.path.join(d, "eph_offsets.csv")))}
        T = float(next(csv.DictReader(open(os.path.join(d, "summary.csv"))))["sim_time_s"])
        for e in csv.DictReader(open(os.path.join(d, "events.csv"))):
            if e["kind"] != "ho" or e["trigger"] != "planned" or e["true_serv_exit_s"] in ("", None):
                continue
            o = off.get(int(e["from_sat"]))
            if o is None:
                continue
            t = float(e["t_s"])
            da = float(o["d0_a_m"]) + (float(o["dT_a_m"]) - float(o["d0_a_m"])) * t / T
            out.append((float(e["pred_serv_tte_s"]) - float(e["true_serv_exit_s"]), -et * da / v))
    return out


def figures(mo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use(os.path.abspath(os.path.join(PAPER, "../templates/ieee.mplstyle")))
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "pdf.fonttype": 42, "ps.fonttype": 42, "legend.fontsize": 7,
                         "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})
    BLUE, ORANGE, GREEN, GRAY, INK = "#2a78d6", "#eb6834", "#1baf7a", "#7b8794", "#1f2933"
    # ---- Fig. 3 (margin)
    rows = fo_rows("los") + fo_rows("suburban")
    x = np.array([r[1] for r in rows])
    y = np.array([r[0] for r in rows])
    b, a = np.polyfit(x, y, 1)
    rr = float(np.corrcoef(x, y)[0, 1])
    put("PropSlopeV", fmt(b, 2))
    put("PropCorrV", fmt(rr, 3))
    put("PropNV", str(len(x)))
    put("NvEtaReal", fmt(eta(470.0), 3))
    fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.2))
    ax[0].plot(x, y, "o", ms=1.6, color=BLUE, alpha=0.35, mew=0)
    lim = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax[0].plot(lim, lim, color=INK, lw=0.9, ls="--", label="equality")
    ax[0].set_xlabel(r"first-order exit error $-\eta\delta_A/v$ (s)")
    ax[0].set_ylabel(r"simulated exit error $\hat{\tau}-\tau$ (s)")
    ax[0].text(0.04, 0.95, f"slope {b:.2f}, r = {rr:.3f}\nn = {len(x)}", transform=ax[0].transAxes, va="top",
               fontsize=7.5)
    ax[0].set_title("(a) MST, both channels", fontsize=8)
    R = mo["los"]
    P, zm = R["P"], R["zm"]
    zz = sorted(zm)
    zf = np.linspace(0, 6.3, 64)
    exits_share = 1.0  # late shares below are per handover; scale the model by the exit share of the runs
    Rl = MM.load("los", range(1, 11))
    exits_share = sum(r["exits"] for r in Rl) / sum(r["n_ho"] for r in Rl)
    ax[1].errorbar(zz, [zm[z][3] for z in zz], yerr=[tci([100.0 * r["late"] / max(1, r["n_ho"]) for r in Rl
                                                         if r["z"] == z])[1] for z in zz],
                   fmt="o", ms=3.4, color=ORANGE, mfc="white", capsize=1.5, lw=0.9, label="simulation, seeds 1-10")
    ax[1].plot(zf, [exits_share * MM.late_model(P, z) for z in zf], color=BLUE, lw=1.1,
               label="analysis, element-set offsets")
    ax[1].plot(zf, [exits_share * 100 * stats.norm.sf(z + (MM.TX - MM.TS / 2) / P["sbar"]) for z in zf],
               color=GRAY, lw=1.0, ls="--", label=r"analysis, Gaussian $Q(z)$")
    if "oos" in R:
        mst = mean_of(os.path.join(NOV, "oos", "los", "mst"), "late")
        tte = mean_of(os.path.join(NOV, "oos", "los", "tte"), "late")
        ax[1].errorbar([0, 2], [mst[0], tte[0]], yerr=[mst[1], tte[1]], fmt="s", ms=3.4, color=GREEN, capsize=1.5,
                       lw=0.9, label="held out: seeds 11-20, +12 h")
    ax[1].set_xlabel(r"margin multiplier $z$")
    ax[1].set_ylabel("late handovers (%)")
    ax[1].set_title("(b) late decisions, line of sight", fontsize=8)
    ax[2].errorbar(zz, [zm[z][0] for z in zz], yerr=[zm[z][1] for z in zz], fmt="o", ms=3.4, color=ORANGE,
                   mfc="white", capsize=1.5, lw=0.9)
    ax[2].plot(zf, [MM.J(P, z) for z in zf], color=BLUE, lw=1.1)
    zs = MM.zstar_gauss(P)
    ze = MM.zstar_numeric(P)
    top = ax[2].get_ylim()[1]
    ax[2].axvline(zs, color=GRAY, ls=":", lw=0.9)
    ax[2].text(zs + 0.12, top * 0.98, r"$z^\star$", fontsize=8, va="top", color=GRAY)
    ax[2].axvline(ze, color=BLUE, ls="--", lw=0.8)
    ax[2].text(ze - 0.12, top * 0.98, r"$\arg\min\,J$", fontsize=8, va="top", ha="right", color=BLUE)
    ax[2].set_xlabel(r"margin multiplier $z$")
    ax[2].set_ylabel("outage (%)")
    ax[2].set_title("(c) outage, line of sight", fontsize=8)
    for a_ in ax:
        a_.grid(True, alpha=0.3)
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.6, 1.0))
    fig.tight_layout(w_pad=1.0, rect=(0, 0, 1, 0.88))
    fig.savefig(os.path.join(FIG, "fig_margin.pdf"))
    plt.close(fig)


def fig_rate():
    """Fig. 2: analysis lines as before; Monte Carlo and ns-3 both read with estimator E1 at 600 s
    (exploratory runs), plus the confirmatory 1200 s points at the reference shell."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(SIM, "analysis_rate"))
    import rate_model as RM
    plt.style.use(os.path.abspath(os.path.join(PAPER, "../templates/ieee.mplstyle")))
    plt.rcParams.update({"font.family": "serif", "pdf.fonttype": 42, "ps.fonttype": 42, "legend.fontsize": 6,
                         "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})
    BLUE, ORANGE, GREEN, INK = "#2a78d6", "#eb6834", "#1baf7a", "#1f2933"
    mcf = json.load(open(os.path.join(HERE, "mc_fig.json")))
    fig, ax = plt.subplots(figsize=(3.45, 2.6))
    Ns = np.geomspace(396, 6336, 13)
    sizes = [396, 792, 1584, 3168, 6336]
    for eQ, col, lab in ((25, BLUE, r"exit rule, $e_Q=25^\circ$"), (40, ORANGE, r"exit rule, $e_Q=40^\circ$")):
        Hc = [60 / RM.poisson_chain(int(n), 550, 53, eQ, K=4000) for n in Ns]
        Ha = [60 / RM.asym(int(n), 550, 53, eQ) if RM.rates(int(n), 550, 53, eQ)["Nvis"] >= 2.0 else np.nan for n in Ns]
        ax.plot(Ns, Hc, color=col, lw=1.1, label=lab)
        ax.plot(Ns, Ha, color=col, lw=0.9, ls="--")
        xs = [n for n in sizes if f"exit_{n}_{eQ}_H600" in mcf]
        ax.plot(xs, [mcf[f"exit_{n}_{eQ}_H600"]["e1"] for n in xs], "o", ms=3.5, mfc="white", color=col)
    Hn = [RM.rates(int(n), 550, 53, 25)["Hnear"] for n in Ns]
    ax.plot(Ns, Hn, color=GREEN, lw=1.1, ls=":", label="nearest-satellite rule")
    ax.plot(sizes, [mcf[f"near_{n}_H600"]["e1"] for n in sizes], "s", ms=3.2, mfc="white", color=GREEN)
    # ns-3, 600 s exploratory, read with E1
    pts = {(1584, 25): os.path.join(R2, "main", "los", "mst"),
           (792, 40): os.path.join(R2, "density", "w792", "los", "mst40"),
           (1584, 40): os.path.join(R2, "main", "los", "mst40"),
           (3168, 40): os.path.join(R2, "density", "w3168", "los", "mst40")}
    for (n, eQ), p in pts.items():
        m, h, _ = mean_of(p, "e1")
        ax.errorbar([n], [m], yerr=[h], fmt="o", ms=4.3, color=BLUE if eQ == 25 else ORANGE, capsize=1.5, lw=0.8)
    for n, p in ((792, os.path.join(R2, "density", "w792", "los", "geof")), (1584, os.path.join(R2, "main", "los", "geof")),
                 (3168, os.path.join(R2, "density", "w3168", "los", "geof"))):
        m, h, _ = mean_of(p, "e1")
        ax.errorbar([n], [m], yerr=[h], fmt="s", ms=4.3, color=GREEN, capsize=1.5, lw=0.8)
    ax.plot([], [], color=INK, lw=0.9, ls="--", label=r"dense limit, $\bar{N}\geq 2$")
    ax.plot([], [], "o", mfc="white", color=INK, ms=3.5, label="Walker Monte Carlo")
    ax.plot([], [], "o", color=INK, ms=4.3, label="ns-3 (MST, MST-40, geodesic)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(sizes)
    ax.set_xticklabels([str(s) for s in sizes])
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    ax.yaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.4, 0.6, 1.0, 1.5, 2.0]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels(["0.2", "0.3", "0.4", "0.6", "1", "1.5", "2"])
    ax.set_ylim(0.18, 2.1)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("satellites in the shell, $N_s$")
    ax.set_ylabel("handovers per terminal-minute")
    ax.legend(loc="upper left", frameon=False, ncol=1)
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(FIG, "fig_rate.pdf"))
    plt.close(fig)
    # ns-3 vs Monte Carlo at 600 s, same estimator
    errs = []
    for (n, eQ), p in pts.items():
        errs.append(100 * abs(mean_of(p, "e1")[0] / mcf[f"exit_{n}_{eQ}_H600"]["e1"] - 1))
    put("NvFigExitErrMax", fmt(max(errs), 0))
    gerr = []
    for n, p in ((792, os.path.join(R2, "density", "w792", "los", "geof")), (1584, os.path.join(R2, "main", "los", "geof")),
                 (3168, os.path.join(R2, "density", "w3168", "los", "geof"))):
        gerr.append(100 * abs(mean_of(p, "e1")[0] / mcf[f"near_{n}_H600"]["e1"] - 1))
    put("NvFigGeoErrMax", fmt(max(gerr), 0))


if __name__ == "__main__":
    res = rate_block()
    if os.path.exists(os.path.join(NOV, "rival_choice.json")):
        table_tune()
        if all(runs(C(e, p)) for e in ("los", "suburban") for p in ("ut", "tted2f", "geof2")):
            table_conf(res)
    attr_block()
    perpass_block()
    mo = margin_block()
    window_block()
    ms = [float(e["margin_s"]) for d in sorted(os.listdir(os.path.join(NOV, "oos", "los", "tte")))
          for e in csv.DictReader(open(os.path.join(NOV, "oos", "los", "tte", d, "events.csv")))
          if e["kind"] == "ho" and e["trigger"] == "planned"]
    put("NvOosLosMargin", fmt(st.median(ms), 2))
    table_margin(mo)
    figures(mo)
    if os.path.exists(os.path.join(HERE, "mc_lat.json")):
        fig_rate2()
    elif os.path.exists(os.path.join(HERE, "mc_fig.json")):
        fig_rate()
    with open(os.path.join(PAPER, "macros_nov.tex"), "w") as f:
        f.write("% written by sim/analysis_nov/analyze_nov.py; do not edit\n")
        for k in sorted(M):
            f.write(f"\\newcommand{{\\{k}}}{{{M[k]}}}\n")
    json.dump(REPORT, open(os.path.join(HERE, "report.json"), "w"), indent=1, default=float)
    print(len(M), "macros")
