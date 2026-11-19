#!/usr/bin/env python3
"""Turn the campaign into the numbers, tables and figures the
paper uses.

Every number quoted in main_v6.tex comes from results_macros.tex, which this
script writes from the run outputs under runs_r2/. Nothing is typed by hand.

Statistics. Every comparison is seed-paired (the same orbits, terminals and
random draws). A difference is reported as the mean paired difference with its
95 % Student-t interval, and tested with the exact two-sided Wilcoxon
signed-rank test. Each table's comparisons form one pre-declared family
(FAMILIES below), Holm-adjusted at alpha = 0.05; the text calls a difference
detected only if its Holm-adjusted p is below 0.05.

    python3 analyze.py
"""
import csv
import glob
import itertools
import json
import math
import os
import statistics as st
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, os.environ.get("RUNS_DIR", "runs_r2"))
PAPER = os.path.abspath(os.path.join(HERE, ".."))
FIG = os.path.join(PAPER, os.environ.get("FIG_DIR", "figures_v6"))
os.makedirs(FIG, exist_ok=True)
ALPHA = 0.05

# Student-t 0.975 quantiles
T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262,
        10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 17: 2.110,
        18: 2.101, 19: 2.093, 20: 2.086, 25: 2.060, 29: 2.045, 30: 2.042}


def tq(df):
    if df in T975:
        return T975[df]
    return 1.96 if df > 30 else T975[min(T975, key=lambda k: abs(k - df))]


def ci(vals):
    vals = [v for v in vals if v is not None and not math.isnan(v)]
    n = len(vals)
    if n == 0:
        return float("nan"), float("nan"), 0
    m = st.mean(vals)
    if n == 1:
        return m, 0.0, 1
    return m, tq(n - 1) * st.stdev(vals) / math.sqrt(n), n


# ----------------------------------------------------------------- loading
_SUM = None


def summaries():
    """{tag: {run: row}} for every run under RUNS."""
    global _SUM
    if _SUM is None:
        _SUM = defaultdict(dict)
        for f in sorted(glob.glob(os.path.join(RUNS, "**", "summary.csv"), recursive=True)):
            tag = os.path.relpath(os.path.dirname(os.path.dirname(f)), RUNS)
            with open(f) as fh:
                row = next(csv.DictReader(fh))
            row["_dir"] = os.path.dirname(f)
            _SUM[tag][int(row["run"])] = row
    return _SUM


def num(row, k):
    try:
        return float(row[k])
    except (KeyError, ValueError, TypeError):
        return float("nan")


# metric -> function(row) giving the value reported
def _per_hour(col):
    return lambda r: num(r, col) / (num(r, "num_ues") * num(r, "sim_time_s") / 3600.0)


METRIC = {
    "ho": lambda r: num(r, "ho_per_ue_min"),
    "ppone": lambda r: 100.0 * num(r, "pp1") / max(1.0, num(r, "ho_success")),
    "ppfive": lambda r: 100.0 * num(r, "pp5") / max(1.0, num(r, "ho_success")),
    "rlf": _per_hour("rlf"),
    "hof": _per_hour("hof"),
    "out": lambda r: 100.0 * num(r, "outage_frac"),
    "snr": lambda r: num(r, "mean_sinr_db"),
    "se": lambda r: num(r, "mean_se"),
    "tos": lambda r: num(r, "median_tos_s"),
    "esc": lambda r: 100.0 * num(r, "escape_share"),
    "rrc": lambda r: num(r, "rrc_per_ue_min"),
    "rrccand": lambda r: num(r, "rrc_cand_per_ue_min"),
    "rrcwin": lambda r: num(r, "rrc_win_per_ue_min"),
    "evals": lambda r: num(r, "tte_evals_per_decision"),
    "wall": lambda r: num(r, "wall_s"),
}
DIG = {"ho": 3, "ppone": 2, "ppfive": 1, "rlf": 2, "hof": 2, "out": 3, "snr": 2, "se": 2, "tos": 0, "esc": 0,
       "rrc": 1, "rrccand": 1, "rrcwin": 1, "evals": 2, "wall": 0}


def vals(tag, metric):
    return {k: METRIC[metric](r) for k, r in summaries().get(tag, {}).items()}


def stat(tag, metric):
    return ci(list(vals(tag, metric).values()))


# ----------------------------------------------------------------- statistics
def wilcoxon_p(d):
    """Exact two-sided Wilcoxon signed-rank p for paired differences d (mid-ranks)."""
    d = [round(x, 9) for x in d]
    d = [x for x in d if x != 0.0]
    n = len(d)
    if n == 0:
        return 1.0
    order = sorted(range(n), key=lambda i: abs(d[i]))
    r = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(d[order[j + 1]]) == abs(d[order[i]]):
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    wp = sum(r[i] for i in range(n) if d[i] > 0)
    tot = n * (n + 1) / 2
    obs = min(wp, tot - wp)
    if n > 20:
        mu = tot / 2
        sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
        return min(1.0, math.erfc(abs(wp - mu) / sd / math.sqrt(2)))
    # exact distribution by dynamic programming over the (half-)ranks
    scale = 2
    rr = [int(round(x * scale)) for x in r]
    dist = {0: 1}
    for x in rr:
        nd = defaultdict(int)
        for s, c in dist.items():
            nd[s] += c
            nd[s + x] += c
        dist = nd
    tt = sum(rr)
    o = int(round(obs * scale))
    cnt = sum(c for s, c in dist.items() if min(s, tt - s) <= o)
    return min(1.0, cnt / 2 ** n)


def paired(tagA, tagB, metric):
    a, b = vals(tagA, metric), vals(tagB, metric)
    k = sorted(set(a) & set(b))
    d = [a[x] - b[x] for x in k]
    if not d:
        return None
    m, h, n = ci(d)
    mb = st.mean([b[x] for x in k])
    return {"diff": m, "ci": h, "n": n, "p": wilcoxon_p(d), "rel": 100.0 * m / mb if mb else float("nan"),
            "lower": sum(1 for x in d if x < -1e-12), "higher": sum(1 for x in d if x > 1e-12),
            "a": st.mean([a[x] for x in k]), "b": mb}


def holm(ps):
    """Holm step-down adjusted p-values, same order as ps."""
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj = [0.0] * m
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[i]))
        adj[i] = run
    return adj


# ----------------------------------------------------------------- macros
macros = {}


def put(name, val):
    assert name.isalpha(), name
    macros[name] = val


def fmt(x, d):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "--"
    s = f"{x:.{d}f}"
    if s.startswith("-") and float(s) == 0.0:
        s = s[1:]
    return s


def pfmt(p):
    if p < 0.001:
        return r"$p<0.001$"
    return f"$p={p:.3f}$"


POLN = {"tteage": "Tteage", "utage": "Utage", "ttefz": "Ttefz", "a3": "Athree", "a3t": "Atuned", "d2": "Dtwo", "t1": "Tone", "t1w1": "Tonew", "elev": "Elev",
        "a3tf": "Atunedf", "d2f": "Dtwof", "t1f": "Tonef", "elevf": "Elevf", "a3tfh": "Atunedfh",
        "d2fh": "Dtwofh", "t1fh": "Tonefh", "elevfh": "Elevfh", "t1f40": "Toneforty", "elevf40": "Elevforty",
        "d2f4050": "Dtwoforty", "eager": "Eager", "tted2": "Ttedtwo", "tted2f": "Ttedtwof", "mst": "Mst",
        "mst40": "Mstforty", "mst40s": "Mstfortys", "mst_prep": "Mstprep", "mst_prepr": "Mstprepr",
        "mst_ute": "Mstute", "tte": "Tte", "tte40": "Tteforty", "ut": "Ut", "ut40": "Utforty",
        "ut40_norank": "Ablrank", "ut40_noprep": "Ablprep", "ut40_noreserve": "Ablreserve",
        "ut40_althold": "Ablhold", "ut40_qesc": "Ablesc", "ut40_oldesc": "Oldesc",
        "geof": "Geof", "wangf": "Wangf", "zabf": "Zabf"}
LAB = {"a3": "A3 (3 dB, 160 ms)", "a3t": "A3-tuned (6 dB, 640 ms)", "d2": "D2", "t1": "T1",
       "t1w1": "T1-1s", "elev": "Elevation", "a3tf": "A3-tuned + net", "d2f": "D2 + net", "t1f": "T1 + net",
       "elevf": "Elevation + net", "a3tfh": "A3-tuned + net + hold", "d2fh": "D2 + net + hold",
       "t1fh": "T1 + net + hold", "elevfh": "Elevation + net + hold", "t1f40": r"T1-40 + net",
       "elevf40": r"Elevation-40 + net", "d2f4050": r"D2-40/50 + net", "eager": "TTE-eager",
       "tted2": "TTE-D2", "tted2f": "TTE-D2 + net", "mst": "MST", "mst40": "MST-40",
       "mst40s": r"MST-40, $T_f$=640\,ms", "mst_prep": r"MST + prep.", "mst_prepr": r"MST + prep.\ + $r$",
       "mst_ute": r"MST + $U$ rank", "tte": "MST + margin", "tte40": "MST-40 + margin", "ut": "UT-CHO",
       "ut40": "UT-CHO-40", "ut40_norank": r"UT-CHO-40, ranked by $\hat{\tau}$",
       "ut40_noprep": "UT-CHO-40, strongest prepared", "ut40_noreserve": "UT-CHO-40, no reserved slots",
       "ut40_althold": "UT-CHO-40 + 5\\,s hold-off", "ut40_qesc": "UT-CHO-40, escape to strongest", "ut40_oldesc": "UT-CHO-40, earlier escape + hold-off",
       "geof": "Geodesic CHO + net", "wangf": "Entropy-weighted + net", "zabf": "Residence priority + net"}
CHN = {"los": "Los", "suburban": "Sub", "urban": "Urb", "denseurban": "Dense"}
MNAME = {"ho": "Ho", "ppone": "Ppone", "ppfive": "Ppfive", "rlf": "Rlf", "hof": "Hof", "out": "Out",
         "snr": "Snr", "se": "Se", "tos": "Tos", "esc": "Esc", "rrc": "Rrc", "rrccand": "Rrccand",
         "rrcwin": "Rrcwin", "evals": "Evals", "wall": "Wall"}


def tag_of(exp, ch, pol):
    return f"{exp}/{ch}/{pol}"


def base_macros():
    """Mean and CI of every metric for every policy of main and urban."""
    s = summaries()
    for tag in sorted(s):
        parts = tag.split("/")
        if len(parts) != 3 or parts[0] not in ("main", "urban") or parts[2] not in POLN:
            continue
        exp, ch, pol = parts
        pre = ("Main" if exp == "main" else "Urb") + CHN[ch] + POLN[pol]
        for mk in METRIC:
            m, h, n = stat(tag, mk)
            put(pre + MNAME[mk], fmt(m, DIG[mk]))
            put(pre + MNAME[mk] + "Ci", fmt(h, DIG[mk]))
        put(pre + "N", str(len(s[tag])))


# ----------------------------------------------------------------- families
# (macro name, tag A, tag B, metric). Each family is Holm-adjusted on its own.
def fam_main():
    F = defaultdict(list)
    L = lambda p: tag_of("main", "los", p)
    S = lambda p: tag_of("main", "suburban", p)
    # LOS handover rate: the exit-triggered rules against MST and the standardized conditions
    for a, b in [("ut", "mst"), ("ut", "t1"), ("ut", "t1w1"), ("ut", "a3t"), ("ut", "d2"), ("ut", "elev"),
                 ("mst", "t1"), ("ut40", "mst40"), ("ut40", "a3t"), ("mst40", "a3t")]:
        F["LosHo"].append((f"Los{POLN[a]}Vs{POLN[b]}Ho", L(a), L(b), "ho"))
    for a, b in [("ut40", "a3t"), ("ut40", "mst40"), ("ut", "mst"), ("mst40", "a3t")]:
        F["LosSnr"].append((f"Los{POLN[a]}Vs{POLN[b]}Snr", L(a), L(b), "snr"))
    # suburban: against the standardized conditions with the same net, and MST-40
    for b in ["a3tf", "d2f", "t1f", "elevf", "mst40"]:
        for mk in ["out", "ho", "ppone"]:
            F["Sub" + MNAME[mk]].append((f"Sub{POLN['ut40']}Vs{POLN[b]}{MNAME[mk]}", S("ut40"), S(b), mk))
    for b in ["mst"]:
        for mk in ["out", "ho"]:
            F["Sub" + MNAME[mk]].append((f"Sub{POLN['ut']}Vs{POLN[b]}{MNAME[mk]}", S("ut"), S(b), mk))
    # mechanisms (Table ctrl): one factor each, LOS and suburban
    for a, b, mk in [("mst_prep", "mst", "ho"), ("mst_prepr", "mst", "ho"), ("mst_ute", "mst", "ho"),
                     ("ut", "mst_prepr", "ho"), ("t1w1", "t1", "ho"), ("tted2", "d2", "ho"), ("mst", "tted2", "ho")]:
        F["CtrlLos"].append((f"Ctrl{POLN[a]}Vs{POLN[b]}{MNAME[mk]}Los", L(a), L(b), mk))
    for a in ["ut40_norank", "ut40_noprep", "ut40_noreserve", "ut40_qesc", "ut40_althold", "ut40_oldesc"]:
        for mk in ["out", "ho", "ppone", "rlf"]:
            F["AblSub" + MNAME[mk]].append((f"Abl{POLN[a]}{MNAME[mk]}Sub", S("ut40"), S(a), mk))
    return F


def fam_urban():
    F = defaultdict(list)
    base = ["a3tf", "d2f", "t1f", "elevf", "tted2f", "a3tfh", "d2fh", "t1fh", "elevfh", "t1f40", "elevf40",
            "d2f4050", "mst40"]
    for ch in ["urban", "denseurban"]:
        U = lambda p: tag_of("urban", ch, p)
        for b in base:
            for mk in ["out", "ho", "hof", "ppone"]:
                F[CHN[ch] + MNAME[mk]].append((f"{CHN[ch]}{POLN['ut40']}Vs{POLN[b]}{MNAME[mk]}", U("ut40"), U(b), mk))
        # MST-40 (the contour and the net) against the standardized conditions
        for b in base[:-1]:
            F[CHN[ch] + "MstOut"].append((f"{CHN[ch]}{POLN['mst40']}Vs{POLN[b]}Out", U("mst40"), U(b), "out"))
        for a in ["ut40_norank", "ut40_noprep", "ut40_noreserve", "ut40_qesc", "ut40_althold", "ut40_oldesc"]:
            for mk in ["out", "ho", "rlf", "hof", "ppone"]:
                F["Abl" + CHN[ch] + MNAME[mk]].append((f"Abl{POLN[a]}{MNAME[mk]}{CHN[ch]}", U("ut40"), U(a), mk))
        for a, b in [("ut", "mst"), ("mst_ute", "mst"), ("mst_prepr", "mst")]:
            for mk in ["out", "ho"]:
                F["Edge" + CHN[ch]].append((f"{CHN[ch]}{POLN[a]}Vs{POLN[b]}{MNAME[mk]}", U(a), U(b), mk))
    return F


RESULTS = {}  # macro name -> paired dict with holm


def run_families(F):
    for fam, items in F.items():
        res = []
        for nm, a, b, mk in items:
            r = paired(a, b, mk)
            if r is not None:
                res.append((nm, mk, r))
        if not res:
            continue
        adj = holm([r["p"] for _, _, r in res])
        put("Fam" + fam + "Size", str(len(res)))
        nsig = nlow = 0
        for (nm, mk, r), pa in zip(res, adj):
            r["holm"] = pa
            r["fam"] = fam
            r["m"] = len(res)
            RESULTS[nm] = r
            d = DIG[mk]
            put("D" + nm, fmt(r["diff"], d))
            put("D" + nm + "Ci", fmt(r["ci"], d))
            put("D" + nm + "Abs", fmt(abs(r["diff"]), d))
            put("D" + nm + "Rel", fmt(abs(r["rel"]), 1))
            put("P" + nm, pfmt(r["p"]))
            put("H" + nm, pfmt(pa))
            put("S" + nm, "detected" if pa < ALPHA else "not detected")
            put("Lw" + nm, str(r["lower"]))
            put("Hi" + nm, str(r["higher"]))
            nsig += pa < ALPHA
            nlow += pa < ALPHA and r["diff"] < 0
        put("Fam" + fam + "Sig", str(nsig))
        put("Fam" + fam + "SigLower", str(nlow))


# ----------------------------------------------------------------- tables
def cell(tag, mk, d=None, withci=False):
    m, h, n = stat(tag, mk)
    if n == 0:
        return "--"
    d = DIG[mk] if d is None else d
    return fmt(m, d) + (f"\\,$\\pm$\\,{fmt(h, d)}" if withci else "")


MAIN_ROWS = [("a3", 0), ("a3t", 0), ("d2", 0), ("t1", 0), ("elev", 0), None,
             ("a3tf", 0), ("d2f", 0), ("t1f", 0), ("elevf", 0), None,
             ("mst", 0), ("mst40", 0), None, ("ut", 1), ("ut40", 1)]


def table_main():
    L = [r"\begin{table*}[!t]",
         r"\caption{Reference comparison, " + macros.get("NSeedMain", "ten") + r" seeds (means; handover rate and "
         r"outage with their 95\,\% CI). Rates per terminal-minute, failures per terminal-hour. ``+ net'': "
         r"the measurement fallback that every exit-triggered rule also carries.}",
         r"\label{tab:main}", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.6pt}",
         r"\begin{tabular}{@{}l" + "r" * 6 + r"@{\hspace{7pt}}" + "r" * 7 + "@{}}", r"\toprule",
         r" & \multicolumn{6}{c}{Line of sight} & \multicolumn{7}{c}{Suburban (TR 38.811 LOS process)}\\",
         r"\cmidrule(lr){2-7}\cmidrule(l){8-14}",
         r"Policy & HO & PP$_1$ & Out.\,\% & SNR & ToS & RRC & HO & PP$_1$ & RLF & HOF & Out.\,\% & SNR & Esc.\,\% \\",
         r"\midrule"]
    for row in MAIN_ROWS:
        if row is None:
            L.append(r"\midrule")
            continue
        key, bold = row
        T = lambda ch: tag_of("main", ch, key)
        c = [cell(T("los"), "ho", 3, True), cell(T("los"), "ppone", 2), cell(T("los"), "out", 3),
             cell(T("los"), "snr", 1), cell(T("los"), "tos", 0), cell(T("los"), "rrc", 1),
             cell(T("suburban"), "ho", 3, True), cell(T("suburban"), "ppone", 2), cell(T("suburban"), "rlf", 2),
             cell(T("suburban"), "hof", 2), cell(T("suburban"), "out", 3, True), cell(T("suburban"), "snr", 1),
             cell(T("suburban"), "esc", 0)]
        name = LAB[key]
        if bold:
            name = r"\textbf{" + name + "}"
        L.append(name + " & " + " & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize HO: completed handovers; PP$_1$: ping-pong, \% of handovers returning to the "
          r"previous cell after a stay below 1\,s (TR 36.839); Out.: outage; SNR in dB; ToS: median time of stay "
          r"(s); RRC: conditional reconfigurations per terminal-minute; Esc.: share of handovers taken by the "
          r"safety net. D2 and T1 carry an A4 leg on the same candidate (threshold $Q$); A3 needs none.}",
          r"\end{table*}"]
    open(os.path.join(PAPER, "table_main.tex"), "w").write("\n".join(L) + "\n")


CTRL_ROWS = ["t1w1", "a3tfh", "d2fh", "t1fh", "elevfh", "tted2", "tted2f", "eager", "mst40s", "mst_prep",
             "mst_prepr", "mst_ute", "ut40_norank", "ut40_noprep", "ut40_noreserve", "ut40_qesc", "ut40_althold", "ut40_oldesc"]


def table_ctrl():
    L = [r"\begin{table}[!t]",
         r"\caption{Single-factor controls (means). Each row changes one element of the policy it names.}",
         r"\label{tab:ctrl}", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.2pt}",
         r"\begin{tabular}{@{}lrrrrrr@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{Line of sight} & \multicolumn{4}{c}{Suburban} \\",
         r"\cmidrule(lr){2-3}\cmidrule(l){4-7}",
         r"Policy & HO & SNR & HO & RLF & Out.\,\% & PP$_1$ \\", r"\midrule"]
    for key in CTRL_ROWS:
        T = lambda ch: tag_of("main", ch, key)
        if not summaries().get(T("los")):
            continue
        L.append(f"{LAB[key]} & {cell(T('los'), 'ho', 3)} & {cell(T('los'), 'snr', 2)} & "
                 f"{cell(T('suburban'), 'ho', 3)} & {cell(T('suburban'), 'rlf', 2)} & "
                 f"{cell(T('suburban'), 'out', 3)} & {cell(T('suburban'), 'ppone', 2)} \\\\")
    L.append(r"\midrule")
    for key in ["mst", "mst40", "ut", "ut40"]:
        T = lambda ch: tag_of("main", ch, key)
        lab = LAB[key] if key.startswith("mst") else r"\textbf{" + LAB[key] + "}"
        L.append(f"{lab} & {cell(T('los'), 'ho', 3)} & {cell(T('los'), 'snr', 2)} & "
                 f"{cell(T('suburban'), 'ho', 3)} & {cell(T('suburban'), 'rlf', 2)} & "
                 f"{cell(T('suburban'), 'out', 3)} & {cell(T('suburban'), 'ppone', 2)} \\\\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize HO per terminal-minute, RLF per terminal-hour, SNR in dB, PP$_1$ in \%. "
          r"TTE-D2: D2 execution with the admitted candidates ranked by predicted time to exit, "
          r"after~\cite{noreen2026location}. $r$: two prepared slots reserved for the strongest satellites.}",
          r"\end{table}"]
    open(os.path.join(PAPER, "table_ctrl.tex"), "w").write("\n".join(L) + "\n")


URB_ROWS = ["a3tf", "d2f", "t1f", "elevf", "tted2f", None,
            "t1f40", "elevf40", "d2f4050", None, "mst", "mst40", None, "ut", "ut40"]
# hold-off controls and UT-CHO-40 ablations go to the supplement (Table S7)
URB_CTRL_ROWS = ["a3tfh", "d2fh", "t1fh", "elevfh", None, "ut40", None,
                 "ut40_norank", "ut40_noprep", "ut40_noreserve", "ut40_qesc", "ut40_althold", "ut40_oldesc"]


def table_urban(rows=None, fname="table_urban.tex", label="tab:urban"):
    rows = URB_ROWS if rows is None else rows
    L = [r"\begin{table*}[!t]",
         r"\caption{Urban and dense-urban environments (TR 38.811 LOS probability, path loss and shadowing; "
         + macros.get("NSeedUrb", "ten") + r" seeds, means). HO per terminal-minute; ping-pong, outage and "
         r"escape share in \%; failures per terminal-hour. $\Delta$Out.: outage of UT-CHO-40 minus that of the row, "
         r"seed-paired mean in percentage points with its 95\,\% CI; $^{\ast}$: Holm-adjusted $p<0.05$ in the "
         r"family of the column.}",
         r"\label{" + label + "}", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.1pt}",
         r"\begin{tabular}{@{}l" + "r" * 7 + r"@{\hspace{7pt}}" + "r" * 7 + "@{}}", r"\toprule",
         r" & \multicolumn{7}{c}{Urban} & \multicolumn{7}{c}{Dense urban} \\",
         r"\cmidrule(lr){2-8}\cmidrule(l){9-15}",
         r"Policy & HO & PP$_1$ & Out. & $\Delta$Out. & HOF & RLF & Esc. & HO & PP$_1$ & Out. & $\Delta$Out. & HOF & RLF & Esc. \\",
         r"\midrule"]
    for key in rows:
        if key is None:
            L.append(r"\midrule")
            continue
        c = []
        for ch in ["urban", "denseurban"]:
            T = tag_of("urban", ch, key)
            nm = f"{CHN[ch]}{POLN['ut40']}Vs{POLN[key]}Out"
            if key.startswith("ut40_"):
                nm = f"Abl{POLN[key]}Out{CHN[ch]}"
            r = RESULTS.get(nm)
            if r is not None:
                sgn = -r["diff"] if key.startswith("ut40_") else r["diff"]
                dcell = f"{fmt(sgn, 2)}\\,$\\pm$\\,{fmt(r['ci'], 2)}" + (r"$^{\ast}$" if r["holm"] < ALPHA else "")
            else:
                dcell = "--"
            c += [cell(T, "ho", 2), cell(T, "ppone", 1), cell(T, "out", 2), dcell, cell(T, "hof", 2),
                  cell(T, "rlf", 2), cell(T, "esc", 0)]
        lab = LAB[key]
        if key in ("ut", "ut40"):
            lab = r"\textbf{" + lab + "}"
        L.append(lab + " & " + " & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize For the ablation rows $\Delta$Out.\ is the ablation minus UT-CHO-40. "
          r"Aircraft and ships see the sky in line of sight and trains the rural table (Section~\ref{sec:setup}).}",
          r"\end{table*}"]
    open(os.path.join(PAPER, fname), "w").write("\n".join(L) + "\n")


def table_headline():
    """Paired differences behind the abstract and the conclusion, with Holm."""
    rows = [("LosUtVsMstHo", r"UT-CHO vs MST, HO, line of sight"),
            ("LosUtVsToneHo", r"UT-CHO vs T1, HO, line of sight"),
            ("LosMstVsToneHo", r"MST vs T1, HO, line of sight"),
            ("LosUtfortyVsMstfortyHo", r"UT-CHO-40 vs MST-40, HO, line of sight"),
            ("SubUtfortyVsMstfortyOut", r"UT-CHO-40 vs MST-40, outage, suburban"),
            ("SubUtfortyVsDtwofOut", r"UT-CHO-40 vs D2 + net, outage, suburban"),
            ("UrbUtfortyVsDtwofOut", r"UT-CHO-40 vs D2 + net, outage, urban"),
            ("UrbUtfortyVsMstfortyOut", r"UT-CHO-40 vs MST-40, outage, urban"),
            ("UrbUtfortyVsTonefortyOut", r"UT-CHO-40 vs T1-40 + net, outage, urban"),
            ("DenseUtfortyVsDtwofOut", r"UT-CHO-40 vs D2 + net, outage, dense urban"),
            ("DenseUtfortyVsMstfortyOut", r"UT-CHO-40 vs MST-40, outage, dense urban"),
            ("UrbUtfortyVsDtwofHof", r"UT-CHO-40 vs D2 + net, HOF, urban")]
    L = [r"\begin{table}[!t]", r"\caption{Headline comparisons: seed-paired mean difference (first minus second) "
         r"with its 95\,\% CI, the share of seeds in which the first is lower, and the exact Wilcoxon $p$ "
         r"before and after the Holm correction of its family.}", r"\label{tab:headline}", r"\centering",
         r"\scriptsize", r"\setlength{\tabcolsep}{2.0pt}", r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
         r"Comparison & Difference & Lower & $p$ & $p_{\mathrm{Holm}}$ \\", r"\midrule"]
    for nm, lab in rows:
        r = RESULTS.get(nm)
        if r is None:
            continue
        mk = {"Ho": "ho", "Out": "out", "Hof": "hof", "Snr": "snr"}[[s for s in ("Hof", "Ho", "Out", "Snr") if nm.endswith(s)][0]]
        d = DIG[mk]
        L.append(f"{lab} & {fmt(r['diff'], d)}\\,$\\pm$\\,{fmt(r['ci'], d)} & {r['lower']}/{r['n']} & "
                 f"{fmt(r['p'], 3) if r['p'] >= 0.001 else '$<$0.001'} & "
                 f"{fmt(r['holm'], 3) if r['holm'] >= 0.001 else '$<$0.001'} \\\\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize HO per terminal-minute, outage in percentage points, HOF per terminal-hour.}",
          r"\end{table}"]
    open(os.path.join(PAPER, "table_headline.tex"), "w").write("\n".join(L) + "\n")


# ----------------------------------------------------------------- derived headline macros
def derived():
    def m(tag, mk):
        return stat(tag, mk)[0]
    L = lambda p: tag_of("main", "los", p)
    S = lambda p: tag_of("main", "suburban", p)
    cut = lambda a, b: 100.0 * (1 - a / b)
    put("CutUtVsMstLos", fmt(cut(m(L("ut"), "ho"), m(L("mst"), "ho")), 1))
    put("CutUtVsToneLos", fmt(cut(m(L("ut"), "ho"), m(L("t1"), "ho")), 0))
    put("CutMstVsToneLos", fmt(cut(m(L("mst"), "ho"), m(L("t1"), "ho")), 0))
    put("CutUtVsTonewLos", fmt(cut(m(L("ut"), "ho"), m(L("t1w1"), "ho")), 0))
    put("CutMstVsTonewLos", fmt(cut(m(L("mst"), "ho"), m(L("t1w1"), "ho")), 0))
    put("CutUtVsAthreeLos", fmt(cut(m(L("ut"), "ho"), m(L("a3"), "ho")), 0))
    put("CutUtVsAtunedLos", fmt(cut(m(L("ut"), "ho"), m(L("a3t"), "ho")), 0))
    put("CutUtfortyVsAtunedLos", fmt(cut(m(L("ut40"), "ho"), m(L("a3t"), "ho")), 0))
    put("CutUtfortyVsMstfortyLos", fmt(cut(m(L("ut40"), "ho"), m(L("mst40"), "ho")), 1))
    # share of the T1 -> UT gap that MST already closes
    t1, mst, ut = m(L("t1"), "ho"), m(L("mst"), "ho"), m(L("ut"), "ho")
    put("ShareMstOfGapLos", fmt(100.0 * (t1 - mst) / (t1 - ut), 0))
    # T1-1s executes at nearly MST's instant: split the T1 -> MST gap into timing and target
    t1w = m(L("t1w1"), "ho")
    put("ShareTargetLos", fmt(100.0 * (t1w - mst) / (t1 - mst), 0))
    # the most frugal standardized setting (any standardized condition or control)
    std = [k for k in ["a3", "a3t", "d2", "t1", "t1w1", "elev", "a3tf", "d2f", "t1f", "elevf"] if
           summaries().get(L(k))]
    best = min(std, key=lambda k: m(L(k), "ho"))
    put("BestStdNameLos", LAB[best].split(" (")[0])
    put("BestStdHoLos", fmt(m(L(best), "ho"), 3))
    # identity: MST + margin == MST at zero ephemeris error, per seed
    for ch in ["los", "suburban"]:
        a, b = vals(tag_of("main", ch, "tte"), "ho"), vals(tag_of("main", ch, "mst"), "ho")
        same = all(abs(a[k] - b[k]) < 1e-12 for k in a) and len(a) == len(b) and len(a) > 0
        put(f"TteEqMst{CHN[ch]}", "identical" if same else "NOT identical")
        a, b = vals(tag_of("main", ch, "tte40"), "ho"), vals(tag_of("main", ch, "mst40"), "ho")
        same = all(abs(a[k] - b[k]) < 1e-12 for k in a) and len(a) == len(b) and len(a) > 0
        put(f"TteFortyEqMstForty{CHN[ch]}", "identical" if same else "NOT identical")
    # blind geometric conditions under blockage
    blind = [m(S(k), "rlf") for k in ["d2", "t1", "elev"]]
    put("BlindRlfMin", fmt(min(blind), 1))
    put("BlindRlfMax", fmt(max(blind), 1))
    nets = [m(S(k), "out") for k in ["a3tf", "d2f", "t1f", "elevf"]]
    put("NetOutMin", fmt(min(nets), 3))
    put("NetOutMax", fmt(max(nets), 3))
    put("SubToneCellChanges", fmt(m(S("t1"), "ho") + m(S("t1"), "rlf") / 60.0, 2))
    # urban ranges
    for ch in ["urban", "denseurban"]:
        U = lambda p: tag_of("urban", ch, p)
        net = ["a3tf", "d2f", "t1f", "elevf", "tted2f"]
        hold = ["a3tfh", "d2fh", "t1fh", "elevfh"]
        matched = ["t1f40", "elevf40", "d2f4050"]
        for grp, keys in [("Net", net), ("Hold", hold), ("Matched", matched), ("Std", net + hold + matched)]:
            o = [m(U(k), "out") for k in keys if summaries().get(U(k))]
            hf = [m(U(k), "hof") for k in keys if summaries().get(U(k))]
            hh = [m(U(k), "ho") for k in keys if summaries().get(U(k))]
            if o:
                put(f"{CHN[ch]}{grp}OutMin", fmt(min(o), 2))
                put(f"{CHN[ch]}{grp}OutMax", fmt(max(o), 2))
                put(f"{CHN[ch]}{grp}HofMin", fmt(min(hf), 2))
                put(f"{CHN[ch]}{grp}HofMax", fmt(max(hf), 2))
                put(f"{CHN[ch]}{grp}HoMin", fmt(min(hh), 2))
                put(f"{CHN[ch]}{grp}HoMax", fmt(max(hh), 2))
        u40 = m(U("ut40"), "hof")
        hfs = [m(U(k), "hof") for k in net + hold + matched if summaries().get(U(k))]
        if hfs:
            rat = [u40 / x for x in hfs if x > 0]
            put(f"{CHN[ch]}HofRatioMin", fmt(min(rat), 1))
            put(f"{CHN[ch]}HofRatioMax", fmt(max(rat), 1))
        # outage advantage in seconds per terminal per 10 min against D2 + net
        r = RESULTS.get(f"{CHN[ch]}UtfortyVsDtwofOut")
        if r:
            put(f"{CHN[ch]}OutSecVsDtwof", fmt(-r["diff"] / 100.0 * 600.0, 2))


# ----------------------------------------------------------------- events
def events_rows(run_dir):
    with open(os.path.join(run_dir, "events.csv")) as fh:
        return list(csv.DictReader(fh))


def hof_mechanism():
    """M2: what the failed targets look like, and which branch chose them."""
    for ch in ["urban", "denseurban"]:
        for key in ["ut40", "mst40", "d2f", "t1f", "a3tf", "ut40_oldesc"]:
            el_f, l3_f, el_s, l3_s, esc_f, nf = [], [], [], [], 0, 0
            for r in summaries().get(tag_of("urban", ch, key), {}).values():
                for e in events_rows(r["_dir"]):
                    if e["kind"] != "ho":
                        continue
                    el, l3 = float(e["target_elev_deg"]), float(e["target_l3_db"])
                    if e["success"] == "0":
                        el_f.append(el)
                        l3_f.append(l3)
                        nf += 1
                        esc_f += e.get("trigger", "") in ("escape", "netesc")
                    else:
                        el_s.append(el)
                        l3_s.append(l3)
            if nf:
                pre = f"Hm{CHN[ch]}{POLN[key]}"
                put(pre + "FailElev", fmt(float(np.median(el_f)), 1))
                put(pre + "FailLthree", fmt(float(np.median(l3_f)), 1))
                put(pre + "OkElev", fmt(float(np.median(el_s)), 1))
                put(pre + "OkLthree", fmt(float(np.median(l3_s)), 1))
                put(pre + "FailEsc", fmt(100.0 * esc_f / nf, 0))


def per_class():
    """M8: per-class handovers and outage per terminal per 600 s, urban."""
    classes = ["static", "pedestrian", "vehicular", "iot", "hst", "aviation", "maritime"]
    out = {}
    for ch in ["urban", "suburban"]:
        exp = "urban" if ch == "urban" else "main"
        for key in ["ut40", "d2f", "mst40", "a3tf", "t1f"]:
            acc = defaultdict(lambda: [0.0, 0.0, 0])
            for r in summaries().get(tag_of(exp, ch, key), {}).values():
                for u in csv.DictReader(open(os.path.join(r["_dir"], "ue_summary.csv"))):
                    a = acc[u["class"]]
                    a[0] += float(u["ho_success"])
                    a[1] += float(u["outage_s"])
                    a[2] += 1
            for c in classes:
                if acc[c][2]:
                    nm = f"Cls{CHN[ch]}{POLN[key]}{c.capitalize()}"
                    put(nm + "Ho", fmt(acc[c][0] / acc[c][2], 1))
                    put(nm + "Out", fmt(acc[c][1] / acc[c][2], 1))
                    put(nm + "N", str(acc[c][2]))
            out[(ch, key)] = dict(acc)
    return out


def table_class(pc):
    classes = [("static", "Static"), ("pedestrian", "Pedestrian"), ("iot", "IoT"), ("vehicular", "Vehicle"),
               ("hst", "Train"), ("aviation", "Aircraft"), ("maritime", "Ship")]
    keys = [("d2f", "D2 + net"), ("mst40", "MST-40"), ("ut40", r"\textbf{UT-CHO-40}")]
    L = [r"\begin{table}[!t]", r"\caption{Urban environment by terminal class: handovers and outage (s) per "
         r"terminal and 600\,s run, pooled over seeds. Aircraft and ships see the sky in line of sight and "
         r"trains the rural table.}", r"\label{tab:class}", r"\centering", r"\scriptsize",
         r"\setlength{\tabcolsep}{2.4pt}", r"\begin{tabular}{@{}lr" + "rr" * len(keys) + r"@{}}", r"\toprule",
         r"Class & $n$ & " + " & ".join(r"\multicolumn{2}{c}{" + lab + "}" for _, lab in keys) + r" \\",
         "".join(f"\\cmidrule(l){{{3 + 2 * i}-{4 + 2 * i}}}" for i in range(len(keys))),
         r" & & " + " & ".join(["HO & Out."] * len(keys)) + r" \\", r"\midrule"]
    for c, lab in classes:
        a0 = pc.get(("urban", "ut40"), {}).get(c)
        if not a0:
            continue
        cells = [str(a0[2])]
        for k, _ in keys:
            a = pc.get(("urban", k), {}).get(c, [0, 0, 0])
            cells += [fmt(a[0] / a[2], 1) if a[2] else "--", fmt(a[1] / a[2], 1) if a[2] else "--"]
        L.append(lab + " & " + " & ".join(cells) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(PAPER, "table_class.tex"), "w").write("\n".join(L) + "\n")


def trigger_macros():
    """Escape share and target elevation from the trigger column."""
    for exp, chs in [("main", ["los", "suburban"]), ("urban", ["urban", "denseurban"])]:
        for ch in chs:
            for key in ["mst", "mst40", "ut", "ut40", "d2f", "t1f", "a3tf", "elevf"]:
                tot = esc = 0
                el = []
                for r in summaries().get(tag_of(exp, ch, key), {}).values():
                    for e in events_rows(r["_dir"]):
                        if e["kind"] != "ho":
                            continue
                        tot += 1
                        esc += e["trigger"] in ("escape", "netesc")
                        el.append(float(e["target_elev_deg"]))
                if tot:
                    put(f"Esc{CHN[ch]}{POLN[key]}", fmt(100.0 * esc / tot, 0))
                    put(f"TgtElev{CHN[ch]}{POLN[key]}", fmt(float(np.median(el)), 0))


# ----------------------------------------------------------------- accuracy (exit prediction)
RE, MU, H550 = 6378137.0, 3.986004418e14, 550e3
OMEGA = 7.2921159e-5


def lam(e_deg, h):
    e = math.radians(e_deg)
    return math.pi / 2 - e - math.asin(RE / (RE + h) * math.cos(e))


def rot_kappa(lat_deg=30.0, inc_deg=53.0, h=H550):
    """v_g / v_g,rel: inertial over Earth-relative ground speed of the sub-satellite point."""
    rs = RE + h
    vg = math.sqrt(MU / rs) * RE / rs
    lat = math.radians(lat_deg)
    s_az = min(1.0, math.cos(math.radians(inc_deg)) / math.cos(lat))  # heading from north
    north, east = vg * math.sqrt(1 - s_az ** 2), vg * s_az
    rel = math.hypot(north, east - OMEGA * RE * math.cos(lat))
    return vg / rel


def sigma_analytic(peak_deg, sa, sc, sr, eq=25.0, h=H550, kappa=None):
    """First-order exit sigma, eq. (sigma), with the exact spherical x_e and the
    Earth-rotation factor kappa = v_g / v_g,rel on every term."""
    rs = RE + h
    v = math.sqrt(MU / rs)
    k = rot_kappa(h=h) if kappa is None else kappa
    v = v / k
    vg = v * RE / rs
    rq = RE * lam(eq, h)
    y = RE * lam(min(peak_deg, 89.9), h)
    if y >= rq:
        return float("nan")
    # spherical: cos(rho_Q/R) = cos(x_e/R) cos(y/R)
    xe = RE * math.acos(min(1.0, math.cos(rq / RE) / math.cos(y / RE)))
    psiq = math.asin(RE / rs * math.cos(math.radians(eq)))
    a = (sa / v) ** 2
    c = (y / xe) ** 2 * (sc * RE / (rs * vg)) ** 2
    r = (rq / xe * (RE / rs) * math.tan(psiq)) ** 2 * (sr / vg) ** 2
    return math.sqrt(a + c + r)


def acc_block():
    res = {}
    for key, sa in [("sa0", 0), ("sa1000", 1000), ("sa5000", 5000), ("sa20000", 20000), ("off5000", 5000)]:
        rows = []
        for f in glob.glob(os.path.join(RUNS, f"acc/{key}/seed*/tte_accuracy.csv")):
            rows += list(csv.DictReader(open(f)))
        rows = [r for r in rows if 0 < float(r["true_tte_s"]) < 890]
        if not rows:
            continue
        pred = np.array([float(r["pred_tte_s"]) for r in rows])
        err = np.array([float(r["pred_tte_s"]) - float(r["true_tte_s"]) for r in rows])
        sig = np.array([float(r["sigma_pred_s"]) for r in rows])
        res[key] = (rows, err, sig, pred)
        tag = {"sa0": "Zero", "sa1000": "One", "sa5000": "Five", "sa20000": "Twenty", "off5000": "OffFive"}[key]
        pos = pred > 0
        put(f"Acc{tag}N", str(len(err)))
        put(f"Acc{tag}Rms", fmt(float(np.sqrt(np.mean(err ** 2))), 2))
        put(f"Acc{tag}RmsPos", fmt(float(np.sqrt(np.mean(err[pos] ** 2))), 2))
        put(f"Acc{tag}NZero", str(int((~pos).sum())))
        put(f"Acc{tag}Ptail", fmt(float(np.percentile(np.abs(err[pos]), 95)), 2))
        if (~pos).any():
            pk = np.array([float(r["peak_elev_deg"]) for r in rows])
            put(f"Acc{tag}ZeroPeakMin", fmt(float(pk[~pos].min()), 0))
            put(f"Acc{tag}ZeroPeakMax", fmt(float(pk[~pos].max()), 0))
        if sa > 0:
            m = pos & (sig > 0)
            z = err[m] / sig[m]
            put(f"Acc{tag}Zsd", fmt(float(np.std(z)), 2))
            put(f"Acc{tag}ZLate", fmt(100.0 * float(np.mean(z > 2.0)), 1))
            put(f"Acc{tag}ZEarly", fmt(100.0 * float(np.mean(z < -2.0)), 1))
            put(f"Acc{tag}SigMed", fmt(float(np.median(sig[m])), 2))
        if "compute_us" in rows[0]:
            us = np.array([float(r["compute_us"]) for r in rows])
            put(f"Acc{tag}Us", fmt(float(np.median(us)), 0))
    if "sa5000" in res:
        rows, err, sig, pred = res["sa5000"]
        pk = np.array([float(r["peak_elev_deg"]) for r in rows])
        an = np.array([sigma_analytic(x, 5000.0, 0.043 * 5000.0, 0.070 * 5000.0) for x in pk])
        m = (sig > 0) & np.isfinite(an) & (pred > 0)
        ratio = sig[m] / an[m]
        put("AccRatioMed", fmt(float(np.median(ratio)), 2))
        put("AccRatioLo", fmt(float(np.percentile(ratio, 5)), 2))
        put("AccRatioHi", fmt(float(np.percentile(ratio, 95)), 2))
        v = math.sqrt(MU / (RE + H550))
        hi = (pk >= 40) & m
        put("AccSigHiPeak", fmt(float(np.median(sig[hi])), 2))
        rel = sig[(pk < 30) & m] / (5000.0 / v)
        put("AccRelLoPeakMed", fmt(float(np.median(rel)), 2))
        put("AccRelLoPeakMax", fmt(float(np.max(rel)), 1))
        put("AlongOnly", fmt(5000.0 / v, 2))
        put("SecPerKm", fmt(1000.0 / v, 3))
        put("SigmaAlongTwenty", fmt(20000.0 / v, 1))
    # constants of the analysis
    n = math.sqrt(MU / (RE + H550) ** 3)
    put("OmegaRatio", fmt(OMEGA / n, 3))
    # Earth rotation: relative to the rotating ground the sub-satellite
    # point moves at v_g,rel, not v_g; for a prograde shell both ascending and
    # descending passes have the same eastward component, so at the service
    # latitude the factor kappa = v_g / v_g,rel is one number per shell
    put("RotKappa", fmt(rot_kappa(), 3))
    put("RotBound", fmt(100.0 * OMEGA / n, 0))
    return res


# ----------------------------------------------------------------- ephemeris experiments
def late_share(run_dir):
    ev = [e for e in events_rows(run_dir) if e["kind"] == "ho"]
    late = sum(1 for e in ev if e["true_serv_exit_s"] not in ("", None) and float(e["true_serv_exit_s"]) <= 0.0)
    return 100.0 * late / max(1, len(ev))


def eph_block(prefix):
    s = summaries()
    out = defaultdict(dict)
    for tag, rows in s.items():
        parts = tag.split("/")
        if parts[0] != prefix:
            continue
        _, lvl, pol = parts
        om, oh, n = ci([METRIC["out"](r) for r in rows.values()])
        hm, hh, _ = ci([METRIC["ho"](r) for r in rows.values()])
        lm, lh, _ = ci([late_share(r["_dir"]) for r in rows.values()])
        out[pol][lvl] = (om, oh, hm, hh, lm, lh, n)
    return out


EPH_NAMES = {"mst": "Mst", "tte": "Tte", "ut": "Ut", "t1": "Tone", "d2": "Dtwo", "elev": "Elev", "tte-under": "Under"}


def eph_macros(out):
    words = {"sa0": "Zero", "sa1000": "One", "sa2000": "Two", "sa5000": "Five", "sa10000": "Ten",
             "sa20000": "Twenty"}
    for pol, d in out.items():
        for lvl, v in d.items():
            k = f"Eph{EPH_NAMES[pol]}{words[lvl]}"
            put(k + "Late", fmt(v[4], 1))
            put(k + "LateCi", fmt(v[5], 1))
            put(k + "Out", fmt(v[0], 3))
            put(k + "Ho", fmt(v[2], 3))
    t = out.get("tte", {})
    if "sa0" in t and "sa20000" in t:
        put("EphTteHoRiseTwenty", fmt(100 * (t["sa20000"][2] / t["sa0"][2] - 1), 0))
    m = out.get("mst", {})
    if "sa5000" in m and m["sa5000"][4] > 0 and "sa0" in m:
        om, _, hm, _, lm = m["sa5000"][:5]
        put("EphMstFiveSecPerLate", fmt(60.0 * (om - m["sa0"][0]) / 100.0 / (hm * lm / 100.0), 2))
    for pol in ["d2", "elev", "t1"]:
        d = out.get(pol, {})
        if "sa0" in d and "sa20000" in d:
            put(f"Eph{EPH_NAMES[pol]}OutRiseTwenty", fmt(d["sa20000"][0] - d["sa0"][0], 3))
            put(f"Eph{EPH_NAMES[pol]}HoRiseTwenty", fmt(100 * (d["sa20000"][2] / d["sa0"][2] - 1), 0))


def emp_macros(out):
    words = {"age0-6": "Six", "age6-12": "Twelve", "age12-24": "Day", "age24-72": "Three"}
    for pol, d in out.items():
        for lvl, v in d.items():
            k = f"Emp{EPH_NAMES[pol]}{words[lvl]}"
            put(k + "Late", fmt(v[4], 1))
            put(k + "Out", fmt(v[0], 3))
            put(k + "Ho", fmt(v[2], 3))


def ephblk_block():
    """M9: measured 0-6 h errors under blockage, every policy, against its exact-ephemeris twin."""
    out = {}
    for ch in ["suburban", "urban"]:
        exp0 = "main" if ch == "suburban" else "urban"
        for pol in ["a3tf", "d2f", "t1f", "elevf", "t1f40", "elevf40", "d2f4050", "mst", "mst40", "tte", "tte40",
                    "ut", "ut40"]:
            tag = f"ephblk/{ch}/{pol}"
            if not summaries().get(tag):
                continue
            twin = tag_of(exp0, ch, pol if pol not in ("tte", "tte40") else {"tte": "mst", "tte40": "mst40"}[pol])
            row = {}
            for mk in ["ho", "out", "hof", "rlf"]:
                row[mk] = stat(tag, mk)
                r = paired(tag, twin, mk)
                row[mk + "_d"] = r
            lm = ci([late_share(r["_dir"]) for r in summaries()[tag].values()])
            row["late"] = lm
            # the late share of the same policy with an exact ephemeris
            row["late_exact"] = ci([late_share(r["_dir"]) for r in summaries().get(twin, {}).values()])
            out[(ch, pol)] = row
            pre = f"Blk{CHN[ch]}{POLN[pol]}"
            put(pre + "Out", fmt(row["out"][0], 3))
            put(pre + "Ho", fmt(row["ho"][0], 3))
            put(pre + "Late", fmt(lm[0], 2))
            put(pre + "LateExact", fmt(row["late_exact"][0], 2))
            if row["out_d"]:
                put(pre + "DOut", fmt(row["out_d"]["diff"], 3))
                put(pre + "DOutCi", fmt(row["out_d"]["ci"], 3))
                put(pre + "DOutP", pfmt(row["out_d"]["p"]))
            if row["ho_d"]:
                put(pre + "DHo", fmt(row["ho_d"]["diff"], 3))
    return out


def realtle_block():
    out = {}
    for ch in ["los", "suburban"]:
        for pol in ["a3tf", "d2f", "t1f", "elevf", "mst", "mst40", "tte", "tte40", "ut", "ut40", "tteage"]:
            tag = f"realtle/{ch}/{pol}"
            if not summaries().get(tag):
                continue
            row = {mk: stat(tag, mk) for mk in ["ho", "out", "hof", "rlf", "snr"]}
            row["late"] = ci([late_share(r["_dir"]) for r in summaries()[tag].values()])
            out[(ch, pol)] = row
            pre = f"Tle{CHN[ch]}{POLN[pol]}"
            for mk in ["ho", "out", "hof", "snr"]:
                put(pre + MNAME[mk], fmt(row[mk][0], DIG[mk]))
            put(pre + "Late", fmt(row["late"][0], 2))
            put(pre + "LateCi", fmt(row["late"][1], 2))
        n = [num(r, "n_sats") for r in summaries().get(f"realtle/{ch}/mst", {}).values()]
        if n:
            put("TleSats", str(int(np.median(n))))
    # measured initial error of the real network copy
    d0 = []
    for f in glob.glob(os.path.join(RUNS, "realtle/los/mst/seed*/eph_offsets.csv")):
        for r in csv.DictReader(open(f)):
            d0.append([float(r["d0_r_m"]), float(r["d0_a_m"]), float(r["d0_c_m"]), float(r["dT_a_m"])])
    if d0:
        d0 = np.array(d0)
        put("TleRunSigA", fmt(float(np.percentile(np.abs(d0[:, 1]), 68.27)) / 1000, 2))
        put("TleRunSigR", fmt(float(np.percentile(np.abs(d0[:, 0]), 68.27)), 0))
        put("TleRunGrowA", fmt(float(np.percentile(np.abs(d0[:, 3] - d0[:, 1]), 68.27)), 0))
    for ch in ["los", "suburban"]:
        r = paired(f"realtle/{ch}/ut", f"realtle/{ch}/mst", "out")
        if r:
            put(f"Tle{CHN[ch]}UtVsMstOut", fmt(r["diff"], 3))
            put(f"Tle{CHN[ch]}UtVsMstOutCi", fmt(r["ci"], 3))
            put(f"Tle{CHN[ch]}UtVsMstOutP", pfmt(r["p"]))
        a = {k: late_share(x["_dir"]) for k, x in summaries().get(f"realtle/{ch}/ut", {}).items()}
        b = {k: late_share(x["_dir"]) for k, x in summaries().get(f"realtle/{ch}/mst", {}).items()}
        ks = sorted(set(a) & set(b))
        if ks:
            put(f"Tle{CHN[ch]}UtVsMstLateP", pfmt(wilcoxon_p([a[k] - b[k] for k in ks])))
            put(f"Tle{CHN[ch]}UtVsMstLateLower", str(sum(a[k] < b[k] for k in ks)))
    return out


# ----------------------------------------------------------------- recent-work baselines (Oct 2026)
RECENT = [("tted2f", "noreen2026location", "Noreen and Schotten"), ("geof", "lee2026geodesic", "Lee and Kim"),
          ("wangf", "wang2026mdp", "Wang et al."), ("zabf", "zaboub2026multi", "Zaboub et al.")]
RECENT_ENV = [("main", "los"), ("main", "suburban"), ("urban", "urban"), ("urban", "denseurban")]


def recent_block():
    """UT-CHO against the four re-implemented recent methods: one Holm family for handover rate
    (line of sight, suburban; reference UT-CHO) and one for outage (urban, dense urban; reference UT-CHO-40)."""
    fam = []
    for key, _, _ in RECENT:
        for exp, ch in RECENT_ENV:
            ref, mk = ("ut", "ho") if exp == "main" else ("ut40", "out")
            r = paired(tag_of(exp, ch, ref), tag_of(exp, ch, key), mk)
            if r:
                fam.append((f"Rec{CHN[ch]}{POLN[key]}", r))
    adj = holm([r["p"] for _, r in fam])
    for (nm, r), pa in zip(fam, adj):
        r["holm"] = pa
        RESULTS[nm] = r
        put(nm + "Diff", fmt(r["diff"], 3))
        put(nm + "Ci", fmt(r["ci"], 3))
        put(nm + "Rel", fmt(r["rel"], 1))
        put(nm + "RelAbs", fmt(abs(r["rel"]), 1))
        put(nm + "P", pfmt(pa))
    # summary macros: range of the relative handover cut against the recent methods
    for ch in ["los", "suburban"]:
        rels = [RESULTS[f"Rec{CHN[ch]}{POLN[k]}"]["rel"] for k, _, _ in RECENT if f"Rec{CHN[ch]}{POLN[k]}" in RESULTS]
        if rels:
            put(f"Rec{CHN[ch]}HoCutMin", fmt(-max(rels), 1))
            put(f"Rec{CHN[ch]}HoCutMax", fmt(-min(rels), 1))
    for ch in ["urban", "denseurban"]:
        ds = [RESULTS[f"Rec{CHN[ch]}{POLN[k]}"]["diff"] for k, _, _ in RECENT if f"Rec{CHN[ch]}{POLN[k]}" in RESULTS]
        if ds:
            put(f"Rec{CHN[ch]}OutCutMin", fmt(-max(ds), 2))
            put(f"Rec{CHN[ch]}OutCutMax", fmt(-min(ds), 2))
    for ch in ["urban", "denseurban"]:
        ds = [-RESULTS[f"Rec{CHN[ch]}{POLN[k]}"]["diff"] for k, _, _ in RECENT
              if k != "geof" and f"Rec{CHN[ch]}{POLN[k]}" in RESULTS]
        if ds:
            put(f"Rec{CHN[ch]}OutWinMin", fmt(min(ds), 2))
            put(f"Rec{CHN[ch]}OutWinMax", fmt(max(ds), 2))
        g, u = stat(tag_of("urban", ch, "geof"), "ho")[0], stat(tag_of("urban", ch, "ut40"), "ho")[0]
        if u == u and u > 0:
            put(f"Rec{CHN[ch]}GeofHoOverUt", fmt(g / u, 1))
    for key, _, _ in RECENT + [("mst", "", ""), ("ut", "", ""), ("ut40", "", "")]:
        for exp, ch in RECENT_ENV:
            for mk in ["ho", "out", "snr", "hof", "rlf"]:
                m, h, n = stat(tag_of(exp, ch, key), mk)
                if n:
                    put(f"Rec{CHN[ch]}{POLN[key]}{MNAME[mk]}", fmt(m, DIG[mk]))
        for ch in ["los", "suburban"]:
            rows = summaries().get(f"realtle/{ch}/{key}")
            if rows:
                m, h, n = ci([late_share(r["_dir"]) for r in rows.values()])
                put(f"Tle{CHN[ch]}{POLN[key]}Late", fmt(m, 2))
                put(f"Tle{CHN[ch]}{POLN[key]}LateCi", fmt(h, 2))


def table_recent():
    L = [r"\begin{table*}[!t]", r"\caption{Comparison With Recent Methods}", r"\label{tab:recent}", r"\centering",
         r"\footnotesize", r"\setlength{\tabcolsep}{2.4pt}",
         r"\begin{tabular}{@{}l" + "r" * 2 + r"@{\hspace{6pt}}" + "r" * 2 + r"@{\hspace{6pt}}" + "r" * 3
         + r"@{\hspace{6pt}}" + "r" * 3 + r"@{\hspace{6pt}}r@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{Line of sight} & \multicolumn{2}{c}{Suburban} & \multicolumn{3}{c}{Urban}"
         r" & \multicolumn{3}{c}{Dense urban} & Real shell \\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}\cmidrule(l){12-12}",
         r"Policy & HO & $\Delta$HO\,\% & HO & $\Delta$HO\,\% & HO & Out.\,\% & $\Delta$Out. & HO & Out.\,\% & $\Delta$Out. & Late\,\% \\",
         r"\midrule"]

    def dcell(ch, key):
        r = RESULTS.get(f"Rec{CHN[ch]}{POLN[key]}")
        if r is None:
            return "--"
        v = -r["rel"] if ch in ("los", "suburban") else -r["diff"]
        d = 1 if ch in ("los", "suburban") else 2
        t = fmt(v, d)
        t = "$-$" + t[1:] if t.startswith("-") else t
        return t + (r"$^{\ast}$" if r["holm"] < ALPHA else "")

    rows = [(k, f"{LAB[k]}~\\cite{{{c}}}") for k, c, _ in RECENT] + [None, ("mst", "MST~\\cite{papapetrou2004satellite}"),
                                                                 None, ("ut", r"\textbf{UT-CHO}"), ("ut40", r"\textbf{UT-CHO-40}")]
    for row in rows:
        if row is None:
            L.append(r"\midrule")
            continue
        key, lab = row
        c = []
        for exp, ch in RECENT_ENV:
            T = tag_of(exp, ch, key)
            if exp == "main":
                c += [cell(T, "ho", 3), dcell(ch, key) if key in dict((k, 1) for k, _, _ in RECENT) else ""]
            else:
                c += [cell(T, "ho", 2), cell(T, "out", 2),
                      dcell(ch, key) if key in dict((k, 1) for k, _, _ in RECENT) else ""]
        rows_t = summaries().get(f"realtle/los/{key}")
        late = fmt(ci([late_share(r["_dir"]) for r in rows_t.values()])[0], 2) if rows_t else "--"
        L.append(lab + " & " + " & ".join(c) + f" & {late} \\\\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize Means over " + macros.get("NSeedMain", "ten") + r" seeds. HO per terminal-minute, outage in \%. "
          r"$\Delta$HO: handover-rate cut of UT-CHO relative to the row, seed-paired mean in \%; $\Delta$Out.: outage of "
          r"the row minus that of UT-CHO-40 in percentage points; $^{\ast}$: Holm-adjusted $p<0.05$ over the "
          r"family of all sixteen comparisons. Late: share of handovers decided after the true exit on the real "
          r"53$^\circ$ shell, line of sight. Every row carries the same measurement safety net.}",
          r"\end{table*}"]
    open(os.path.join(PAPER, "table_recent.tex"), "w").write("\n".join(L) + "\n")


def table_eph_std(blk, tle):
    """Ephemeris error against every policy: measured 0-6 h errors under blockage, and the real shell."""
    rows = [("a3tf", "A3-tuned + net"), ("d2f", "D2 + net"), ("t1f", "T1 + net"), ("elevf", "Elevation + net"),
            ("mst", "MST"), ("tte", "MST + margin"), ("tteage", "MST + age margin"), ("ut", r"\textbf{UT-CHO}"), ("mst40", "MST-40"),
            ("tte40", "MST-40 + margin"),
            ("ut40", r"\textbf{UT-CHO-40}")]
    L = [r"\begin{table*}[!t]",
         r"\caption{Ephemeris error against every policy, ten seeds. Left: each satellite's error drawn from the "
         r"measured 0 to 6\,h element-set errors; $\Delta$Out.\ is the outage minus that of the same policy with an "
         r"exact ephemeris, seed-paired (percentage points). Right: the real 53$^\circ$ Starlink shell, truth from "
         r"the supplemental element sets and the network's copy from the public ones, both SGP4.}",
         r"\label{tab:ephstd}", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.4pt}",
         r"\begin{tabular}{@{}lrrrrrr@{\hspace{8pt}}rrrrr@{}}", r"\toprule",
         r" & \multicolumn{3}{c}{Suburban, measured error} & \multicolumn{3}{c}{Urban, measured error} & "
         r"\multicolumn{2}{c}{Real shell, LOS} & \multicolumn{3}{c}{Real shell, suburban} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}\cmidrule(l){10-12}",
         r"Policy & Late\,\% & Out.\,\% & $\Delta$Out. & Late\,\% & Out.\,\% & $\Delta$Out. & Late\,\% & Out.\,\% & HO & Late\,\% & Out.\,\% \\",
         r"\midrule"]
    for key, lab in rows:
        c = []
        for ch in ["suburban", "urban"]:
            b = blk.get((ch, key))
            if b:
                d = b["out_d"]
                c += [fmt(b["late"][0], 2), fmt(b["out"][0], 3), (fmt(d["diff"], 3) if d else "--")]
            else:
                c += ["--"] * 3
        t = tle.get(("los", key))
        c += [fmt(t["late"][0], 2), fmt(t["out"][0], 3)] if t else ["--"] * 2
        t = tle.get(("suburban", key))
        c += [fmt(t["ho"][0], 3), fmt(t["late"][0], 2), fmt(t["out"][0], 3)] if t else ["--"] * 3
        L.append(lab + " & " + " & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize Late: handovers decided after the serving satellite left the coverage contour on "
          r"its true orbit. The D2 references, the T1 epoch and the elevation the elevation trigger evaluates all "
          r"come from the network's erroneous copy; A3 uses none of them. Under blockage every late share "
          r"except MST's and UT-CHO's equals that of the same policy with an exact ephemeris to the printed "
          r"precision; these are mostly escapes after the exit. MST + margin and MST + age margin ran on the real shell "
          r"only (--).}",
          r"\end{table*}"]
    open(os.path.join(PAPER, "table_ephstd.tex"), "w").write("\n".join(L) + "\n")


# ----------------------------------------------------------------- lower bound, sensitivity, density
def lb_block():
    out = {}
    for tag, rows in summaries().items():
        if not tag.startswith("lb/"):
            continue
        name = tag.split("/")[1]
        gaps, at = [], []
        for r in rows.values():
            for u in csv.DictReader(open(os.path.join(r["_dir"], "ue_summary.csv"))):
                lb = float(u["lb_handovers"])
                # one definition: completed handovers plus failed ones, against the bound
                ho = float(u["ho_success"]) + float(u["hof"])
                gaps.append(ho - lb)
                at.append(ho <= lb)
        out[name] = (ci(gaps), ci([num(r, "lb_handovers_per_ue") for r in rows.values()]),
                     float(np.mean(at)) if at else float("nan"), len(rows))
    lab = {"a3t": "Atuned", "d2": "Dtwo", "t1": "Tone", "elev": "Elev", "mst": "Mst", "mst_prepr": "Mstprepr",
           "ut": "Ut", "mst_k16": "Mstk"}
    for k, v in out.items():
        if k in lab:
            put(f"Lb{lab[k]}Gap", fmt(v[0][0], 2))
            put(f"Lb{lab[k]}GapCi", fmt(v[0][1], 2))
            put(f"Lb{lab[k]}AtBound", fmt(100 * v[2], 0))
    if out:
        v = next(iter(out.values()))
        put("LbBound", fmt(v[1][0], 2))
        put("LbN", str(v[3]))
        allb = [num(r, "lb_handovers_per_ue") for t, rows in summaries().items() if t.startswith("lb/")
                for r in rows.values()]
        put("LbUniform", "uniform" if (max(allb) - min(allb)) < 1e-9 else "varies")
    return out


def sens_block():
    out = {}
    for tag, rows in summaries().items():
        p = tag.split("/")
        if p[0] != "sens":
            continue
        out[(p[1], p[2])] = {mk: stat(tag, mk) for mk in ["ho", "snr", "out"]}
    return out


def table_lb_sens(lb, sens):
    names = [("a3t", "A3-tuned"), ("d2", "D2"), ("t1", "T1"), ("elev", "Elevation"), ("mst", "MST"),
             ("ut", r"\textbf{UT-CHO}"), ("ut40", r"\textbf{UT-CHO-40}")]
    L = [r"\begin{table}[!t]", r"\caption{Optimality gap and sensitivity (ten seeds; gap with its 95\,\% CI).}",
         r"\label{tab:sens}", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.4pt}",
         r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
         r"Policy & HO above & \multicolumn{3}{c}{HO/terminal-min (SNR, dB)} \\", r"\cmidrule(l){3-5}",
         r"       & bound/term. & $e_{\min}=10^\circ$ & $e_{\min}=40^\circ$ & 1200~km \\", r"\midrule"]
    for k, lab in names:
        g = lb.get(k)
        gs = f"{g[0][0]:.2f}\\,$\\pm$\\,{g[0][1]:.2f}" if g else "--"
        c = []
        for case in ["minel10", "minel40", "leo1200"]:
            v = sens.get((case, k))
            c.append(f"{v['ho'][0]:.2f} ({v['snr'][0]:.1f})" if v else "--")
        L.append(f"{lab} & {gs} & " + " & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\\[2pt]{\scriptsize Bound: offline greedy cover on the true orbits, static terminals, 600\,s; "
          r"handovers counted as completed plus failed.}", r"\end{table}"]
    open(os.path.join(PAPER, "table_sens.tex"), "w").write("\n".join(L) + "\n")
    cases = {"minel10": "Ten", "minel40": "Forty", "leo1200": "Polar"}
    lab = {"a3t": "Atuned", "d2": "Dtwo", "t1": "Tone", "elev": "Elev", "mst": "Mst", "mst40": "Mstforty",
           "ut": "Ut", "ut40": "Utforty"}
    for (case, k), v in sens.items():
        if k in lab:
            b = f"Sens{cases[case]}{lab[k]}"
            put(b + "Ho", fmt(v["ho"][0], 2))
            put(b + "Snr", fmt(v["snr"][0], 1))
            put(b + "Out", fmt(v["out"][0], 2))


def density_macros():
    for tag, rows in summaries().items():
        p = tag.split("/")
        if p[0] != "density":
            continue
        sh, ch, name = p[1], p[2], p[3]
        nm = f"Dens{'Half' if sh == 'w792' else 'Double'}{CHN[ch]}{POLN[name]}"
        put(nm + "Ho", fmt(stat(tag, "ho")[0], 2))
        put(nm + "Out", fmt(stat(tag, "out")[0], 2))
    for sh, w in [("w792", "Half"), ("w3168", "Double")]:
        for ch in ["los", "suburban"]:
            r = paired(f"density/{sh}/{ch}/ut40", f"density/{sh}/{ch}/t1f", "ho")
            if r:
                put(f"Dens{w}{CHN[ch]}UtfortyVsTonefHo", fmt(r["diff"], 2))
                put(f"Dens{w}{CHN[ch]}UtfortyVsTonefHoP", pfmt(r["p"]))


def losmodel_block():
    """M16/M8: LOS process and terminal environment, urban, against the reference model."""
    out = {}
    for v in ["mk50", "mk200", "mk50a3", "noclass"]:
        for pol in ["a3tf", "d2f", "t1f", "elevf", "t1f40", "mst40", "ut", "ut40"]:
            tag = f"losmodel/{v}/{pol}"
            if not summaries().get(tag):
                continue
            out[(v, pol)] = {mk: stat(tag, mk) for mk in ["ho", "out", "hof", "rlf"]}
            nm = f"Lm{ {'mk50': 'Mfifty', 'mk200': 'Mtwohundred', 'mk50a3': 'Mthree', 'noclass': 'Noclass'}[v]}{POLN[pol]}"
            put(nm + "Out", fmt(out[(v, pol)]["out"][0], 2))
            put(nm + "Ho", fmt(out[(v, pol)]["ho"][0], 2))
            put(nm + "Hof", fmt(out[(v, pol)]["hof"][0], 2))
        # does UT-CHO-40 keep its outage rank against D2 + net under this model?
        r = paired(f"losmodel/{v}/ut40", f"losmodel/{v}/d2f", "out")
        vn = {'mk50': 'Mfifty', 'mk200': 'Mtwohundred', 'mk50a3': 'Mthree', 'noclass': 'Noclass'}[v]
        if r:
            put(f"Lm{vn}UtfortyVsDtwofOut", fmt(r["diff"], 3))
            put(f"Lm{vn}UtfortyVsDtwofOutCi", fmt(r["ci"], 3))
            put(f"Lm{vn}UtfortyVsDtwofOutP", pfmt(r["p"]))
        r = paired(f"losmodel/{v}/ut40", f"losmodel/{v}/mst40", "out")
        if r:
            put(f"Lm{vn}UtfortyVsMstfortyOut", fmt(r["diff"], 3))
            put(f"Lm{vn}UtfortyVsMstfortyOutP", pfmt(r["p"]))
    return out


def table_losmodel(lm):
    if not lm:
        return
    pols = [("a3tf", "A3-tuned + net"), ("d2f", "D2 + net"), ("t1f", "T1 + net"),
            ("mst40", "MST-40"), ("ut40", r"\textbf{UT-CHO-40}")]
    vs = [("ref", "Reference"), ("noclass", "All in clutter"), ("mk50", "Markov 50\\,m"),
          ("mk200", "Markov 200\\,m"), ("mk50a3", "Markov 50\\,m, 3$^\\circ$")]
    L = [r"\begin{table}[!t]", r"\caption{Urban outage (\%) under other line-of-sight processes, ten seeds. "
         r"Reference: redraw after 50\,m or a new 10$^\circ$ bin, aircraft and ships in line of sight. "
         r"All in clutter: every class in the urban table. Markov: first-order process with the given "
         r"correlation distance and an angular scale of 10$^\circ$ unless stated.}", r"\label{tab:losmodel}",
         r"\centering", r"\scriptsize", r"\setlength{\tabcolsep}{2.0pt}", r"\begin{tabular}{@{}l" + "r" * len(vs) + r"@{}}",
         r"\toprule", "Policy & " + " & ".join(l for _, l in vs) + r" \\", r"\midrule"]
    for k, lab in pols:
        c = []
        for v, _ in vs:
            if v == "ref":
                m = stat(tag_of("urban", "urban", k), "out")[0]
            else:
                m = lm.get((v, k), {}).get("out", (float("nan"),))[0]
            c.append(fmt(m, 2))
        L.append(lab + " & " + " & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(PAPER, "table_losmodel.tex"), "w").write("\n".join(L) + "\n")


def reest_block():
    for ch in ["suburban", "urban"]:
        exp0 = "main" if ch == "suburban" else "urban"
        for pol in ["d2", "t1", "d2f", "t1f", "a3tf", "mst40", "ut40"]:
            tag = f"reest/{ch}/{pol}"
            if not summaries().get(tag):
                continue
            nm = f"Re{CHN[ch]}{POLN[pol]}"
            put(nm + "Out", fmt(stat(tag, "out")[0], 3))
            r = paired(tag, tag_of(exp0, ch, pol), "out")
            if r:
                put(nm + "DOut", fmt(r["diff"], 3))
        r = paired(f"reest/{ch}/ut40", f"reest/{ch}/d2f", "out")
        if r:
            put(f"Re{CHN[ch]}UtfortyVsDtwofOut", fmt(r["diff"], 3))
            put(f"Re{CHN[ch]}UtfortyVsDtwofOutP", pfmt(r["p"]))
        r = paired(f"reest/{ch}/ut40", f"reest/{ch}/mst40", "out")
        if r:
            put(f"Re{CHN[ch]}UtfortyVsMstfortyOut", fmt(r["diff"], 3))
            put(f"Re{CHN[ch]}UtfortyVsMstfortyOutP", pfmt(r["p"]))


# ----------------------------------------------------------------- TLE measurement
def tle_block():
    f = os.path.join(HERE, "tle", "tle_error_samples.csv")
    if not os.path.exists(f):
        return None
    all_rows = list(csv.DictReader(open(f)))
    rows = [r for r in all_rows if r["dt_min"] == "0"]
    E = np.array([[float(r["age_h"]), float(r["err_r_m"]), float(r["err_a_m"]), float(r["err_c_m"])] for r in rows])
    keep = np.all(np.abs(E[:, 1:]) <= 1e5, axis=1)
    put("TleOutliers", str(int((~keep).sum())))
    E = E[keep]
    s68 = lambda x: float(np.percentile(np.abs(x), 68.27))
    put("TleN", f"{len(E):,}".replace(",", "{,}"))
    put("TleSigR", fmt(s68(E[:, 1]), 0))
    put("TleSigA", fmt(s68(E[:, 2]) / 1000, 1))
    put("TleSigC", fmt(s68(E[:, 3]), 0))
    put("TleRatioCA", fmt(s68(E[:, 3]) / s68(E[:, 2]), 3))
    put("TleRatioRA", fmt(s68(E[:, 1]) / s68(E[:, 2]), 3))
    bins = [(0, 6, "Six"), (6, 12, "Twelve"), (12, 24, "Day"), (24, 72, "Three")]
    out = []
    for lo, hi, nm in bins:
        m = (E[:, 0] >= lo) & (E[:, 0] < hi)
        put(f"TleSigA{nm}", fmt(s68(E[m, 2]) / 1000, 1))
        put(f"TleN{nm}", str(int(m.sum())))
        put(f"TleBiasA{nm}", fmt(float(np.median(E[m, 2])) / 1000, 2))
        put(f"TlePosA{nm}", fmt(100.0 * float(np.mean(E[m, 2] > 0)), 0))
        out.append((lo, hi, s68(E[m, 2]), int(m.sum())))
    # along-track growth rate over 30 min against the drift a radial (semi-major-axis) error implies
    by = defaultdict(dict)
    for r in all_rows:
        by[r["catnr"]][r["dt_min"]] = r
    rates, rad = [], []
    for c, d in by.items():
        if "0" in d and "30" in d:
            a0, a30 = float(d["0"]["err_a_m"]), float(d["30"]["err_a_m"])
            if abs(a0) <= 1e5 and abs(a30) <= 1e5:
                rates.append((a30 - a0) / 1800.0)
                rad.append(float(d["0"]["err_r_m"]))
    n = math.sqrt(MU / (RE + H550) ** 3)
    if rates:
        put("TleRateA", fmt(s68(np.array(rates)), 2))
        put("TleRateFromR", fmt(1.5 * n * s68(np.array(rad)), 2))
    return E, out


# ----------------------------------------------------------------- figures
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.style.use(os.path.abspath(os.path.join(PAPER, "../templates/ieee.mplstyle")))
NAVY, BLUE, GOLD, GRAY, STEEL, RED = "#0E3A6B", "#2E78B5", "#C98A0B", "#4A4A4A", "#7FA3CC", "#C23B22"
# One color, marker and dash per policy family in every result figure (figures.md CVD palette, fixed order).
# UT-CHO takes the house orange that also marks "this work" in the diagrams.
C_BLUE, C_ORANGE, C_GREEN, C_GOLD = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
C_PURPLE, C_TEAL, C_RED, C_GREY, C_INK = "#8a5cd1", "#0f9ba8", "#d64545", "#7b8794", "#1f2933"
C_LIGHT = "#9fc3ee"
PSTY = {  # family: (color, marker, linestyle)
    "ut": (C_ORANGE, "D", "-"), "ut40": (C_ORANGE, "D", "--"), "mst": (C_BLUE, "o", "--"),
    "mst40": (C_BLUE, "o", "-."), "tte": (C_TEAL, "s", "-"), "tte-under": (C_RED, "x", "-."),
    "t1": (C_PURPLE, "v", ":"), "d2": (C_GREEN, "^", "-."), "elev": (C_GREY, "s", ":"), "a3": (C_GOLD, "o", "-")}
plt.rcParams.update({"font.family": "serif", "pdf.fonttype": 42, "ps.fonttype": 42, "legend.fontsize": 7,
                     "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})


def fig_accuracy(res):
    if not res:
        return
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.5))
    cols = {"sa0": C_INK, "sa1000": C_BLUE, "sa5000": C_TEAL, "sa20000": C_GREY}
    lss = {"sa0": "-", "sa1000": "--", "sa5000": "-.", "sa20000": ":"}
    for key in ["sa0", "sa1000", "sa5000", "sa20000"]:
        if key not in res:
            continue
        rows, err, sig, pred = res[key]
        x = np.sort(np.abs(err[pred > 0]))
        y = np.arange(1, len(x) + 1) / len(x)
        sa = int(key[2:])
        lab = "no error" if sa == 0 else f"$\\sigma_A$ = {sa / 1000:g} km"
        ax[0].plot(np.maximum(x, 1e-3), y, color=cols[key], ls=lss[key], lw=1.2, label=lab)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("|predicted - realized exit| (s)")
    ax[0].set_ylabel("CDF")
    ax[0].legend(loc="lower right", frameon=False)
    ax[0].set_title(r"(a) prediction error, samples with $\hat{\tau}>0$", fontsize=8)
    if "sa5000" in res:
        rows, err, sig, pred = res["sa5000"]
        sa, sc, sr = 5000.0, 0.043 * 5000.0, 0.070 * 5000.0
        v = math.sqrt(MU / (RE + H550))
        pk = np.array([float(r["peak_elev_deg"]) for r in rows])
        m = (sig > 0) & (pred > 0)
        ax[1].scatter(pk[m], sig[m] / (sa / v), s=5, color=C_LIGHT, alpha=0.6, lw=0, zorder=2, rasterized=True,
                      label="central difference (simulator)")
        g = np.linspace(25.5, 89.5, 300)
        ax[1].plot(g, [sigma_analytic(x, sa, sc, sr) / (sa / v) for x in g], color=C_INK, lw=1.3, zorder=3,
                   label="first-order analysis")
        ax[1].axhline(1.0, color=C_GREY, ls=":", lw=0.8)
        ytop = 2.5
        put("AccClipped", str(int(np.sum(sig[m] / (sa / v) > ytop))))
        ax[1].set_ylim(0.9, ytop)
        ax[1].set_xlabel("pass peak elevation (deg)")
        ax[1].set_ylabel(r"$\sigma_\tau / (\sigma_A/v)$")
        ax[1].legend(loc="upper right", frameon=False)
        ax[1].set_title(r"(b) exit sigma against pass geometry, $\sigma_A$ = 5 km", fontsize=8)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(os.path.join(FIG, "fig_accuracy.pdf"))
    plt.close(fig)


FAMS = [("a3", "A3 (offset, TTT)", C_GOLD, "o", "-"), ("elev", "Elevation + net (floor)", C_GREY, "s", ":"),
        ("d2", "D2 + net (thresholds)", C_GREEN, "^", "-."), ("t1", "T1 + net (exit contour)", C_PURPLE, "v", ":"),
        ("t1c", "T1 + net, prepared inside the contour", C_PURPLE, "x", "--"),
        ("mst", "MST (exit contour)", C_BLUE, "o", "--"), ("ut", "UT-CHO (exit contour)", C_ORANGE, "D", "-")]


def pareto_block():
    pts = defaultdict(list)
    for tag, rows in summaries().items():
        p = tag.split("/")
        if p[0] != "pareto":
            continue
        ch, name = p[1], p[2]
        fam = name.split("_")[0]
        if name.startswith("t1_w"):
            fam = "t1w"
        hm, hh, n = ci([METRIC["ho"](r) for r in rows.values()])
        sm, sh, _ = ci([METRIC["snr"](r) for r in rows.values()])
        om, oh, _ = ci([METRIC["out"](r) for r in rows.values()])
        pts[(ch, fam)].append((hm, sm, hh, sh, om, name, n, oh))
    return pts


def frontier(p, key=1, better=max):
    out, best = [], None
    for x in sorted(p):
        if best is None or (better(x[key], best) == x[key] and x[key] != best):
            out.append(x)
            best = x[key]
    return out


XMAX = 2.5


def fig_pareto(pts):
    if not pts:
        return
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.0))
    panels = [("los", 1, max, "mean SNR (dB)", "(a) line of sight: SNR"),
              ("suburban", 4, min, "outage (%)", "(b) suburban LOS process: outage")]
    for k, (ch, key, better, ylab, title) in enumerate(panels):
        for fam, lab, col, mk, ls in FAMS:
            p = sorted(x for x in pts.get((ch, fam), []) if x[0] <= XMAX)
            if not p:
                continue
            yerr = [x[3] if key == 1 else x[7] for x in p]
            ax[k].errorbar([x[0] for x in p], [x[key] for x in p], xerr=[x[2] for x in p], yerr=yerr,
                           fmt=mk, color=col, ms=3.6, lw=0.6, capsize=1.2,
                           mfc="white" if fam not in ("ut", "mst") else col, label=lab,
                           zorder=3 if fam == "ut" else 2)
            f = frontier(p, key, better)
            ax[k].plot([x[0] for x in f], [x[key] for x in f], color=col, ls=ls, lw=0.9)
        ax[k].set_xscale("log")
        lo = 0.15 if k == 0 else 0.6
        ax[k].set_xlim(lo, XMAX)
        ticks = [0.2, 0.3, 0.5, 1.0, 2.0] if k == 0 else [0.7, 1.0, 1.5, 2.0]
        ax[k].set_xticks(ticks)
        ax[k].set_xticklabels([f"{t:g}" for t in ticks])
        ax[k].minorticks_off()
        ax[k].set_xlabel("handovers per terminal-minute")
        ax[k].set_ylabel(ylab)
        ax[k].set_title(title, fontsize=8)
        ax[k].grid(True, which="major")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 1.0),
               handletextpad=0.4, columnspacing=1.0)
    fig.tight_layout(w_pad=1.5, rect=(0, 0, 1, 0.84))
    fig.savefig(os.path.join(FIG, "fig_pareto.pdf"))
    plt.close(fig)


def pareto_macros(pts):
    if not pts:
        return
    std = lambda ch: [x for (c, f), v in pts.items() if c == ch and f in ("a3", "elev", "d2", "t1", "t1c", "t1w")
                      for x in v]
    fam = lambda ch, f: {x[5]: x for x in pts.get((ch, f), [])}
    lo = std("los")
    put("FrontMinStdHoLos", fmt(min(x[0] for x in lo), 2))
    put("FrontMinUtHoLos", fmt(min(x[0] for x in fam("los", "ut").values()), 2))
    put("FrontMinMstHoLos", fmt(min(x[0] for x in fam("los", "mst").values()), 2))
    best = None
    for name, x in fam("los", "ut").items():
        cand = [o[1] for o in lo if o[0] <= x[0]]
        if cand:
            g = x[1] - max(cand)
            if best is None or g > best[0]:
                best = (g, x)
    if best:
        put("FrontGainLos", fmt(best[0], 2))
        put("FrontGainHoLos", fmt(best[1][0], 2))
        put("FrontGainContourLos", best[1][5].split("_e")[1])
    so, ut = std("suburban"), fam("suburban", "ut")
    if ut:
        b = min(ut.values(), key=lambda x: x[4])
        put("FrontOutUtSub", fmt(b[4], 3))
        put("FrontHoUtSub", fmt(b[0], 2))
        put("FrontContourSub", b[5].split("_e")[1])
        cand = [o for o in so if o[0] <= b[0]]
        if cand:
            o = min(cand, key=lambda x: x[4])
            put("FrontOutStdSub", fmt(o[4], 3))
    t1 = fam("los", "t1")
    for e in [45, 50]:
        if f"t1_e{e}" in t1:
            put(f"FrontTone{ 'Fortyfive' if e == 45 else 'Fifty'}Ho", fmt(t1[f"t1_e{e}"][0], 2))
    t1c = fam("los", "t1c")
    for e in [40, 45, 50]:
        if f"t1c_e{e}" in t1c:
            put(f"FrontTonec{ {40: 'Forty', 45: 'Fortyfive', 50: 'Fifty'}[e]}Ho", fmt(t1c[f"t1c_e{e}"][0], 2))
    w = fam("los", "t1w")
    if w:
        put("FrontToneWinMinHo", fmt(min(x[0] for x in w.values()), 3))
        put("FrontToneWinMaxHo", fmt(max(x[0] for x in w.values()), 3))
    a3s = fam("suburban", "a3")
    for nm, key in [("Six", "a3_o6_t160"), ("Eight", "a3_o8_t160")]:
        if key in a3s:
            put(f"FrontAthree{nm}FastOutSub", fmt(a3s[key][4], 3))
            put(f"FrontAthree{nm}FastHoSub", fmt(a3s[key][0], 2))
    a3 = fam("los", "a3")
    if "a3_o10_t640" in a3:
        put("FrontAthreeTenOutLos", fmt(a3["a3_o10_t640"][4], 2))
    ue = fam("los", "ut")
    if "ut_e50" in ue:
        put("FrontUtFiftyHo", fmt(ue["ut_e50"][0], 2))


def fig_eph(out, outemp):
    if not out:
        return
    fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.6))
    sty = [(k, lab) + PSTY[k] for k, lab in [("mst", "MST"), ("tte", "MST + margin"),
           ("tte-under", "MST + margin, quarter budget"), ("ut", "UT-CHO"), ("t1", "T1"), ("d2", "D2"),
           ("elev", "Elevation")]]
    levels = ["sa0", "sa1000", "sa2000", "sa5000", "sa10000", "sa20000"]
    for j, (pol, lab, col, mk, ls) in enumerate(sty):
        d = out.get(pol)
        if not d:
            continue
        xs = [x for x in levels if x in d]
        off = (j - 3) * 0.06
        xk = [levels.index(x) + off for x in xs]
        ax[0].errorbar(xk, [d[x][4] for x in xs], yerr=[d[x][5] for x in xs], color=col, marker=mk, ls=ls, ms=3.4,
                       lw=1.0, capsize=1.2, label=lab, mfc=col if pol in ("tte", "ut") else "white")
        ax[1].errorbar(xk, [d[x][0] for x in xs], yerr=[d[x][1] for x in xs], color=col, marker=mk, ls=ls, ms=3.4,
                       lw=1.0, capsize=1.2, mfc=col if pol in ("tte", "ut") else "white")
    for a in ax[:2]:
        a.set_xticks(range(len(levels)))
        a.set_xticklabels(["0", "1", "2", "5", "10", "20"])
        a.set_xlabel(r"along-track error $\sigma_A$ (km), categorical")
        a.set_ylim(bottom=0)
        a.grid(True, axis="y")
    ax[0].set_ylabel("late handovers (%)")
    ax[1].set_ylabel("outage (%)")
    ax[0].set_title("(a) Gaussian error: late", fontsize=8)
    ax[1].set_title("(b) Gaussian error: outage", fontsize=8)
    bins = ["age0-6", "age6-12", "age12-24", "age24-72"]
    for j, (pol, lab, col, mk, ls) in enumerate(sty):
        d = outemp.get(pol)
        if not d:
            continue
        xs = [b for b in bins if b in d]
        off = (j - 3) * 0.06
        ax[2].errorbar([bins.index(b) + off for b in xs], [d[b][4] for b in xs], yerr=[d[b][5] for b in xs],
                       color=col, marker=mk, ls=ls, ms=3.4, lw=1.0, capsize=1.2,
                       mfc=col if pol in ("tte", "ut") else "white")
    ax[2].set_xticks(range(len(bins)))
    ax[2].set_xticklabels(["0-6", "6-12", "12-24", "24-72"])
    ax[2].set_xlabel("element-set age (h)")
    ax[2].set_ylabel("late handovers (%)")
    ax[2].set_ylim(bottom=0)
    ax[2].grid(True, axis="y")
    ax[2].set_title("(c) measured error: late", fontsize=8)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(w_pad=1.0, rect=(0, 0, 1, 0.84))
    fig.savefig(os.path.join(FIG, "fig_eph.pdf"))
    plt.close(fig)


def fig_tle(res):
    if not res:
        return
    E, bins = res
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.4))
    for k, (lab, col, ls) in enumerate([("radial", C_TEAL, "--"), ("along track", C_BLUE, "-"),
                                        ("cross track", C_GREY, "-.")]):
        x = np.sort(np.abs(E[:, k + 1]))
        ax[0].plot(x, np.arange(1, len(x) + 1) / len(x), color=col, ls=ls, lw=1.2, label=lab)
    ax[0].set_xscale("log")
    ax[0].set_xlim(1, 1e5)
    ax[0].set_xlabel("absolute error (m)")
    ax[0].set_ylabel("empirical CDF")
    ax[0].set_title("(a) error by axis, all element-set ages", fontsize=8)
    ax[0].legend(loc="upper left", frameon=False)
    ax[0].grid(True, which="major")
    ax[1].scatter(E[:, 0], np.abs(E[:, 2]) / 1000, s=1.4, color=C_LIGHT, alpha=0.45, lw=0, rasterized=True)
    for j, (lo, hi, s, n) in enumerate(bins):
        ax[1].plot([lo, hi], [s / 1000, s / 1000], color=C_INK, lw=1.6,
                   label="68.3rd percentile per age bin" if j == 0 else None)
    ax[1].legend(loc="lower right", frameon=False)
    ax[1].set_yscale("log")
    ax[1].set_xlim(0, 72)
    ax[1].set_ylim(0.01, 100)
    ax[1].set_xlabel("element-set age (h)")
    ax[1].set_ylabel("|along-track error| (km)")
    ax[1].set_title("(b) along-track error against age", fontsize=8)
    ax[1].grid(True, which="major")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(os.path.join(FIG, "fig_tle.pdf"), dpi=600)
    plt.close(fig)


def ue_se(exp, ch, key):
    v = []
    for r in summaries().get(tag_of(exp, ch, key), {}).values():
        for u in csv.DictReader(open(os.path.join(r["_dir"], "ue_summary.csv"))):
            if u["mean_se"] not in ("", "nan"):
                v.append(float(u["mean_se"]))
    return np.array(v)


UE_CDF_POL = [("a3tf", "A3-tuned + net", C_GOLD, "-"), ("t1f", "T1 + net", C_PURPLE, ":"),
              ("d2f", "D2 + net", C_GREEN, "-."), ("mst", "MST", C_BLUE, "--"),
              ("ut", "UT-CHO", C_ORANGE, "-"), ("ut40", "UT-CHO-40", C_ORANGE, "--")]


def ue_block():
    for ch in ["los", "suburban"]:
        for key in ["a3tf", "t1f", "d2f", "elevf", "mst", "mst40", "ut", "ut40"]:
            v = ue_se("main", ch, key)
            if len(v):
                put(f"Se{CHN[ch]}{POLN[key]}Pfive", fmt(float(np.percentile(v, 5)), 2))
                put(f"Se{CHN[ch]}{POLN[key]}Med", fmt(float(np.median(v)), 2))


def fig_ue_cdf():
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.4), sharey=True)
    for k, (ch, title) in enumerate([("los", "(a) line of sight"), ("suburban", "(b) suburban LOS process")]):
        for key, lab, col, ls in UE_CDF_POL:
            v = np.sort(ue_se("main", ch, key))
            if len(v) == 0:
                continue
            ax[k].plot(v, np.arange(1, len(v) + 1) / len(v), color=col, ls=ls, lw=1.2, label=lab)
        ax[k].set_xlabel("per-terminal spectral efficiency (bit/s/Hz)")
        ax[k].set_title(title, fontsize=8)
        ax[k].grid(True)
    ax[0].set_ylabel("empirical CDF")
    ax[0].legend(loc="upper left", frameon=False)
    fig.tight_layout(w_pad=1.2)
    fig.savefig(os.path.join(FIG, "fig_ue_cdf.pdf"))
    plt.close(fig)


def fig_trace():
    pols = [("a3t", "A3-tuned", C_GOLD), ("t1", "T1", C_PURPLE), ("mst", "MST", C_BLUE), ("ut", "UT-CHO", C_ORANGE),
            ("ut40", "UT-CHO-40", C_ORANGE)]
    have = [p for p in pols if os.path.exists(os.path.join(RUNS, "trace", p[0], "seed01", "trace.csv"))]
    if not have:
        return
    fig, ax = plt.subplots(len(have), 1, figsize=(3.45, 0.85 * len(have) + 0.5), sharex=True)
    ax = np.atleast_1d(ax)
    for k, (key, lab, col) in enumerate(have):
        tr = defaultdict(list)
        serv = []
        for r in csv.DictReader(open(os.path.join(RUNS, "trace", key, "seed01", "trace.csv"))):
            t, sat, el = float(r["t_s"]), int(r["sat"]), float(r["elev_deg"])
            tr[sat].append((t, el))
            if r["serving"] == "1":
                serv.append((t, el, sat))
        for sat, pts in tr.items():
            pts.sort()
            seg = [pts[0]]
            for p in pts[1:]:
                if p[0] - seg[-1][0] > 1.5:
                    ax[k].plot(*zip(*seg), color="#B9C3CC", lw=0.6)
                    seg = []
                seg.append(p)
            ax[k].plot(*zip(*seg), color="#B9C3CC", lw=0.6)
        serv.sort()
        hos = sum(1 for a, b in zip(serv, serv[1:]) if a[2] != b[2])
        put("TraceHo" + POLN[key], str(hos))
        seg = [serv[0]]
        for p in serv[1:]:
            if p[2] != seg[-1][2] or p[0] - seg[-1][0] > 1.5:
                ax[k].plot([x[0] for x in seg], [x[1] for x in seg], color=col, lw=1.6)
                seg = []
            seg.append(p)
        ax[k].plot([x[0] for x in seg], [x[1] for x in seg], color=col, lw=1.6)
        ax[k].set_ylim(22, 92)
        ax[k].set_yticks([30, 60, 90])
        ax[k].text(0.015, 0.90, f"{lab}: {hos} handovers", transform=ax[k].transAxes, va="top", fontsize=7,
                   bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))
        ax[k].grid(True, axis="y")
    ax[len(have) // 2].set_ylabel("elevation (deg)")
    ax[-1].set_xlabel("time (s)")
    ax[-1].set_xlim(0, 600)
    fig.tight_layout(h_pad=0.25)
    fig.savefig(os.path.join(FIG, "fig_trace.pdf"))
    plt.close(fig)


def cost_macros():
    for ch in ["los", "suburban"]:
        for key in ["a3t", "mst", "ut", "ut40"]:
            tag = tag_of("main", ch, key)
            if summaries().get(tag):
                put(f"Wall{CHN[ch]}{POLN[key]}", fmt(float(np.median(list(vals(tag, "wall").values()))), 0))
                put(f"Evals{CHN[ch]}{POLN[key]}", fmt(float(np.median(list(vals(tag, "evals").values()))), 2))
    for ch in ["suburban"]:
        for key in ["ut", "ut40"]:
            for tag in [f"ephblk/{ch}/{key}"]:
                if summaries().get(tag):
                    put(f"EvalsBlk{CHN[ch]}{POLN[key]}", fmt(float(np.median(list(vals(tag, "evals").values()))), 2))


def rank_macros():
    for ch in ["los", "suburban"]:
        for key in ["mst", "mst_prep", "mst_prepr", "ut"]:
            tot = top = 0
            for r in summaries().get(tag_of("main", ch, key), {}).values():
                for e in events_rows(r["_dir"]):
                    if e["kind"] == "ho" and e.get("target_rank_all", "") not in ("", "-1"):
                        tot += 1
                        top += int(float(e["target_rank_all"]) == 1)
            if tot:
                put(f"TopRank{CHN[ch]}{POLN[key]}", fmt(100.0 * top / tot, 0))


# ----------------------------------------------------------------- interference
INTERF_V = [("ref", "None", "None"), ("k3c176", "Kthree", r"$M{=}3$"), ("k1c176", "Kone", r"$M{=}1$"),
            ("k1c10", "Konehi", r"$M{=}1$, 10\,dB")]
INTERF_POL = ["a3tf", "d2f", "t1f", "elevf", "t1", "d2f4050", "mst", "mst40", "ut", "ut40"]


def itag(ch, v, pol):
    if v == "ref":
        return tag_of("main" if ch == "los" else "urban", ch, pol)
    return f"interf/{ch}/{v}/{pol}"


def interf_block():
    """Co-channel interference sensitivity: reuse 3 and 1 at 17.6 dB coupling, reuse 1 at 10 dB."""
    if not any(t.startswith("interf/") for t in summaries()):
        return False
    F = defaultdict(list)
    for ch in ["los", "urban"]:
        C = CHN[ch]
        for v, V, _ in INTERF_V:
            for pol in INTERF_POL:
                tag = itag(ch, v, pol)
                if not summaries().get(tag):
                    continue
                for mk in ["ho", "out", "snr", "rlf", "hof", "esc"]:
                    put(f"If{C}{V}{POLN[pol]}{MNAME[mk]}", fmt(stat(tag, mk)[0], DIG[mk]))
            if v == "ref":
                continue
            if ch == "los":
                for a, b, mk in [("ut", "mst", "ho"), ("mst", "t1", "ho"), ("ut40", "a3tf", "ho"),
                                 ("ut40", "a3tf", "snr"), ("mst", "a3tf", "snr")]:
                    F[f"If{C}{V}"].append((f"If{C}{V}{POLN[a]}Vs{POLN[b]}{MNAME[mk]}", itag(ch, v, a),
                                           itag(ch, v, b), mk))
            else:
                for b in ["d2f", "t1f", "elevf", "d2f4050", "mst40"]:
                    F[f"If{C}{V}"].append((f"If{C}{V}UtfortyVs{POLN[b]}Out", itag(ch, v, "ut40"),
                                           itag(ch, v, b), "out"))
        # how much interference costs each policy in SNR: SINR(K=1) - SNR(ref), seed-paired
        for pol in INTERF_POL:
            r = paired(itag(ch, "k1c176", pol), itag(ch, "ref", pol), "snr")
            if r:
                put(f"If{C}Drop{POLN[pol]}", fmt(abs(r["diff"]), 2))
            r = paired(itag(ch, "k1c10", pol), itag(ch, "ref", pol), "snr")
            if r:
                put(f"If{C}Drophi{POLN[pol]}", fmt(abs(r["diff"]), 2))
    run_families(F)
    # the same comparisons at the reference (no interference), for the table's reading
    return True


def table_interf():
    if not any(t.startswith("interf/") for t in summaries()):
        return
    rows = [("a3tf", "A3-tuned + net"), ("d2f", "D2 + net"), ("t1f", "T1 + net"), ("elevf", "Elevation + net"),
            ("t1", "T1"), ("d2f4050", "D2-40/50 + net"), ("mst", "MST"), ("mst40", "MST-40"),
            ("ut", r"\textbf{UT-CHO}"), ("ut40", r"\textbf{UT-CHO-40}")]
    hdr = " & ".join(lab for _, _, lab in INTERF_V)
    L = [r"\begin{table*}[!t]",
         r"\caption{Co-channel interference sensitivity, ten seeds. Satellites are assigned cyclically, by index, to $M$ carriers (reuse factor $M$); "
         r"every other covering satellite on the serving carrier adds its received main-lobe power less a coupling "
         r"of 17.6\,dB (first sidelobe of the TR~38.811 aperture pattern) or 10\,dB, at full load. "
         r"Only satellites above $e_{\min}$ interfere, each at its main-lobe power less a fixed coupling rather "
         r"than through its actual off-axis gain; the model is therefore a sensitivity, not a calibrated "
         r"interference budget. None: orthogonal resources (the main campaign). Quality is then an SINR. "
         r"Under line of sight T1 is shown without the fallback, as in the main table.}",
         r"\label{tab:interf}", r"\centering", r"\scriptsize", r"\setlength{\tabcolsep}{2.5pt}",
         r"\begin{tabular}{@{}l" + "rrrr" * 3 + r"@{}}", r"\toprule",
         r" & \multicolumn{4}{c}{LOS: HO/terminal-min} & \multicolumn{4}{c}{LOS: mean SINR (dB)} & "
         r"\multicolumn{4}{c}{Urban: outage (\%)} \\",
         r"\cmidrule(lr){2-5}\cmidrule(lr){6-9}\cmidrule(l){10-13}",
         "Policy & " + hdr + " & " + hdr + " & " + hdr + r" \\", r"\midrule"]
    for pol, lab in rows:
        c = []
        for ch, mk, d in [("los", "ho", 3), ("los", "snr", 2), ("urban", "out", 2)]:
            for v, _, _ in INTERF_V:
                tag = itag(ch, v, pol)
                m = stat(tag, mk)[0] if summaries().get(tag) else float("nan")
                c.append(fmt(m, d))
        L.append(lab + " & " + " & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    open(os.path.join(PAPER, "table_interf.tex"), "w").write("\n".join(L) + "\n")


def seed_macros():
    s = summaries()
    words = {5: "five", 10: "ten", 20: "twenty"}
    for exp, nm in [("main", "Main"), ("urban", "Urb"), ("eph", "Eph"), ("pareto", "Pareto")]:
        ns = [len(v) for t, v in s.items() if t.startswith(exp + "/")]
        if ns:
            n = int(np.median(ns))
            put(f"NSeed{nm}", words.get(n, str(n)))
            put(f"NSeed{nm}Num", str(n))


# ----------------------------------------------------------------- corrections and control runs
def r3_macros():
    """Macros for the corrections and the two control runs."""
    # issue 13: the inertial offset projected on the Earth-relative track, and its normal part
    rs = RE + H550
    n = math.sqrt(MU / rs ** 3)
    vg = n * RE
    lat, inc = math.radians(30.0), math.radians(53.0)
    s_az = math.cos(inc) / math.cos(lat)
    north, east = vg * math.sqrt(1 - s_az ** 2), vg * s_az
    reln, rele = north, east - OMEGA * RE * math.cos(lat)
    rel = math.hypot(reln, rele)
    proj = (north * reln + east * rele) / (vg * rel)
    put("RotKappaProj", fmt(proj * vg / rel, 3))
    put("RotNormalPct", fmt(100.0 * math.sqrt(max(0.0, 1 - proj ** 2)), 1))
    # issue 16: UT-CHO-40 and MST-40 under line of sight, seed by seed
    a, b = vals("main/los/ut40", "ho"), vals("main/los/mst40", "ho")
    ks = sorted(set(a) & set(b))
    if ks:
        d = [abs(a[k] - b[k]) for k in ks]
        put("SeedMaxLosUtfortyVsMstfortyHo", fmt(max(d), 3))
        put("SeedSameLosUtfortyVsMstfortyHo", str(sum(x < 1e-12 for x in d)))
    # issue 2: ping-pong at the 1 s minimum time of stay, urban and dense urban
    urb = [stat(tag_of("urban", "urban", p), "ppone")[0] for p in URB_ROWS if p and summaries().get(tag_of("urban", "urban", p))]
    if urb:
        put("UrbPponeMax", fmt(max(urb), 1))
    den = [stat(tag_of("urban", "denseurban", p), "ppone")[0] for p in ["elevf40", "d2f4050", "mst40", "ut40"]]
    put("DensePponeFortyMin", fmt(min(den), 1))
    put("DensePponeFortyMax", fmt(max(den), 1))
    # issue 6: conditional reconfigurations of the standardized conditions under line of sight
    rr = {p: stat(tag_of("main", "los", p), "rrc")[0] for p in ["a3", "a3t", "d2", "t1", "elev", "a3tf", "d2f",
                                                                  "t1f", "elevf"]}
    put("StdRrcLosMin", fmt(min(rr.values()), 1))
    put("StdRrcLosMax", fmt(max(rr.values()), 1))
    # issues 3 and 10: interference, SINR drop per reuse factor and urban outage change at reuse 1
    for v, V in [("k3c176", "Kthree"), ("k1c176", "Kone")]:
        dr = []
        for pol in INTERF_POL:
            r = paired(itag("los", v, pol), itag("los", "ref", pol), "snr")
            if r:
                dr.append(abs(r["diff"]))
        if dr:
            put(f"IfLosDrop{V}Min", fmt(min(dr), 2))
            put(f"IfLosDrop{V}Max", fmt(max(dr), 2))
    do = []
    for pol in INTERF_POL:
        if pol == "t1":
            continue
        r = paired(itag("urban", "k1c176", pol), itag("urban", "ref", pol), "out")
        if r:
            do.append(r["diff"])
    if do:
        put("IfUrbKoneOutChgMin", fmt(-max(do), 2))
        put("IfUrbKoneOutChgMax", fmt(-min(do), 2))
    # issue 4: real-shell control, MST + margin at the coverage edge
    F = defaultdict(list)
    for ch in ["los", "suburban"]:
        C = CHN[ch]
        for a, b in [("tte", "mst"), ("ut", "tte")]:
            if not (summaries().get(f"realtle/{ch}/{a}") and summaries().get(f"realtle/{ch}/{b}")):
                continue
            F["TleCtrl"].append((f"Tle{C}{POLN[a]}Vs{POLN[b]}Out", f"realtle/{ch}/{a}", f"realtle/{ch}/{b}", "out"))
            F["TleCtrl"].append((f"Tle{C}{POLN[a]}Vs{POLN[b]}Ho", f"realtle/{ch}/{a}", f"realtle/{ch}/{b}", "ho"))
    if F:
        run_families(F)
    for ch in ["los", "suburban"]:
        C = CHN[ch]
        for a, b in [("tte", "mst"), ("ut", "tte")]:
            A = {k: late_share(x["_dir"]) for k, x in summaries().get(f"realtle/{ch}/{a}", {}).items()}
            B = {k: late_share(x["_dir"]) for k, x in summaries().get(f"realtle/{ch}/{b}", {}).items()}
            ks = sorted(set(A) & set(B))
            if ks:
                d = [A[k] - B[k] for k in ks]
                m, h, _ = ci(d)
                nm = f"Tle{C}{POLN[a]}Vs{POLN[b]}Late"
                put(nm, fmt(m, 2))
                put(nm + "Ci", fmt(h, 2))
                put(nm + "P", pfmt(wilcoxon_p(d)))
                put(nm + "Lower", str(sum(x < -1e-12 for x in d)))
                put(nm + "Higher", str(sum(x > 1e-12 for x in d)))
        # what kind of decision is late (planned exit or escape) for MST, MST + margin and UT-CHO
        for pol in ["mst", "tte", "ut"]:
            tot = pl = 0
            for r in summaries().get(f"realtle/{ch}/{pol}", {}).values():
                for e in events_rows(r["_dir"]):
                    if e["kind"] == "ho" and e["true_serv_exit_s"] not in ("", None) and float(e["true_serv_exit_s"]) <= 0:
                        tot += 1
                        pl += e["trigger"] == "planned"
            if tot:
                put(f"Tle{C}{POLN[pol]}LateN", str(tot))
                put(f"Tle{C}{POLN[pol]}LatePlannedN", str(pl))
        gl = [abs(st.mean(late_share(x["_dir"]) for x in summaries()[f"realtle/{ch}/ut"].values())
              - st.mean(late_share(x["_dir"]) for x in summaries()[f"realtle/{ch}/tte"].values()))
          for ch in ["los", "suburban"] if summaries().get(f"realtle/{ch}/tte")]
    if gl:
        put("TleLosUtVsTteLateAbsMax", fmt(max(gl), 2))
    # issue 11: UT-CHO at the coverage edge under blockage with measured error, against MST
    for ch in ["suburban", "urban"]:
        C = CHN[ch]
        if not summaries().get(f"ephblk/{ch}/ut"):
            continue
        r = paired(f"ephblk/{ch}/ut", f"ephblk/{ch}/mst", "out")
        if r:
            put(f"Blk{C}UtVsMstOut", fmt(r["diff"], 3))
            put(f"Blk{C}UtVsMstOutCi", fmt(r["ci"], 3))
            put(f"Blk{C}UtVsMstOutP", pfmt(r["p"]))
        A = {k: late_share(x["_dir"]) for k, x in summaries().get(f"ephblk/{ch}/ut", {}).items()}
        B = {k: late_share(x["_dir"]) for k, x in summaries().get(f"ephblk/{ch}/mst", {}).items()}
        ks = sorted(set(A) & set(B))
        if ks:
            d = [A[k] - B[k] for k in ks]
            put(f"Blk{C}UtVsMstLateP", pfmt(wilcoxon_p(d)))
            put(f"Blk{C}UtVsMstLateLower", str(sum(x < -1e-12 for x in d)))
    # issue 1: the largest gap between a late share with error and with an exact ephemeris, under
    # blockage, over every policy but those executing at the coverage edge
    gaps = []
    for ch in ["suburban", "urban"]:
        exp0 = "main" if ch == "suburban" else "urban"
        for pol in ["a3tf", "d2f", "t1f", "elevf", "mst40", "tte40", "ut40"]:
            tag = f"ephblk/{ch}/{pol}"
            twin = tag_of(exp0, ch, {"tte40": "mst40"}.get(pol, pol))
            if summaries().get(tag) and summaries().get(twin):
                la = st.mean(late_share(r["_dir"]) for r in summaries()[tag].values())
                lb = st.mean(late_share(r["_dir"]) for r in summaries()[twin].values())
                gaps.append(abs(la - lb))
    if gaps:
        put("BlkLateExactMaxGap", fmt(max(gaps), 3))


# ----------------------------------------------------------------- age-aware margin
# Specified in advance: families, metrics and
# the non-inferiority margin below are fixed there.
AM_DELTA = 0.5  # non-inferiority margin on the late share (percentage points)
AM_POL = {"mst": "Mst", "tte": "Fix", "ttefz": "Fz", "tteage": "Age", "ut": "Ut", "utage": "Utage"}
AM_SET = {"AmLos": "realtle/los", "AmSub": "realtle/suburban", "AmGau": "ephage/los"}
AM_TABLE = os.path.join(HERE, "tle", "age_sigma_table.csv")


def am_seed_metrics(run_dir):
    late = ho = 0
    aband = []
    margins = []
    for e in events_rows(run_dir):
        if e["kind"] != "ho":
            continue
        ho += 1
        x = e["true_serv_exit_s"]
        if x in ("", None):
            continue
        margins.append(float(e["margin_s"]))
        if float(x) <= 0.0:
            late += 1
        if e["trigger"] == "planned":
            aband.append(max(0.0, float(x)))
    return {"late": 100.0 * late / max(1, ho), "aband": st.mean(aband) if aband else float("nan"),
            "_margins": margins}


_AM = {}


def am_vals(tag, mk):
    if tag not in _AM:
        _AM[tag] = {k: am_seed_metrics(r["_dir"]) for k, r in summaries().get(tag, {}).items()}
    if mk in ("late", "aband"):
        return {k: v[mk] for k, v in _AM[tag].items()}
    return vals(tag, mk)


def am_ni_p(d, delta):
    """One-sided paired t test of H0: mean(d) >= delta against mean(d) < delta."""
    n = len(d)
    if st.stdev(d) == 0.0:
        return 0.0 if st.mean(d) < delta else 1.0
    t = (st.mean(d) - delta) / (st.stdev(d) / math.sqrt(n))
    try:
        from scipy import stats as _ss
        return float(_ss.t.cdf(t, n - 1))
    except ImportError:  # pragma: no cover
        raise


def am_family(name, tests):
    """tests: (macro, tagA, tagB, metric, kind), kind 'ni' (late, delta) or 'sup' (two-sided Wilcoxon)."""
    res = []
    for nm, a, b, mk, kind in tests:
        A, B = am_vals(a, mk), am_vals(b, mk)
        ks = sorted(set(A) & set(B))
        if len(ks) < 2:
            continue
        d = [A[k] - B[k] for k in ks]
        m, h, n = ci(d)
        p = am_ni_p(d, AM_DELTA) if kind == "ni" else wilcoxon_p(d)
        res.append((nm, kind, m, h, n, p, sum(x < -1e-12 for x in d), sum(x > 1e-12 for x in d), wilcoxon_p(d)))
    if not res:
        return {}
    adj = holm([r[5] for r in res])
    out = {}
    put(name + "Size", str(len(res)))
    nrej = 0
    for (nm, kind, m, h, n, p, lo, hi, pw), pa in zip(res, adj):
        dg = 2 if nm.endswith("Late") else 3
        put(nm, fmt(m, dg))
        put(nm + "Abs", fmt(abs(m), dg))
        put(nm + "Ci", fmt(h, dg))
        put(nm + "Up", fmt(m + h, dg))
        put(nm + "P", pfmt(p))
        put(nm + "H", pfmt(pa))
        put(nm + "W", pfmt(pw))
        put(nm + "Lower", str(lo))
        put(nm + "Higher", str(hi))
        put(nm + "S", "shown" if pa < ALPHA else "not shown")
        nrej += pa < ALPHA
        out[nm] = {"kind": kind, "diff": m, "ci": h, "n": n, "p": p, "holm": pa, "lower": lo, "higher": hi}
    put(name + "Rej", str(nrej))
    return out


def agemargin_block():
    """The pre-specified age-aware margin experiment."""
    have = [s for s, pre in AM_SET.items() if summaries().get(pre + "/tteage")]
    if not have:
        return {}
    res = {"means": {}}
    for s, pre in AM_SET.items():
        for pol, P in AM_POL.items():
            tag = f"{pre}/{pol}"
            if not summaries().get(tag):
                continue
            row = {}
            for mk in ["late", "aband", "ho", "out", "hof"]:
                row[mk] = ci(list(am_vals(tag, mk).values()))
                put(f"{s}{P}{MNAME.get(mk, mk.capitalize())}", fmt(row[mk][0], 2 if mk in ("late", "aband", "hof") else 3))
            ms = [x for v in _AM.get(tag, {}).values() for x in v["_margins"]]
            if ms:
                q = np.percentile(ms, [5, 50, 95])
                row["margin"] = q
                put(f"{s}{P}MarginLo", fmt(q[0], 2))
                put(f"{s}{P}MarginMed", fmt(q[1], 2))
                put(f"{s}{P}MarginHi", fmt(q[2], 2))
                put(f"{s}{P}MarginMean", fmt(float(np.mean(ms)), 2))
            res["means"][tag] = row
    R = "realtle/{}/{}"
    G = "ephage/los/{}"
    res["P"] = am_family("AmFamP", [
        ("AmLosAgeVsFixLate", R.format("los", "tteage"), R.format("los", "tte"), "late", "ni"),
        ("AmLosAgeVsFixAband", R.format("los", "tteage"), R.format("los", "tte"), "aband", "sup"),
        ("AmSubAgeVsFixLate", R.format("suburban", "tteage"), R.format("suburban", "tte"), "late", "ni"),
        ("AmSubAgeVsFixAband", R.format("suburban", "tteage"), R.format("suburban", "tte"), "aband", "sup")])
    S1, S2, S3 = [], [], []
    for ch, C in [("los", "AmLos"), ("suburban", "AmSub")]:
        S1 += [(f"{C}UtageVsUtLate", R.format(ch, "utage"), R.format(ch, "ut"), "late", "ni"),
               (f"{C}UtageVsUtAband", R.format(ch, "utage"), R.format(ch, "ut"), "aband", "sup"),
               (f"{C}UtageVsUtHo", R.format(ch, "utage"), R.format(ch, "ut"), "ho", "sup"),
               (f"{C}UtageVsUtOut", R.format(ch, "utage"), R.format(ch, "ut"), "out", "sup")]
        S2 += [(f"{C}AgeVsFixHo", R.format(ch, "tteage"), R.format(ch, "tte"), "ho", "sup"),
               (f"{C}AgeVsFixOut", R.format(ch, "tteage"), R.format(ch, "tte"), "out", "sup"),
               (f"{C}AgeVsFixHof", R.format(ch, "tteage"), R.format(ch, "tte"), "hof", "sup")]
        S3 += [(f"{C}AgeVsFzLate", R.format(ch, "tteage"), R.format(ch, "ttefz"), "late", "sup"),
               (f"{C}AgeVsFzAband", R.format(ch, "tteage"), R.format(ch, "ttefz"), "aband", "sup"),
               (f"{C}FzVsFixLate", R.format(ch, "ttefz"), R.format(ch, "tte"), "late", "sup"),
               (f"{C}FzVsFixAband", R.format(ch, "ttefz"), R.format(ch, "tte"), "aband", "sup")]
    res["S1"] = am_family("AmFamSone", S1)
    res["S2"] = am_family("AmFamStwo", S2)
    res["S3"] = am_family("AmFamSthree", S3)
    res["G"] = am_family("AmFamG", [
        ("AmGauAgeVsFixLate", G.format("tteage"), G.format("tte"), "late", "ni"),
        ("AmGauAgeVsFixAband", G.format("tteage"), G.format("tte"), "aband", "sup"),
        ("AmGauUtageVsUtLate", G.format("utage"), G.format("ut"), "late", "ni"),
        ("AmGauUtageVsUtAband", G.format("utage"), G.format("ut"), "aband", "sup")])
    res["G2"] = am_family("AmFamGtwo", [
        ("AmGauAgeVsFzLate", G.format("tteage"), G.format("ttefz"), "late", "sup"),
        ("AmGauAgeVsFzAband", G.format("tteage"), G.format("ttefz"), "aband", "sup")])
    P = res["P"]
    ok = len(P) == 4 and all(r["holm"] < ALPHA for r in P.values()) and \
        all(r["diff"] < 0 for k, r in P.items() if r["kind"] == "sup")
    put("AmVerdict", "supported" if ok else "not supported")
    res["verdict"] = ok
    # late share by the age bin of the serving satellite's element set (real shell, pooled seeds)
    edges = [float(r["age_lo_h"]) for r in csv.DictReader(open(AM_TABLE))]
    for ch, C in [("los", "AmLos"), ("suburban", "AmSub")]:
        for pol in ["tte", "tteage"]:
            young = [0, 0]
            old = [0, 0]  # 12 h and older
            tot_n = [0]
            for r in summaries().get(R.format(ch, pol), {}).values():
                af = os.path.join(r["_dir"], "sat_ages.csv")
                ages = None
                if os.path.exists(af):
                    ages = {int(a["sat"]): float(a["age0_h"]) for a in csv.DictReader(open(af))}
                else:
                    tw = os.path.join(RUNS, R.format(ch, "tteage"), os.path.basename(r["_dir"]), "sat_ages.csv")
                    if os.path.exists(tw):
                        ages = {int(a["sat"]): float(a["age0_h"]) for a in csv.DictReader(open(tw))}
                if ages is None:
                    continue
                for e in events_rows(r["_dir"]):
                    if e["kind"] != "ho" or e["true_serv_exit_s"] in ("", None):
                        continue
                    ag = ages[int(e["from_sat"])]
                    tot_n[0] += 1
                    if 6.0 <= ag < 12.0:
                        continue
                    b = young if ag < 6.0 else old
                    b[0] += 1
                    b[1] += float(e["true_serv_exit_s"]) <= 0.0
            P2 = AM_POL[pol]
            if young[0]:
                put(f"{C}{P2}LateYoung", fmt(100.0 * young[1] / young[0], 2))
                put(f"{C}{P2}ShareYoung", fmt(100.0 * young[0] / tot_n[0], 0))
            if old[0]:
                put(f"{C}{P2}LateOld", fmt(100.0 * old[1] / old[0], 2))
    cj = json.load(open(os.path.join(RUNS, "realtle/los/tte/seed01/cmd.json")))["cmd"]
    put("TleBudgetA", fmt([float(a.split("=")[1]) for a in cj if a.startswith("--assumedSigmaAlongM=")][0] / 1000, 2))
    put("AmZ", "2.326")
    put("AmDelta", fmt(AM_DELTA, 1))
    put("AmBins", str(len(edges)))
    return res


def table_agemargin(res):
    """Supplement Table S6: the pre-specified age-aware margin experiment."""
    if not res:
        return
    M = res["means"]
    lab = [("mst", "MST (no margin)"), ("tte", "MST + margin (fixed)"), ("ttefz", "MST + margin, all-age, $z$=2.33"),
           ("tteage", "MST + age margin"), ("ut", "UT-CHO (fixed)"), ("utage", "UT-CHO + age margin")]
    L = [r"\begin{table*}[!t]",
         r"\caption{Pre-registered age-aware margin, ten seeds. Late: share of handovers decided after the true exit; "
         r"Ab.: serving time given up per planned handover (s); HO per terminal-minute; Out.\ in \%; "
         r"$m$: 5th to 95th percentile of the margin (s). Gaussian: each satellite's age bin drawn with the measured "
         r"frequencies, its copy offset by a Gaussian error with that bin's sigmas. The fixed rows budget one sigma for "
         r"every satellite, $z=2$: along track \TleBudgetA\,km on the real shell (as in the main text) and the "
         r"measured all-age \TleSigA\,km in the Gaussian setting; the other margin rows use $z=\AmZ$.}",
         r"\label{tab:agemargin}", r"\centering", r"\footnotesize", r"\setlength{\tabcolsep}{2.2pt}",
         r"\begin{tabular}{@{}lrrrrrrrrrrrrr@{}}", r"\toprule",
         r" & \multicolumn{5}{c}{Real shell, LOS} & \multicolumn{4}{c}{Real shell, suburban} & "
         r"\multicolumn{4}{c}{Gaussian with age, LOS} \\",
         r"\cmidrule(lr){2-6}\cmidrule(lr){7-10}\cmidrule(l){11-14}",
         r"Policy & Late & Ab. & HO & Out. & $m$ & Late & Ab. & HO & Out. & Late & Ab. & HO & $m$ \\", r"\midrule"]
    for pol, name in lab:
        row = []
        for s, pre in AM_SET.items():
            r = M.get(f"{pre}/{pol}")
            if not r:
                row += ["--"] * (4 if s != "AmLos" else 5)
                continue
            mm = r.get("margin")
            mr = (fmt(mm[0], 2) + "--" + fmt(mm[2], 2)) if mm is not None else "--"
            if s == "AmLos":
                row += [fmt(r["late"][0], 2), fmt(r["aband"][0], 2), fmt(r["ho"][0], 3), fmt(r["out"][0], 3), mr]
            elif s == "AmSub":
                row += [fmt(r["late"][0], 2), fmt(r["aband"][0], 2), fmt(r["ho"][0], 3), fmt(r["out"][0], 3)]
            else:
                row += [fmt(r["late"][0], 2), fmt(r["aband"][0], 2), fmt(r["ho"][0], 3), mr]
        L.append(name + " & " + " & ".join(row) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    open(os.path.join(PAPER, "table_agemargin.tex"), "w").write("\n".join(L) + "\n")


def fit_width(name, width):
    """Scale a generated table to the text or column width (it is wider at its natural size)."""
    p = os.path.join(PAPER, name)
    t = open(p).read()
    if "\\resizebox" not in t:
        i = t.index("\\begin{tabular}")
        j = t.index("\\end{tabular}") + len("\\end{tabular}")
        t = t[:i] + "\\resizebox{" + width + "}{!}{" + t[i:j] + "}" + t[j:]
        open(p, "w").write(t)


# Short TAES-style captions; the detail moves to the note under each table.
CAPTIONS = {
    "table_main.tex": ("Reference Comparison With an Exact Ephemeris",
        r"Ten seeds; means, with the 95\,\% CI for HO and outage. HO: completed handovers per terminal-minute; "
        r"PP$_1$: ping-pong, \% of handovers returning to the previous cell after a stay below 1\,s (TR 36.839); "
        r"Out.: outage; SNR in dB; ToS: median time of stay (s); RRC: conditional reconfigurations per "
        r"terminal-minute; RLF and HOF per terminal-hour; Esc.: share of handovers taken by the safety net. "
        r"+ net: the safety net that every exit-triggered rule also carries. D2 and T1 carry an A4 leg on the "
        r"same candidate (threshold $Q$); A3 needs none."),
    "table_ctrl.tex": ("Single-Factor Controls",
        r"Ten seeds, means; each row changes one element of the policy it names. HO per terminal-minute, RLF per "
        r"terminal-hour, SNR in dB, PP$_1$ in \%. TTE-D2: D2 execution with the admitted candidates ranked by "
        r"predicted time to exit, after~\cite{noreen2026location}. $r$: two prepared slots reserved for the "
        r"strongest satellites."),
    "table_urban.tex": ("Urban and Dense-Urban Environments",
        r"TR 38.811 LOS probability, path loss and shadowing; ten seeds, means. HO per terminal-minute; PP$_1$, "
        r"outage and escape share in \%; HOF and RLF per terminal-hour. $\Delta$Out.: outage of UT-CHO-40 minus "
        r"that of the row, seed-paired mean in percentage points with its 95\,\% CI; $^{\ast}$: Holm-adjusted "
        r"$p<0.05$ in the family of the column. Aircraft and ships see the sky in LOS and trains the rural table."),
    "table_urban_ctrl.tex": ("Urban and Dense-Urban Controls and Ablations",
        r"Columns as in the urban table of the main text (ten seeds, means). $\Delta$Out.: for the hold-off rows, "
        r"outage of UT-CHO-40 minus that of the row; for the ablation rows, the ablation minus UT-CHO-40; "
        r"seed-paired mean in percentage points with its 95\,\% CI; $^{\ast}$: Holm-adjusted $p<0.05$."),
    "table_headline.tex": ("Headline Paired Comparisons",
        r"Seed-paired mean difference (first minus second) with its 95\,\% CI, the share of seeds in which the "
        r"first is lower, and the exact Wilcoxon $p$ before and after the Holm correction of its family. HO per "
        r"terminal-minute, outage in percentage points, HOF per terminal-hour."),
    "table_class.tex": ("Urban Environment by Terminal Class",
        r"Handovers and outage (s) per terminal and 600\,s run, pooled over seeds. Aircraft and ships see the sky "
        r"in LOS and trains the rural table."),
    "table_ephstd.tex": ("Ephemeris Error Against Every Policy",
        r"Ten seeds. Left: each satellite's error drawn from the measured 0 to 6\,h element-set errors; "
        r"$\Delta$Out.: outage minus that of the same policy with an exact ephemeris, seed-paired (percentage "
        r"points). Right: the real 53$^\circ$ Starlink shell, truth from the supplemental element sets and the "
        r"network's copy from the public ones, both SGP4. Late: handovers decided after the serving satellite left "
        r"the coverage edge on its true orbit. The D2 references, the T1 epoch and the elevation that Elevation "
        r"evaluates come from the network's copy; A3 uses none of them. Under blockage every late share except "
        r"those of MST and UT-CHO equals that of the same policy with an exact ephemeris to the printed precision; "
        r"these are mostly escapes after the exit. MST + margin and MST + age margin ran on the real shell only (--)."),
    "table_sens.tex": ("Optimality Gap and Sensitivity",
        r"Ten seeds; gap with its 95\,\% CI. Bound: offline greedy cover on the true orbits, static terminals, "
        r"600\,s; handovers counted as completed plus failed."),
    "table_losmodel.tex": (r"Urban Outage (\%) Under Other LOS Processes",
        r"Ten seeds. Reference: redraw after 50\,m or a new 10$^\circ$ bin, aircraft and ships in LOS. All in "
        r"clutter: every class in the urban table. Markov: first-order process with the given correlation "
        r"distance and an angular scale of 10$^\circ$ unless stated."),
    "table_interf.tex": ("Co-Channel Interference Sensitivity",
        r"Ten seeds. Satellites are assigned cyclically, by index, to $M$ carriers (reuse factor $M$); every other "
        r"covering satellite above $e_{\min}$ on the serving carrier adds its main-lobe power less a coupling of "
        r"17.6\,dB (first sidelobe of the TR~38.811 aperture pattern) or 10\,dB, at full load, rather than its "
        r"actual off-axis gain, so the model is a sensitivity, not a calibrated interference budget. None: "
        r"orthogonal resources (the main campaign); otherwise quality is an SINR. Under LOS, T1 is shown without "
        r"the safety net, as in the main table."),
    "table_agemargin.tex": ("Pre-Registered Age-Aware Margin", None),
}


def retitle(fname, title, note):
    path = os.path.join(PAPER, fname)
    t = open(path).read()
    i = t.index(r"\caption{") + len(r"\caption{")
    depth, j = 1, i
    while depth:
        depth += {"{": 1, "}": -1}.get(t[j], 0)
        j += 1
    old = t[i:j - 1]
    t = t[:i] + title + t[j - 1:]
    if note is None:  # keep the old caption text as the note
        note = old
    lines = [ln for ln in t.split("\n") if not ln.startswith(r"\\[2pt]{\scriptsize")]
    k = max(n for n, ln in enumerate(lines) if ln.startswith(r"\end{table"))
    lines.insert(k, r"\\[2pt]{\scriptsize " + note + "}")
    open(path, "w").write("\n".join(lines))



# ----------------------------------------------------------------- Proposition 1 (Oct 2026)
ZS = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 6]
MU_E = 3.986004418e14
R_E = 6378137.0


def ztag(ch, z):
    return f"realtle/{ch}/tte" if z == 2 else f"zsweep/{ch}/z{z}"


def _fo_rows(ch, pol):
    """Planned exits of the real-shell runs: observed exit error, along-track first-order error
    and the per-pass sigma (from the margin). Returns (eps_obs, eps_fo, sigma, run_key)."""
    eta = rot_kappa()
    out = []
    for k, r in summaries().get(f"realtle/{ch}/{pol}", {}).items():
        d = r["_dir"]
        f = os.path.join(d, "eph_offsets.csv")
        if not os.path.exists(f):
            continue
        cmd = json.load(open(os.path.join(d, "cmd.json")))["cmd"]
        alt = float(next((a.split("=")[1] for a in cmd if a.startswith("--altitudeKm=")), "550")) * 1e3
        z = float(next((a.split("=")[1] for a in cmd if a.startswith("--marginZ=")), "2"))
        v = math.sqrt(MU_E / (R_E + alt))
        T = num(r, "sim_time_s")
        off = {int(o["sat"]): o for o in csv.DictReader(open(f))}
        for e in events_rows(d):
            if e["kind"] != "ho" or e.get("trigger") != "planned" or e["true_serv_exit_s"] in ("", None):
                continue
            o = off.get(int(e["from_sat"]))
            if o is None:
                continue
            t = float(e["t_s"])
            da = float(o["d0_a_m"]) + (float(o["dT_a_m"]) - float(o["d0_a_m"])) * t / T
            sig = (float(e["margin_s"]) - 0.2) / z if z > 0 else float("nan")
            out.append((float(e["pred_serv_tte_s"]) - float(e["true_serv_exit_s"]), -eta * da / v, sig, k))
    return out


def _run_counts(run_dir):
    ev = [e for e in events_rows(run_dir) if e["kind"] == "ho"]
    late = sum(1 for e in ev if e["true_serv_exit_s"] not in ("", None) and float(e["true_serv_exit_s"]) <= 0.0)
    planned = sum(1 for e in ev if e.get("trigger") == "planned")
    intr = [float(e["interruption_ms"]) / 1e3 for e in ev if e.get("interruption_ms") not in ("", None)]
    return len(ev), late, planned, (st.mean(intr) if intr else float("nan"))


def prop_block():
    res = {}
    rows = _fo_rows("los", "mst") + _fo_rows("suburban", "mst")
    if len(rows) > 10:
        x = np.array([r[1] for r in rows])
        y = np.array([r[0] for r in rows])
        b, a = np.polyfit(x, y, 1)
        res["fo"] = (x, y, b, float(np.corrcoef(x, y)[0, 1]))
        put("PropSlope", fmt(b, 2))
        put("PropCorr", fmt(float(np.corrcoef(x, y)[0, 1]), 3))
        put("PropN", str(len(rows)))
        put("PropRmse", fmt(float(np.sqrt(np.mean((y - x) ** 2))), 2))
    for ch in ["los", "suburban"]:
        zr = {}
        for z in ZS:
            S = summaries().get(ztag(ch, z))
            if not S:
                continue
            late = ci([late_share(r["_dir"]) for r in S.values()])
            zr[z] = {"late": late, "out": stat(ztag(ch, z), "out"), "ho": stat(ztag(ch, z), "ho")}
        if len(zr) < 4:
            continue
        # cost fit over every run: outage seconds per terminal = a + c_h N_ho + c_l N_late (per terminal)
        X, Y, ints, sig_all = [], [], [], []
        for z in zr:
            for r in summaries()[ztag(ch, z)].values():
                n, late, planned, intr = _run_counts(r["_dir"])
                U = num(r, "num_ues")
                X.append([1.0, n / U, late / U])
                Y.append(num(r, "outage_frac") * num(r, "sim_time_s"))
                ints.append(intr)
        # c_h is the measured mean interruption of a handover; a and c_l are fitted with c_h held there
        c_h = float(np.nanmean(ints))
        X, Y = np.array(X), np.array(Y)
        coef, *_ = np.linalg.lstsq(X[:, [0, 2]], Y - c_h * X[:, 1], rcond=None)
        coef = [float(coef[0]), c_h, float(coef[1])]
        c_l = float(coef[2])
        Tsim = float(np.mean([num(r, "sim_time_s") for z in zr for r in summaries()[ztag(ch, z)].values()]))
        sig = [r[2] for r in _fo_rows(ch, "tte") if not math.isnan(r[2])]
        sbar = float(np.median(sig)) if sig else float("nan")
        zz = np.array(sorted(zr))
        H = np.array([zr[z]["ho"][0] for z in zz])
        Lz = np.array([zr[z]["late"][0] / 100.0 for z in zz])
        # H(z) = H0 (1 + k sigma z)(1 + beta l(z)) ~ H0 + H0 k sigma z + H0 beta l(z): a late decision
        # adds beta handovers (the escape and the re-entry it causes)
        A_ = np.column_stack([np.ones_like(zz, dtype=float), zz * sbar, Lz])
        cf, *_ = np.linalg.lstsq(A_, H, rcond=None)
        H0, k, beta = float(cf[0]), float(cf[1] / cf[0]), float(cf[2] / cf[0])
        c_eff = c_l + beta * c_h
        R = c_eff / (math.sqrt(2 * math.pi) * k * sbar * c_h) if k > 0 else float("nan")
        zstar = math.sqrt(2 * math.log(R)) if R > 1 else float("nan")
        # exact root of the stationarity condition
        from math import erfc, exp, pi, sqrt
        Qf = lambda t: 0.5 * erfc(t / sqrt(2))
        ph = lambda t: exp(-t * t / 2) / sqrt(2 * pi)
        hfun = lambda t: k * sbar * (c_h + c_eff * Qf(t)) - c_eff * ph(t) * (1 + k * sbar * t)
        lo, hi = 0.0, 12.0
        zex = float("nan")
        if k > 0 and hfun(lo) < 0 < hfun(hi):
            for _ in range(80):
                mid = (lo + hi) / 2
                lo, hi = (mid, hi) if hfun(mid) < 0 else (lo, mid)
            zex = (lo + hi) / 2
        # empirical-tail analysis from the measured offsets (planned exits of the z = 2 runs)
        fo = _fo_rows(ch, "tte")
        n0, l0, p0, _ = zip(*[_run_counts(r["_dir"]) for r in summaries()[ztag(ch, 0)].values()])
        pshare = sum(p0) / max(1, sum(n0))
        zfine = np.linspace(0, 6.5, 66)
        emp = [100.0 * pshare * np.mean([1.0 if e > t * sg - 0.05 else 0.0 for _, e, sg, _ in fo if not math.isnan(sg)])
               for t in zfine] if fo else []
        gau = [100.0 * pshare * Qf(t) for t in zfine]
        # escapes after a link loss are late whatever the margin: measured floor at the largest margin
        zmax = max(zr)
        floor = max(0.0, zr[zmax]["late"][0] - (float(np.interp(zmax, zfine, emp)) if emp else 0.0))
        emp = [e + floor for e in emp]
        gau = [g + floor for g in gau]
        a0 = float(coef[0])
        # outage share = a/T + H(z)/60 (c_h + c_l L(z)), H(z) = H0 (1 + k sigma z + beta L(z))
        J = lambda P: [100 * (a0 / Tsim + H0 * (1 + k * sbar * t + beta * p / 100.0) / 60.0 * (c_h + c_l * p / 100.0))
                       for t, p in zip(zfine, P)]
        zemp = float(zfine[int(np.argmin(J(emp)))]) if emp else float("nan")
        zsim = min(zr, key=lambda z: zr[z]["out"][0])
        res[ch] = dict(zr=zr, c_h=c_h, c_l=c_l, k=k, beta=beta, c_eff=c_eff, sbar=sbar, zstar=zstar, zex=zex, zemp=zemp, zsim=zsim,
                       zfine=zfine, emp=emp, gau=gau, Jemp=J(emp) if emp else [], Jgau=J(gau), pshare=pshare,
                       floor=floor)
        C = CHN[ch]
        put(f"Prop{C}Ch", fmt(c_h * 1e3, 0))
        put(f"Prop{C}Cl", fmt(c_l, 2))
        put(f"Prop{C}K", fmt(k * 1e3, 2))
        put(f"Prop{C}Beta", fmt(beta, 2))
        put(f"Prop{C}Ceff", fmt(c_eff, 2))
        put(f"Prop{C}Sigma", fmt(sbar, 2))
        put(f"Prop{C}Zstar", fmt(zstar, 2))
        put(f"Prop{C}Zexact", fmt(zex, 2))
        put(f"Prop{C}Zemp", fmt(zemp, 1))
        put(f"Prop{C}Zsim", fmt(zsim, 1))
        put(f"Prop{C}Floor", fmt(floor, 2))
        offs = [e for _, e, sg, _ in fo if not math.isnan(sg)]
        if offs:
            put(f"Prop{C}LateSide", fmt(100.0 * sum(1 for e in offs if e > 0) / len(offs), 0))
        if emp:
            Je = res[ch]["Jemp"]
            dl = [abs(float(np.interp(z, zfine, emp)) - zr[z]["late"][0]) for z in zr]
            dj = [abs(float(np.interp(z, zfine, Je)) - zr[z]["out"][0]) for z in zr]
            inci = sum(1 for z in zr if abs(float(np.interp(z, zfine, Je)) - zr[z]["out"][0]) <= zr[z]["out"][1])
            put(f"Prop{C}LateErr", fmt(max(dl), 1))
            put(f"Prop{C}OutErr", fmt(max(dj), 3))
            put(f"Prop{C}OutInCi", str(inci))
            put(f"Prop{C}NZ", str(len(zr)))
            put(f"Prop{C}A", fmt(a0, 2))
        for z in zr:
            zn = {0: "Zero", 0.5: "Half", 1: "One", 1.5: "OneHalf", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
                  6: "Six"}[z]
            put(f"Zs{C}{zn}Late", fmt(zr[z]["late"][0], 2))
            put(f"Zs{C}{zn}Out", fmt(zr[z]["out"][0], 3))
            put(f"Zs{C}{zn}OutCi", fmt(zr[z]["out"][1], 3))
            put(f"Zs{C}{zn}Ho", fmt(zr[z]["ho"][0], 3))
        r2 = paired(ztag(ch, zsim), ztag(ch, 2), "out") if zsim != 2 else None
        if r2:
            put(f"Zs{C}OptVsTwoOut", fmt(r2["diff"], 3))
            put(f"Zs{C}OptVsTwoOutCi", fmt(r2["ci"], 3))
            put(f"Zs{C}OptVsTwoOutP", pfmt(r2["p"]))
    if "los" in res:
        put("PropZstar", fmt(res["los"]["zstar"], 1))
    return res


def fig_margin(res):
    if "fo" not in res or "los" not in res:
        return
    fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.45))
    x, y, b, rr = res["fo"]
    ax[0].plot(x, y, "o", ms=1.6, color="#2a78d6", alpha=0.35, mew=0)
    lim = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax[0].plot(lim, lim, color="#1f2933", lw=0.9, ls="--", label="equality")
    ax[0].set_xlabel(r"first-order exit error $-\eta\delta_A/v$ (s)")
    ax[0].set_ylabel(r"measured exit error $\hat{\tau}-\tau$ (s)")
    ax[0].text(0.04, 0.95, f"slope {b:.2f}, r = {rr:.3f}\nn = {len(x)}", transform=ax[0].transAxes, va="top",
               fontsize=7)
    ax[0].set_title("(a) real shell, MST", fontsize=8)
    R = res["los"]
    zr = R["zr"]
    zz = sorted(zr)
    ax[1].errorbar(zz, [zr[z]["late"][0] for z in zz], yerr=[zr[z]["late"][1] for z in zz], fmt="o", ms=3.4,
                   color="#eb6834", mfc="white", capsize=1.5, lw=0.9, label="simulation")
    ax[1].plot(R["zfine"], R["emp"], color="#2a78d6", lw=1.1, label="analysis, measured offsets")
    ax[1].plot(R["zfine"], R["gau"], color="#7b8794", lw=1.0, ls="--", label=r"analysis, Gaussian $Q(z)$")
    ax[1].set_xlabel(r"margin multiplier $z$")
    ax[1].set_ylabel("late handovers (%)")
    ax[1].set_title("(b) late decisions", fontsize=8)
    ax[2].errorbar(zz, [zr[z]["out"][0] for z in zz], yerr=[zr[z]["out"][1] for z in zz], fmt="o", ms=3.4,
                   color="#eb6834", mfc="white", capsize=1.5, lw=0.9)
    ax[2].plot(R["zfine"], R["Jemp"], color="#2a78d6", lw=1.1)
    ax[2].plot(R["zfine"], R["Jgau"], color="#7b8794", lw=1.0, ls="--")
    top = ax[2].get_ylim()[1]
    for zv, c, lab, ha in [(R["zstar"], "#1f2933", r"$z^\star$", "left"), (R["zemp"], "#2a78d6", r"$z_{\mathrm{emp}}$", "right")]:
        if not math.isnan(zv):
            ax[2].axvline(zv, color=c, ls=":", lw=0.9)
            ax[2].text(zv + (0.12 if ha == "left" else -0.12), top * 0.97, lab, color=c, fontsize=7, ha=ha, va="top")
    ax[2].set_xlabel(r"margin multiplier $z$")
    ax[2].set_ylabel("outage (%)")
    ax[2].set_title("(c) outage, line of sight", fontsize=8)
    for a in ax:
        a.grid(True, alpha=0.3)
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.62, 1.0))
    fig.tight_layout(w_pad=1.0, rect=(0, 0, 1, 0.88))
    fig.savefig(os.path.join(FIG, "fig_margin.pdf"))
    plt.close(fig)

if __name__ == "__main__":
    seed_macros()
    base_macros()
    run_families(fam_main())
    run_families(fam_urban())
    put("Alpha", fmt(ALPHA, 2))
    derived()
    table_main()
    table_ctrl()
    table_urban()
    table_urban(URB_CTRL_ROWS, "table_urban_ctrl.tex", "tab:urbanctrl")
    table_headline()
    trigger_macros()
    hof_mechanism()
    pc = per_class()
    table_class(pc)
    res = acc_block()
    fig_accuracy(res)
    pp = pareto_block()
    fig_pareto(pp)
    pareto_macros(pp)
    eo = eph_block("eph")
    eph_macros(eo)
    ee = eph_block("ephemp")
    emp_macros(ee)
    fig_eph(eo, ee)
    blk = ephblk_block()
    tle = realtle_block()
    recent_block()
    table_recent()
    fig_margin(prop_block())
    table_eph_std(blk, tle)
    lbo = lb_block()
    so = sens_block()
    table_lb_sens(lbo, so)
    density_macros()
    lm = losmodel_block()
    table_losmodel(lm)
    reest_block()
    interf_block()
    table_interf()
    fig_trace()
    fig_tle(tle_block())
    ue_block()
    fig_ue_cdf()
    cost_macros()
    rank_macros()
    r3_macros()
    am = agemargin_block()
    table_agemargin(am)
    for name, w in [("table_main.tex", r"\textwidth"), ("table_urban.tex", r"\textwidth"), ("table_urban_ctrl.tex", r"\textwidth"),
                    ("table_ephstd.tex", r"\textwidth"), ("table_class.tex", r"\columnwidth"),
                    ("table_headline.tex", r"\columnwidth"),
                    ("table_ctrl.tex", r"\columnwidth"), ("table_losmodel.tex", r"\columnwidth"),
                    ("table_interf.tex", r"\textwidth"), ("table_agemargin.tex", r"\textwidth"),
                    ("table_recent.tex", r"\textwidth")]:
        fit_width(name, w)
    for name, (title, note) in CAPTIONS.items():
        if os.path.exists(os.path.join(PAPER, name)):
            retitle(name, title, note)
    with open(os.path.join(PAPER, "results_macros.tex"), "w") as f:
        f.write("% generated by sim/analyze.py from " + os.path.basename(RUNS) + "; do not edit\n")
        for k in sorted(macros):
            f.write(f"\\newcommand{{\\{k}}}{{{macros[k]}}}\n")
    print(len(macros), "macros")
