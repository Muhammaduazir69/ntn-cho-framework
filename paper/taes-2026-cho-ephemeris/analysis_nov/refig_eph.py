"""Redraw Fig. S9 (fig_eph) with the panel label 'element-set error', same data and code."""
import os, sys
SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SIM)
os.chdir(SIM)
import analyze as A
A.FIG = os.path.join(os.path.dirname(SIM), "figures_v8")
src = open(os.path.join(SIM, "analyze.py")).read()
eo, ee = A.eph_block("eph"), A.eph_block("ephemp")
import matplotlib.pyplot as plt
_orig = plt.Axes.set_title
def _t(self, label, *a, **k):
    return _orig(self, label.replace("measured error", "element-set error"), *a, **k)
plt.Axes.set_title = _t
A.fig_eph(eo, ee)
print("ok")
