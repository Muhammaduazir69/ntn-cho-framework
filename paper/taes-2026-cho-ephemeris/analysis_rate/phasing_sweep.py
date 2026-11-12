#!/usr/bin/env python3
"""Nearest-satellite handover rate of Walker shells across the phasing factor F, against
Proposition 2 (eq. near), which describes the phasing average. Writes phasing.json.

    python3 phasing_sweep.py
"""
import json
import os

import numpy as np

import rate_model as M
import walker_mc as W

HERE = os.path.dirname(os.path.abspath(__file__))
SHELLS = [(792, 36), (1584, 72), (3168, 72)]
EPOCHS, H = 6, 2400.0


def nearest_rate(N, P, F, rng):
    raan, u0, n, r, i = W.shell(N, P, F, 550, 53)
    out = []
    for _ in range(EPOCHS):
        u = u0 + n * rng.uniform(0, 86400)
        t = np.arange(0, H, 1.0)
        E = W.elev(t, raan, u, n, r, i, np.radians(30), rng.uniform(0, 2 * np.pi))
        out.append(60 * np.count_nonzero(np.diff(np.argmax(E, axis=0))) / H)
    return float(np.mean(out))


if __name__ == "__main__":
    res = {}
    for N, P in SHELLS:
        rng = np.random.default_rng(N)
        Fs = list(range(0, P, max(1, P // 12)))
        rates = [nearest_rate(N, P, F, rng) for F in Fs]
        an = M.rates(N, 550, 53, 25)["Hnear"]
        res[str(N)] = dict(P=P, F=Fs, rates=rates, mean=float(np.mean(rates)), lo=min(rates), hi=max(rates),
                           an=float(an))
        print(N, P, f"phasing mean {np.mean(rates):.3f} [{min(rates):.3f}, {max(rates):.3f}] vs analysis {an:.3f}",
              flush=True)
    json.dump(res, open(os.path.join(HERE, "phasing.json"), "w"), indent=1)


def exit_spread(N=3168, P=72, Fs=(1, 17, 35), eQ=25, epochs=8):
    """Exit-triggered rate of the same shell at three phasings (it hardly moves)."""
    out = {}
    for F in Fs:
        rng = np.random.default_rng(5)
        res = []
        for _ in range(epochs):
            _, rr = W.run(N, P, F, 550, 53, eQ, H, rng, 30, dt=1.0)
            res += rr[1:]
        out[F] = 60 / float(np.mean(res))
    return out
