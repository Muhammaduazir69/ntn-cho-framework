#!/usr/bin/env python3
"""Exact-geometry Monte Carlo of exit-triggered handover (MST and contour rules) under line of
sight on a Walker-Delta shell: circular orbits, Earth rotation, spherical Earth for elevation.
Used to separate the Poisson approximation of the handover-rate analysis from horizon effects.

    python3 walker_mc.py N P F h inc eQ horizon_s runs [lat]
"""
import sys
import numpy as np

R, MU, WE = 6378.137, 398600.4418, 7.2921159e-5


def shell(N, P, F, h, inc):
    S = N // P
    p = np.repeat(np.arange(P), S)
    s = np.tile(np.arange(S), P)
    raan = 2 * np.pi * p / P
    u0 = 2 * np.pi * s / S + 2 * np.pi * F * p / N
    r = R + h
    return raan, u0, np.sqrt(MU / r ** 3), r, np.radians(inc)


def elev(t, raan, u0, n, r, i, lat, lon0):
    u = u0[:, None] + n * t[None, :]
    O = raan[:, None]
    x = r * (np.cos(O) * np.cos(u) - np.sin(O) * np.sin(u) * np.cos(i))
    y = r * (np.sin(O) * np.cos(u) + np.cos(O) * np.sin(u) * np.cos(i))
    z = r * np.sin(u) * np.sin(i)
    th = lon0 + WE * t[None, :]
    ux, uy, uz = np.cos(lat) * np.cos(th), np.cos(lat) * np.sin(th), np.sin(lat) * np.ones_like(th)
    dx, dy, dz = x - R * ux, y - R * uy, z - R * uz
    d = np.sqrt(dx * dx + dy * dy + dz * dz)
    return np.degrees(np.arcsin((dx * ux + dy * uy + dz * uz) / d))


def run(N, P, F, h, inc, eQ, H, rng, lat=30.0, dt=0.5, first="strongest"):
    raan, u0, n, r, i = shell(N, P, F, h, inc)
    shift = rng.uniform(0, 86400)
    u0 = u0 + n * shift
    lon0 = rng.uniform(0, 2 * np.pi)
    tmax = H + 1500
    t = np.arange(0, tmax, dt)
    E = elev(t, raan, u0, n, r, i, np.radians(lat), lon0)
    vis = E >= eQ
    # remaining time inside the contour, per satellite and sample
    rem = np.zeros_like(E)
    run_ = np.zeros(E.shape[0])
    for k in range(len(t) - 1, -1, -1):
        run_ = np.where(vis[:, k], run_ + dt, 0.0)
        rem[:, k] = run_
    k = 0
    if first == "strongest":
        s = int(np.argmax(E[:, 0]))
    else:
        s = int(np.argmax(rem[:, 0]))
    ho, res, kk = 0, [], 0
    while True:
        kexit = k + int(rem[s, k] / dt)
        # stay past the contour until some satellite is inside it
        while kexit < len(t) and not vis[:, kexit].any():
            kexit += 1
        if t[min(kexit, len(t) - 1)] >= H:
            break
        c = int(np.argmax(rem[:, kexit]))
        res.append(t[kexit] - t[k])
        s, k = c, kexit
        ho += 1
    return ho, res


if __name__ == "__main__":
    N, P, F = map(int, sys.argv[1:4])
    h, inc, eQ, H = map(float, sys.argv[4:8])
    runs = int(sys.argv[8])
    lat = float(sys.argv[9]) if len(sys.argv) > 9 else 30.0
    rng = np.random.default_rng(1)
    hos, ress = [], []
    for _ in range(runs):
        ho, res = run(N, P, F, h, inc, eQ, H, rng, lat)
        hos.append(ho)
        ress += res[1:]
    print(f"rate {60*np.mean(hos)/H:.4f} per min  (+-{60*1.96*np.std(hos)/np.sqrt(runs)/H:.4f}); "
          f"mean residence after first {np.mean(ress):.1f} s -> {60/np.mean(ress):.4f} per min")
