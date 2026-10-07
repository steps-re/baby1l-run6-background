"""0D well-mixed release model with an optional sweep-gas step in k_top.

Same ODE as release/r0.py and the lab's libra_toolbox Model:

    V dc/dt = TBR * Gamma(t) - A_top k_top(t) c - A_wall k_wall c - lambda V c

k_top(t) is piecewise constant over gas periods (one value before a gas
switch, another after). With no switch this is exactly R0Model, which
check_against_r0() confirms.

Geometry is the lab's (libra_toolbox Model at cb41e93): radius 7 cm,
V = 1 L, height = V / (pi r^2), wall = 2 pi (r + 0.06 in) h + pi r^2.
"""

from __future__ import annotations

import math

import numpy as np

HALF_LIFE_S = 4500.0 * 86400.0  # release/r0.py
LAMBDA = math.log(2.0) / HALF_LIFE_S

R_M = 0.07
V_M3 = 1.0e-3
A_TOP = math.pi * R_M**2
H_M = V_M3 / A_TOP
L_WALL_M = 0.06 * 0.0254
A_WALL = 2.0 * math.pi * (R_M + L_WALL_M) * H_M + A_TOP


def cumulative_release(k_tops, k_wall, source_rate, irradiations, switch_s, times_top, times_wall,
                       decay=True, pause=None):
    """Cumulative particles released to the top (IV) and wall (OV).

    k_tops: sequence of k_top values, one per gas period (len 1 or 2).
    source_rate: TBR * Gamma, particles/s during irradiation.
    irradiations: list of (start_s, stop_s).
    switch_s: gas switch time in s, or None.
    Times <= 0 give zero release.
    pause: optional (start_s, stop_s) during which nothing is released (k_top =
    k_wall = 0), e.g. while the salt is frozen. Production and decay continue.
    """
    lam = LAMBDA if decay else 0.0
    tt = np.asarray(times_top, float)
    tw = np.asarray(times_wall, float)
    edges = {0.0}
    for a, b in irradiations:
        edges.update((max(a, 0.0), max(b, 0.0)))
    if switch_s is not None and len(k_tops) > 1:
        edges.add(float(switch_s))
    if pause is not None:
        edges.update((max(float(pause[0]), 0.0), max(float(pause[1]), 0.0)))
    allt = np.concatenate([tt, tw])
    edges.update(float(x) for x in allt if x > 0)
    edges = sorted(edges)

    def ktop_at(mid):
        if switch_s is not None and len(k_tops) > 1 and mid >= switch_s:
            return k_tops[1]
        return k_tops[0]

    c = 0.0
    qt = qw = 0.0
    record = {0.0: (0.0, 0.0)}
    for a, b in zip(edges[:-1], edges[1:]):
        dt = b - a
        mid = 0.5 * (a + b)
        s = 0.0
        for ia, ib in irradiations:
            if ia <= mid < ib:
                s = source_rate
        kt = ktop_at(mid)
        kw = k_wall
        if pause is not None and pause[0] <= mid < pause[1]:
            kt = kw = 0.0
        alpha = (A_TOP * kt + A_WALL * kw) / V_M3 + lam
        sv = s / V_M3
        if alpha > 0:
            e = math.exp(-alpha * dt)
            cinf = sv / alpha
            intc = cinf * dt + (c - cinf) * (1.0 - e) / alpha
            c = cinf + (c - cinf) * e
        else:
            intc = c * dt + 0.5 * sv * dt * dt
            c = c + sv * dt
        qt += A_TOP * kt * intc
        qw += A_WALL * kw * intc
        record[b] = (qt, qw)
    top = np.array([record[float(x)][0] if x > 0 else 0.0 for x in tt])
    wall = np.array([record[float(x)][1] if x > 0 else 0.0 for x in tw])
    return top, wall


def check_against_r0():
    """Max relative gap to release/r0.py R0Model for a constant k_top."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from release.r0 import IrradiationWindow, R0Model

    irr = [(0.0, 35580.0), (36000.0, 43620.0)]
    t = np.array([1.1, 2.9, 6.1, 9.2, 14.1, 19.3, 25.1, 32.2, 35.1]) * 86400.0
    kt, kw, src = 9.0e-8, 1.5e-9, 2.5e-3 * 8.65e7
    m = R0Model(V_M3, A_TOP, kt, A_WALL, kw, TBR=2.5e-3,
                windows=[IrradiationWindow(a, b, 8.65e7) for a, b in irr], decay=True,
                decay_constant=LAMBDA)
    ref = m.solve(t)
    top, wall = cumulative_release([kt], kw, src, irr, None, t, t)
    return float(max(np.max(np.abs(top / ref["cumulative_top_release"] - 1)),
                     np.max(np.abs(wall / ref["cumulative_wall_release"] - 1))))
