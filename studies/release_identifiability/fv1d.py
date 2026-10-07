"""1D salt column with finite diffusivity, exact in time.

    dc/dt = D d2c/dz2 + S(t)/V - lambda c,      0 < z < H
    top    (z = H): -D dc/dz = k_top c        (release to IV)
    bottom (z = 0):  D dc/dz = k_b c,  k_b = k_wall A_wall / A_top   (to OV)

Putting the whole wall loss on the bottom face keeps the total wall
conductance of the 0D model, so the large-D limit is exactly model0d (checked
in check_large_D). A real crucible loses through its side wall too. That is a
2D effect this 1D rung does not carry.

Space: N uniform finite volumes, half-cell resistance at each Robin face
(face flux = c_cell / (h/2D + 1/k)). Time: the semi-discrete system
dc/dt = M c + s is solved exactly on each interval of constant source and
k_top via the eigendecomposition of the symmetric matrix M.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

import model0d as M0

R_GAS = 8.314462618


def d_calderoni_T(T_K):
    """T in FLiBe, Calderoni et al. 2008 (study brief form):
    D = 9.3e-7 exp(-42 kJ/mol / RT) m2/s."""
    return 9.3e-7 * math.exp(-42000.0 / (R_GAS * T_K))


@lru_cache(maxsize=256)
def _eig(D, k_top, k_b, n, lam):
    h = M0.H_M / n
    g_t = 1.0 / (h / (2 * D) + 1.0 / k_top) if k_top > 0 else 0.0
    g_b = 1.0 / (h / (2 * D) + 1.0 / k_b) if k_b > 0 else 0.0
    a = D / h**2
    Mx = np.zeros((n, n))
    idx = np.arange(n)
    Mx[idx[:-1], idx[1:]] = a
    Mx[idx[1:], idx[:-1]] = a
    Mx[idx, idx] = -2 * a
    Mx[0, 0] = -a - g_b / h
    Mx[-1, -1] = -a - g_t / h
    if n == 1:
        Mx[0, 0] = -(g_b + g_t) / h
    Mx -= lam * np.eye(n)
    w, Q = np.linalg.eigh(Mx)
    return w, Q, g_t, g_b, h


def cumulative_release(D, k_tops, k_wall, source_rate, irradiations, switch_s, times_top, times_wall,
                       n=60, decay=True):
    lam = M0.LAMBDA if decay else 0.0
    k_b = k_wall * M0.A_WALL / M0.A_TOP
    tt = np.asarray(times_top, float)
    tw = np.asarray(times_wall, float)
    edges = {0.0}
    for a, b in irradiations:
        edges.update((max(a, 0.0), max(b, 0.0)))
    if switch_s is not None and len(k_tops) > 1:
        edges.add(float(switch_s))
    edges.update(float(x) for x in np.concatenate([tt, tw]) if x > 0)
    edges = sorted(edges)
    c = np.zeros(n)
    qt = qw = 0.0
    rec = {0.0: (0.0, 0.0)}
    for a, b in zip(edges[:-1], edges[1:]):
        dt = b - a
        mid = 0.5 * (a + b)
        s = sum(source_rate for ia, ib in irradiations if ia <= mid < ib)
        kt = k_tops[1] if (switch_s is not None and len(k_tops) > 1 and mid >= switch_s) else k_tops[0]
        w, Q, g_t, g_b, h = _eig(float(D), float(kt), float(k_b), n, lam)
        a0 = Q.T @ c
        bvec = Q.T @ np.full(n, s / M0.V_M3)
        ainf = -bvec / w
        e = np.exp(w * dt)
        inta = ainf * dt + (a0 - ainf) * (e - 1.0) / w
        c = Q @ (ainf + (a0 - ainf) * e)
        intc = Q @ inta
        qt += M0.A_TOP * g_t * intc[-1]
        qw += M0.A_TOP * g_b * intc[0]
        rec[b] = (qt, qw)
    top = np.array([rec[float(x)][0] if x > 0 else 0.0 for x in tt])
    wall = np.array([rec[float(x)][1] if x > 0 else 0.0 for x in tw])
    return top, wall


def check_large_D():
    """Max relative gap to the 0D model at D = 1e-3 m2/s (Biot ~ 1e-5)."""
    irr = [(0.0, 35580.0), (36000.0, 43620.0)]
    t = np.array([1.1, 2.9, 6.1, 9.2, 14.1, 19.3, 25.1, 32.2, 35.1]) * 86400.0
    args = ([1.0e-7], 1.0e-9, 2.2e5, irr, None, t, t)
    a = cumulative_release(1e-3, *args)
    b = M0.cumulative_release(*args)
    return float(max(np.max(np.abs(a[0] / b[0] - 1)), np.max(np.abs(a[1] / b[1] - 1))))


def check_grid(D=1e-10):
    """Relative change in cumulative top release from n=60 to n=240 at small D."""
    irr = [(0.0, 35580.0), (36000.0, 43620.0)]
    t = np.array([1.1, 2.9, 6.1, 9.2, 14.1, 19.3, 25.1, 32.2, 35.1]) * 86400.0
    a = cumulative_release(D, [1e-7], 1e-9, 2.2e5, irr, None, t, t, n=60)[0]
    b = cumulative_release(D, [1e-7], 1e-9, 2.2e5, irr, None, t, t, n=240)[0]
    return float(np.max(np.abs(a / b - 1)))
