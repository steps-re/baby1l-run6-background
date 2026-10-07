"""Likelihood, fits, profiles and Fisher information for the 0D release model.

Data are per-sample increments (Bq collected by each bubbler change), not the
cumulative curve, so counting errors are independent except where two samples
share a blank. That sharing is in the counting covariance from
build_increments.py.

Noise model, per stream s (IV, OV):
    Sigma = C_counting + diag(sig_s^2)
C_counting is the Poisson covariance from the LSC files. sig_s (Bq) is a
fitted additive floor per stream for everything counting does not cover:
bubbler efficiency, vial handling, and model error. With noise="counting" the
floors are off. Two alternatives were tried and rejected (see README): a floor
proportional to the data (lets the fit discount large increments and chase
near-zero ones; run 4 k_top fell to 2e-9 m/s) and one proportional to the
model (local minima; run 4 stuck 960 NLL units above the additive optimum).

Parameters are log10 of the physical quantities. Parametrization "A" is
(k_top per gas period, k_wall, [TBR scale], sig_IV, [sig_OV]). Parametrization "B"
swaps (k_top of the first period, k_wall) for (alpha, logit phi), where
alpha = (A_top k_top + A_wall k_wall)/V is the total loss rate and
phi = A_top k_top / (A_top k_top + A_wall k_wall) is the top share.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

import model0d as M

HERE = Path(__file__).resolve().parent
INC = json.loads((HERE / "increments.json").read_text())
# quench-matched background deltas (quench.py); absent until quench.py has run
QUENCH = json.loads((HERE / "quench.json").read_text()) if (HERE / "quench.json").exists() else {"runs": {}}
BQ_PER_PARTICLE = INC["Bq_per_particle"]

BOUNDS = {
    "log10_ktop": (-11.0, -3.0),
    "log10_ktop_p2": (-11.0, -3.0),
    "log10_kwall": (-13.0, -4.0),
    # TBR scale: measured TBR is 0.92-1.93x the modelled TBR across runs 1-4 and
    # the neutron rate errors are 2-20%, so a factor of 3.16 either way is a
    # generous nuisance range. Wider lets fits trade a 5x production for a
    # noise floor that swallows the curve (seen for run 4 at small D).
    "log10_tbr_scale": (-0.5, 0.5),
    "log10_sig_IV": (-4.0, 1.0),
    "log10_sig_OV": (-4.0, 1.0),
    "log10_alpha": (-10.0, -2.0),
    "logit_phi": (-14.0, 14.0),
}
CHI2_95 = 3.841458820694124
CHI2_68 = 0.9892518671359535  # 1 dof, 68.27%


@dataclass
class Dataset:
    name: str
    run: int
    t_max_day: float | None = None  # keep samples with t <= t_max
    step: bool = False  # k_top steps at the gas switch
    streams: tuple = ("IV", "OV")
    tbr: str = "fixed"  # "fixed" at the lab's measured TBR, "lab" at the lab's TBR_used_in_model, or "free"
    noise: str = "fitted"  # "fitted" or "counting"
    background: str = "published"  # "published" (the lab's numbers) or "quench" (quench.json, runs 3 and 4)
    pause: tuple | None = None  # (start_s, stop_s) with no release, see model0d.cumulative_release
    # filled in __post_init__
    y: np.ndarray = field(init=False, repr=False)
    C: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        r = INC["runs"][str(self.run)]
        self.meta = r
        order = r["total_vector_order"]
        keep, ys, sid, tt = [], [], [], {}
        for s in self.streams:
            st = r["streams"][s]
            t = np.array(st["times_day"])
            inc = np.array(st["inc_total_Bq"])
            if self.background == "quench":
                q = QUENCH["runs"][str(self.run)]
                if q.get("anchored_delta_Bq") is None:
                    raise ValueError(f"run {self.run}: no quench-matched background ({q.get('reason')})")
                inc = inc + np.array(q["anchored_delta_Bq"][s])
            m = np.ones_like(t, bool) if self.t_max_day is None else (t <= self.t_max_day)
            tt[s] = t[m] * 86400.0
            for j in np.where(m)[0]:
                keep.append(order.index(f"{s}[{j}]"))
                ys.append(inc[j])
                sid.append(s)
        C = np.array(r["cov_total_Bq2"])
        self.C = C[np.ix_(keep, keep)]
        self.y = np.array(ys)
        self.sid = np.array(sid)
        self.times = tt
        self.irr = [tuple(x) for x in r["irradiations_s"]]
        self.gamma = r["neutron_rate_n_s"]
        # "lab" uses the TBR the lab's own model used (processed_data.json TBR_used_in_model):
        # the measured TBR for runs 1, 2, 4 and 6, the modelled TBR for run 3.
        self.tbr_meas = r["TBR_used_in_model"] if self.tbr == "lab" else r["measured_TBR"]
        self.switch = r["gas_switch_s"] if self.step else None
        self.names_A = (["log10_ktop"] + (["log10_ktop_p2"] if self.step else []) + ["log10_kwall"]
                        + (["log10_tbr_scale"] if self.tbr == "free" else [])
                        + ([f"log10_sig_{s}" for s in self.streams] if self.noise == "fitted" else []))
        # names_B: alpha, phi, then every A parameter other than k_top (period 1) and k_wall
        rest = [n for n in self.names_A if n not in ("log10_ktop", "log10_kwall")]
        self.names_B = ["log10_alpha", "logit_phi"] + rest

    def n(self):
        return len(self.y)

    # --- parameter handling -------------------------------------------------
    def names(self, par="A"):
        return self.names_A if par == "A" else self.names_B

    def physical(self, theta, par="A"):
        d = dict(zip(self.names(par), theta))
        if par == "B":
            alpha = 10 ** d["log10_alpha"]
            phi = 1.0 / (1.0 + math.exp(-d["logit_phi"]))
            kt = phi * alpha * M.V_M3 / M.A_TOP
            kw = (1 - phi) * alpha * M.V_M3 / M.A_WALL
        else:
            kt = 10 ** d["log10_ktop"]
            kw = 10 ** d["log10_kwall"]
        ktops = [kt] + ([10 ** d["log10_ktop_p2"]] if self.step else [])
        s = 10 ** d.get("log10_tbr_scale", 0.0)
        f = {st: 10 ** d.get(f"log10_sig_{st}", -np.inf) for st in self.streams}
        return ktops, kw, s, f

    def to_B(self, thetaA):
        d = dict(zip(self.names_A, thetaA))
        kt, kw = 10 ** d["log10_ktop"], 10 ** d["log10_kwall"]
        at, aw = M.A_TOP * kt, M.A_WALL * kw
        d["log10_alpha"] = math.log10((at + aw) / M.V_M3)
        phi = at / (at + aw)
        phi = min(max(phi, 1e-15), 1 - 1e-15)
        d["logit_phi"] = math.log(phi / (1 - phi))
        return np.array([d[n] for n in self.names_B])

    def to_A(self, thetaB):
        ktops, kw, _, _ = self.physical(thetaB, "B")
        d = dict(zip(self.names_B, thetaB))
        d["log10_ktop"] = math.log10(ktops[0])
        d["log10_kwall"] = math.log10(max(kw, 1e-300))
        return np.array([d[n] for n in self.names_A])

    # --- model and likelihood -------------------------------------------------
    def predict(self, theta, par="A"):
        ktops, kw, s, _ = self.physical(theta, par)
        src = s * self.tbr_meas * self.gamma
        top, wall = M.cumulative_release(ktops, kw, src, self.irr, self.switch,
                                         self.times.get("IV", []), self.times.get("OV", []),
                                         pause=self.pause)
        out = []
        for st, cum in (("IV", top), ("OV", wall)):
            if st in self.streams:
                out.append(np.diff(np.concatenate([[0.0], cum])) * BQ_PER_PARTICLE)
        return np.concatenate(out)

    def sigma(self, theta, par="A"):
        """Full noise covariance at theta (counting + fitted floors)."""
        if self.noise != "fitted":
            return self.C.copy()
        _, _, _, f = self.physical(theta, par)
        fr = np.array([f[s] for s in self.sid])
        return self.C + np.diag(fr**2)

    def nll(self, theta, par="A", y=None):
        y = self.y if y is None else y
        mu = self.predict(theta, par)
        S = self.sigma(theta, par)
        r = y - mu
        try:
            cf = cho_factor(S, lower=True)
        except np.linalg.LinAlgError:
            return 1e12
        return 0.5 * float(r @ cho_solve(cf, r)) + float(np.sum(np.log(np.diag(cf[0])))) \
            + 0.5 * len(r) * math.log(2 * math.pi)

    def chi2(self, theta, par="A", y=None):
        y = self.y if y is None else y
        r = y - self.predict(theta, par)
        return float(r @ np.linalg.solve(self.C, r))

    def bounds(self, par="A"):
        return [BOUNDS[n] for n in self.names(par)]

    def starts(self):
        """Deterministic multi-start grid in parametrization A, not centred on
        any fitted value: k_top spans 1e-8..1e-5, k_wall 1e-11..1e-8."""
        out = []
        for kt in (-8.0, -7.0, -6.0, -5.0):
            for kw in (-11.0, -9.0):
                th = []
                for n in self.names_A:
                    th.append({"log10_ktop": kt, "log10_ktop_p2": kt, "log10_kwall": kw,
                               "log10_tbr_scale": 0.0, "log10_sig_IV": -1.0,
                               "log10_sig_OV": -1.0}[n])
                out.append(np.array(th))
        return out


POLISH = {"on": True}  # Powell polish after L-BFGS-B; analysis turns it off in loops


def _min(fun, x0, bounds, polish=None):
    r = minimize(fun, x0, method="L-BFGS-B", bounds=bounds,
                 options={"maxiter": 3000, "ftol": 1e-12, "gtol": 1e-8})
    if not (POLISH["on"] if polish is None else polish):
        return r
    # polish with Powell from the L-BFGS-B point (robust to flat noise directions)
    r2 = minimize(fun, r.x, method="Powell", bounds=bounds,
                  options={"maxiter": 20000, "xtol": 1e-6, "ftol": 1e-12})
    return r2 if r2.fun <= r.fun else r


def fit(ds: Dataset, y=None, starts=None):
    best = None
    fun = lambda th: ds.nll(th, "A", y)  # noqa: E731
    for x0 in (starts or ds.starts()):
        r = _min(fun, x0, ds.bounds("A"))
        if best is None or r.fun < best.fun:
            best = r
    return best.x, float(best.fun)


def constrained_fit(ds: Dataset, par, name, value, x_start, y=None, extra_starts=()):
    """Minimise the NLL with parameter `name` (in parametrization `par`) fixed.
    `name` and `value` may also be sequences to fix several parameters."""
    names = ds.names(par)
    fixed_names = [name] if isinstance(name, str) else list(name)
    fixed_vals = [value] if isinstance(name, str) else list(np.atleast_1d(value))
    fix = {names.index(n): v for n, v in zip(fixed_names, fixed_vals)}
    free = [j for j in range(len(names)) if j not in fix]
    bnds = [ds.bounds(par)[j] for j in free]

    def full(z):
        th = np.empty(len(names))
        for j, v in fix.items():
            th[j] = v
        th[free] = z
        return th

    fun = lambda z: ds.nll(full(z), par, y)  # noqa: E731
    best = None
    for xs in [x_start, *extra_starts]:
        z0 = np.clip(np.asarray(xs, float)[free], [b[0] for b in bnds], [b[1] for b in bnds])
        r = _min(fun, z0, bnds)
        if best is None or r.fun < best.fun:
            best = r
    return full(best.x), float(best.fun)


def profile(ds: Dataset, thetaA, nll_min, name, par=None, half_width=3.0, n=31, y=None, refine=True):
    """Profile likelihood of one parameter. Returns grid, 2*dNLL, and CIs.

    1. Sweep a grid outwards from the optimum with warm starts.
    2. Repair: wherever the profile drops while moving away from the optimum
       (a sign of a missed minimum), refit that point from the full start grid.
    3. Refine: find each 95% edge by root-finding on constrained fits, not by
       interpolating between grid points.
    """
    par = par or ("B" if name in ("log10_alpha", "logit_phi") else "A")
    th_hat = thetaA if par == "A" else ds.to_B(thetaA)
    i = ds.names(par).index(name)
    lo_b, hi_b = BOUNDS[name]
    c = th_hat[i]
    lo, hi = max(lo_b, c - half_width), min(hi_b, c + half_width)
    grid = np.unique(np.concatenate([np.linspace(lo, hi, n), [c]]))
    prof = np.full(len(grid), np.nan)
    ic = int(np.argmin(np.abs(grid - c)))
    thetas = [None] * len(grid)
    for direction in (range(ic, len(grid)), range(ic, -1, -1)):
        prev = th_hat
        for k in direction:
            th, v = constrained_fit(ds, par, name, grid[k], prev, y, extra_starts=(th_hat,))
            if not np.isfinite(prof[k]) or v < prof[k]:
                prof[k], thetas[k] = v, th
            prev = thetas[k]
    all_starts = [s_ if par == "A" else ds.to_B(s_) for s_ in ds.starts()]
    repaired = 0
    for direction in (range(ic + 1, len(grid)), range(ic - 1, -1, -1)):
        prev_k = ic
        for k in direction:
            if prof[k] < prof[prev_k] - 1e-6:
                th, v = constrained_fit(ds, par, name, grid[k], thetas[prev_k], y,
                                        extra_starts=(th_hat, *all_starts))
                if v < prof[k]:
                    prof[k], thetas[k] = v, th
                repaired += 1
            prev_k = k
    m = min(nll_min, float(np.nanmin(prof)))
    d2 = 2.0 * (prof - m)
    pinned = [_pinned(ds, par, i, th, y) for th in thetas]
    ci95 = _crossings(grid, d2, CHI2_95, pinned)
    if refine:
        ci95 = _refine(ds, par, name, grid, d2, thetas, m, y, ci95, CHI2_95)
    return {"name": name, "par": par, "grid": grid.tolist(), "two_delta_nll": d2.tolist(),
            "pinned_nuisance": pinned, "repaired_points": repaired,
            "mle": float(grid[int(np.argmin(d2))]) if nll_min > np.nanmin(prof) + 1e-6 else float(c),
            "ci95": ci95,
            "ci68": _crossings(grid, d2, CHI2_68, pinned),
            "bound_lo": lo_b, "bound_hi": hi_b}


def _refine(ds, par, name, grid, d2, thetas, m, y, ci, level):
    """Root-find each closed 95% edge between its bracketing grid points."""
    from scipy.optimize import brentq
    out = dict(ci)
    for side in ("lo", "hi"):
        x = ci[side]
        if x is None:
            continue
        k_out = int(np.searchsorted(grid, x)) if side == "hi" else int(np.searchsorted(grid, x)) - 1
        k_in = k_out - 1 if side == "hi" else k_out + 1
        if not (0 <= k_out < len(grid) and 0 <= k_in < len(grid)):
            continue
        start = thetas[k_in]

        def g(v):
            _, val = constrained_fit(ds, par, name, v, start, y, extra_starts=(thetas[k_out],))
            return 2 * (val - m) - level

        try:
            a_, b_ = sorted((grid[k_in], grid[k_out]))
            out[side] = float(brentq(g, a_, b_, xtol=1e-4))
        except ValueError:
            pass  # bracket lost its sign change after refitting: keep the interpolated edge
    return out


def _pinned(ds, par, i, th, y=None, step=0.1, tol=1e-3):
    """Nuisance parameters held by a bound at this profile point: on the bound,
    and moving them `step` decades back inside raises the NLL (the bound binds).
    A parameter on a bound along a flat direction does not bind and is not
    listed. Noise floors are ignored (a floor at ~0 limits nothing). k_wall on
    its lower bound means k_wall = 0, the physical limit."""
    out = []
    f0 = ds.nll(th, par, y)
    for j, n in enumerate(ds.names(par)):
        if j == i or n.startswith("log10_sig_"):
            continue
        lo, hi = BOUNDS[n]
        if th[j] - lo < 1e-3:
            d = +step
        elif hi - th[j] < 1e-3:
            d = -step
        else:
            continue
        t2 = np.array(th, float)
        t2[j] += d
        if ds.nll(t2, par, y) - f0 > tol:
            out.append(n)
    return out


def _crossings(grid, d2, level, pinned):
    """Where 2*dNLL crosses `level` on each side of the minimum.

    A side counts as closed by the data only if (a) no nuisance parameter is
    binding at a bound at the grid points either side of the crossing, and
    (b) the profile does not dip back below `level` further out. A binding
    k_wall lower bound is the physical limit k_wall = 0 (reported, not
    discounted). A binding TBR-scale bound is the prior range and does not
    count as information from the data.
    """
    i0 = int(np.argmin(d2))
    lo = hi = None
    lo_by = hi_by = []
    lo_re = hi_re = False
    for k in range(i0, 0, -1):
        if d2[k - 1] >= level > d2[k] or (d2[k - 1] >= level and d2[k] < level):
            lo = float(np.interp(level, [d2[k], d2[k - 1]], [grid[k], grid[k - 1]]))
            lo_by = sorted(set(pinned[k - 1]) | set(pinned[k]))
            lo_re = bool(np.any(np.asarray(d2[: k - 1]) < level))
            break
    for k in range(i0, len(grid) - 1):
        if d2[k + 1] >= level > d2[k] or (d2[k + 1] >= level and d2[k] < level):
            hi = float(np.interp(level, [d2[k], d2[k + 1]], [grid[k], grid[k + 1]]))
            hi_by = sorted(set(pinned[k + 1]) | set(pinned[k]))
            hi_re = bool(np.any(np.asarray(d2[k + 2:]) < level))
            break
    prior = lambda by: [n for n in by if n != "log10_kwall"]  # noqa: E731
    lo_data = lo is not None and not prior(lo_by) and not lo_re
    hi_data = hi is not None and not prior(hi_by) and not hi_re
    if lo_data and hi_data:
        status = "identifiable"
    elif not lo_data and not hi_data:
        status = "not identifiable"
    elif hi_data:
        status = "upper bound only"
    else:
        status = "lower bound only"
    return {"lo": lo, "hi": hi, "status": status, "lo_set_by_bound": lo_by, "hi_set_by_bound": hi_by,
            "lo_reenters": lo_re, "hi_reenters": hi_re}


def fisher(ds: Dataset, thetaA, h=1e-3):
    """Observed Fisher information (Hessian of the NLL) in parametrization A."""
    k = len(thetaA)
    H = np.zeros((k, k))
    f0 = ds.nll(thetaA)
    for i in range(k):
        for j in range(i, k):
            ei = np.zeros(k); ei[i] = h  # noqa: E702
            ej = np.zeros(k); ej[j] = h  # noqa: E702
            if i == j:
                H[i, i] = (ds.nll(thetaA + ei) - 2 * f0 + ds.nll(thetaA - ei)) / h**2
            else:
                H[i, j] = H[j, i] = (ds.nll(thetaA + ei + ej) - ds.nll(thetaA + ei - ej)
                                     - ds.nll(thetaA - ei + ej) + ds.nll(thetaA - ei - ej)) / (4 * h * h)
    return H


def fisher_summary(ds: Dataset, thetaA):
    H = fisher(ds, thetaA)
    names = ds.names_A
    at_bound = [n for n, v in zip(names, thetaA)
                if abs(v - BOUNDS[n][0]) < 1e-3 or abs(v - BOUNDS[n][1]) < 1e-3]
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        cov, se = None, np.full(len(names), np.nan)
    # physics block after marginalising the noise terms: inverse of the cov block
    phys = [i for i, n in enumerate(names) if not n.startswith("log10_sig_")]
    eig = None
    if cov is not None:
        blk = np.linalg.inv(cov[np.ix_(phys, phys)])
        w, v = np.linalg.eigh(blk)
        eig = {"eigenvalues": w.tolist(), "eigenvectors": v.T.tolist(),
               "names": [names[i] for i in phys],
               "condition_number": float(w[-1] / w[0]) if w[0] > 0 else float("inf")}
    return {"names": names, "se_log10": se.tolist(), "wald95_halfwidth_log10": (1.96 * se).tolist(),
            "hessian_min_eig": float(np.linalg.eigvalsh(H)[0]), "params_at_bound": at_bound,
            "physics_block": eig}
