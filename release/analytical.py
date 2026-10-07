"""Closed-form solution of the R1 Dirichlet release problem (SPEC.md section 6).

This is the reference that every code on the release ladder must match to 0.1%
before it touches data. It solves

    dc/dt = D d2c/dz2 + S(t),        0 < z < L,  t > 0
    dc/dz = 0            at z = 0    (bottom of the salt column, no flux)
    c = 0                at z = L    (top surface, Dirichlet)
    c(z, 0) = 0

with S(t) piecewise constant in time: a set of irradiation windows
(t_start, t_stop, S) during which the volumetric source is S, and zero outside.
Every quantity is per unit top area. Multiply by the column's top area to get
totals.

Derivation
----------
The operator d2/dz2 with dc/dz(0) = 0 and c(L) = 0 has the orthogonal
eigenfunctions

    phi_n(z) = cos(mu_n z),   mu_n = (2n + 1) pi / (2L),   n = 0, 1, 2, ...

with integral_0^L phi_n^2 dz = L/2. Expanding c(z, t) = sum_n a_n(t) phi_n(z)
and the uniform source S = sum_n S_n phi_n(z), where

    S_n = (2/L) S integral_0^L cos(mu_n z) dz = 2 S (-1)^n / (L mu_n),

each mode obeys da_n/dt = -k_n a_n + S_n with k_n = D mu_n^2. For a source of
unit strength switched on at t = 0 and left on, a_n(t) = (S_n/k_n)(1 - exp(-k_n t)).
A window (t0, t1, S) is the on-response at t0 minus the on-response at t1,
scaled by S, and windows add (the problem is linear). This is the
superposition the module uses.

The steady parts of the series are summed in closed form so that truncation
only touches the decaying exponentials. Using sum_n 1/mu_n^2 = L^2/2 and
sum_n 1/mu_n^4 = L^4/6 (both from the Dirichlet-Neumann eigenvalue sums), the
unit-strength on-responses at elapsed time tau = t - t_on > 0 are

    c_on(z, tau) = (L^2 - z^2) / (2D)
                   - sum_n [2 (-1)^n / (L D mu_n^3)] exp(-k_n tau) cos(mu_n z)
    J_on(tau)    = L - sum_n [2 / (L mu_n^2)] exp(-k_n tau)        (top flux, -D dc/dz at z = L)
    Q_on(tau)    = L tau - L^3 / (3D) + sum_n [2 / (L D mu_n^4)] exp(-k_n tau)   (cumulative top release)
    I_on(tau)    = L^3 / (3D) - sum_n [2 / (L D mu_n^4)] exp(-k_n tau)          (inventory, integral of c over z)

and all four are zero for tau <= 0. The first term of c_on is the steady
parabola S (L^2 - z^2)/(2D). The mass balance I_on + Q_on = L tau holds mode by
mode, so the analytical solution conserves mass exactly. The checks below test
the implementation of each formula against the others and against a
finite-difference solution that shares none of this algebra.

Convergence
-----------
n_modes is the convergence parameter. For tau > 0 the neglected tail is bounded
by exp(-k_N tau) times a slowly varying factor, so the series converges
exponentially once k_N tau is a few units. The only slow point is tau -> 0+,
where the flux tail decays like 1/N and the cumulative release tail like 1/N^3.
At tau = 0 exactly the module returns the exact value zero instead of the
truncated series.

Checks (run this file)
----------------------
1. Finite differences: a Crank-Nicolson solution with Rannacher start-up on a
   fine grid, written in numpy below, with none of the eigenfunction algebra.
2. Mass balance: the source integral equals the inventory (spatial quadrature of
   c) plus the cumulative release, to 1e-6 relative. The cumulative release is
   also compared with the time quadrature of the flux.

Units are SI throughout: m, s, particles per m3.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_banded


@dataclass(frozen=True)
class Window:
    """One irradiation window: volumetric source `rate` (1/m3/s) on (start, stop]."""

    start: float
    stop: float
    rate: float

    def __post_init__(self):
        if self.stop <= self.start:
            raise ValueError(f"window stop {self.stop} must exceed start {self.start}")


def source_at(t, windows):
    """Piecewise-constant source S(t). A window is on for start < t <= stop.

    The half-open convention matches an implicit time stepper whose steps land
    on the window edges: the step that ends at `stop` is inside the window, the
    step that ends at the next `start` is outside it.
    """
    t = np.asarray(t, dtype=float)
    s = np.zeros_like(t)
    for w in windows:
        s = np.where((t > w.start) & (t <= w.stop), w.rate, s)
    return s


def source_integral(t, windows):
    """Integral of S(t') dt' from 0 to t (1/m3), per unit volume."""
    t = np.asarray(t, dtype=float)
    total = np.zeros_like(t)
    for w in windows:
        total = total + w.rate * np.clip(np.minimum(t, w.stop) - w.start, 0.0, None)
    return total


class SlabSolution:
    """Closed-form series solution described in the module docstring.

    Args:
        D: diffusivity (m2/s)
        L: column depth (m)
        windows: list of Window (start s, stop s, rate 1/m3/s)
        n_modes: number of cosine modes kept in the exponential sums
    """

    def __init__(self, D, L, windows, n_modes=2000):
        if D <= 0 or L <= 0:
            raise ValueError("D and L must be positive")
        self.D = float(D)
        self.L = float(L)
        self.windows = list(windows)
        self.n_modes = int(n_modes)
        n = np.arange(self.n_modes)
        self.mu = (2 * n + 1) * np.pi / (2 * self.L)
        self.k = self.D * self.mu**2
        self.sign = np.where(n % 2 == 0, 1.0, -1.0)

    # unit-strength on-responses, tau is an array of elapsed times

    def _exp(self, tau):
        """exp(-k_n tau) for tau > 0, zero row where tau <= 0. Shape (len(tau), n_modes)."""
        tau = np.atleast_1d(np.asarray(tau, dtype=float))
        e = np.exp(-np.outer(np.clip(tau, 0.0, None), self.k))
        e[tau <= 0.0, :] = 0.0
        return e

    def _on_mask(self, tau):
        return (np.atleast_1d(np.asarray(tau, dtype=float)) > 0.0).astype(float)

    def _flux_on(self, tau):
        e = self._exp(tau)
        return self._on_mask(tau) * self.L - e @ (2.0 / (self.L * self.mu**2))

    def _cumulative_on(self, tau):
        tau_a = np.atleast_1d(np.asarray(tau, dtype=float))
        e = self._exp(tau_a)
        steady = self.L * np.clip(tau_a, 0.0, None) - self._on_mask(tau_a) * self.L**3 / (
            3 * self.D
        )
        return steady + e @ (2.0 / (self.L * self.D * self.mu**4))

    def _inventory_on(self, tau):
        e = self._exp(tau)
        return self._on_mask(tau) * self.L**3 / (3 * self.D) - e @ (
            2.0 / (self.L * self.D * self.mu**4)
        )

    def _concentration_on(self, z, tau):
        """Shape (len(tau), len(z))."""
        z = np.atleast_1d(np.asarray(z, dtype=float))
        e = self._exp(tau)
        parabola = (self.L**2 - z**2) / (2 * self.D)
        coef = 2.0 * self.sign / (self.L * self.D * self.mu**3)
        cosines = np.cos(np.outer(self.mu, z))  # (n_modes, len(z))
        return np.outer(self._on_mask(tau), parabola) - (e * coef) @ cosines

    def _superpose(self, fn, t, *args):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        out = None
        for w in self.windows:
            part = w.rate * (fn(*args, t - w.start) - fn(*args, t - w.stop))
            out = part if out is None else out + part
        if out is None:
            shape = fn(*args, t).shape
            out = np.zeros(shape)
        return out

    # public API

    def concentration(self, z, t):
        """c(z, t) in 1/m3. Returns shape (len(t), len(z)), squeezed for scalars."""
        t_a = np.atleast_1d(np.asarray(t, dtype=float))
        out = None
        for w in self.windows:
            part = w.rate * (
                self._concentration_on(z, t_a - w.start) - self._concentration_on(z, t_a - w.stop)
            )
            out = part if out is None else out + part
        if out is None:
            out = np.zeros((t_a.size, np.atleast_1d(z).size))
        return out.squeeze()

    def flux_top(self, t):
        """Release flux through z = L, -D dc/dz, in 1/m2/s. Positive outward."""
        return self._superpose(lambda tau: self._flux_on(tau), t).squeeze()

    def cumulative_release(self, t):
        """Integral of flux_top from 0 to t, in 1/m2."""
        return self._superpose(lambda tau: self._cumulative_on(tau), t).squeeze()

    def inventory(self, t):
        """Integral of c over the column at time t, in 1/m2."""
        return self._superpose(lambda tau: self._inventory_on(tau), t).squeeze()

    def source_integral(self, t):
        """Source fed into the column up to t, per unit area: L * integral S dt (1/m2)."""
        return self.L * source_integral(t, self.windows)


# Finite-difference cross-check (shares no algebra with the series)


def finite_difference(D, L, windows, t_out, n_cells=2000, dt=60.0):
    """Crank-Nicolson finite-difference solution on a uniform grid.

    Second order in space and time. The first steps after t = 0 and after each
    window edge are taken as implicit-Euler half steps (Rannacher start-up) so
    the step source does not ring. dt must divide every window edge and every
    output time, which the caller arranges by rounding them to dt.

    Returns dict with z, t_out, c (len(t_out) x len(z)), flux, cumulative release
    (trapezoid of the flux), inventory (trapezoid of c), and a second cumulative
    release from mass balance.
    """
    windows = list(windows)
    edges = sorted({w.start for w in windows} | {w.stop for w in windows})
    t_out = np.asarray(t_out, dtype=float)
    for edge in edges + list(t_out):
        if abs(edge / dt - round(edge / dt)) > 1e-9:
            raise ValueError(f"time {edge} is not a multiple of dt={dt}")

    h = L / n_cells
    z = np.linspace(0.0, L, n_cells + 1)
    n_unk = n_cells  # unknowns c_0 .. c_{N-1}, c_N = 0 is the Dirichlet node
    # second-difference operator with the zero-flux ghost node at z = 0
    main = -2.0 * np.ones(n_unk)
    upper = np.ones(n_unk - 1)
    lower = np.ones(n_unk - 1)
    upper[0] = 2.0  # ghost c_{-1} = c_1
    scale = D / h**2

    def banded(theta, step):
        """Banded form of I - theta step D A."""
        ab = np.zeros((3, n_unk))
        ab[0, 1:] = -theta * step * scale * upper
        ab[1, :] = 1.0 - theta * step * scale * main
        ab[2, :-1] = -theta * step * scale * lower
        return ab

    def apply_A(c):
        out = main * c
        out[:-1] += upper * c[1:]
        out[1:] += lower * c[:-1]
        return scale * out

    def flux_top(c):
        # second-order one-sided derivative at z = L with c_N = 0
        dcdz = (-4.0 * c[-1] + c[-2]) / (2 * h)
        return -D * dcdz

    def inventory(c):
        full = np.append(c, 0.0)
        return np.trapezoid(full, z)

    n_steps = int(round(t_out[-1] / dt))
    c = np.zeros(n_unk)
    t = 0.0
    times = [0.0]
    fluxes = [flux_top(c)]
    profiles = {0.0: np.append(c, 0.0)} if 0.0 in t_out else {}
    ab_cn = banded(0.5, dt)
    ab_ie_half = banded(1.0, dt / 2)
    steps_since_switch = 0
    for _ in range(n_steps):
        t_next = t + dt
        s_val = float(source_at(t_next, windows))  # constant on (t, t_next]
        if any(abs(t - edge) < 1e-9 for edge in edges):
            steps_since_switch = 0
        if steps_since_switch < 2:
            # Rannacher: two implicit-Euler half steps in place of one CN step
            for _half in range(2):
                rhs = c + (dt / 2) * s_val
                c = solve_banded((1, 1), ab_ie_half, rhs)
        else:
            rhs = c + 0.5 * dt * apply_A(c) + dt * s_val
            c = solve_banded((1, 1), ab_cn, rhs)
        steps_since_switch += 1
        t = t_next
        times.append(t)
        fluxes.append(flux_top(c))
        if np.any(np.isclose(t, t_out, rtol=0, atol=1e-6)):
            profiles[float(t_out[np.argmin(np.abs(t_out - t))])] = np.append(c, 0.0)

    times = np.asarray(times)
    fluxes = np.asarray(fluxes)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (fluxes[1:] + fluxes[:-1]) * np.diff(times))])
    idx = [int(round(tt / dt)) for tt in t_out]
    c_out = np.array([profiles[float(tt)] for tt in t_out])
    inv_out = np.array([inventory(row[:-1]) for row in c_out])
    q_out = cumulative[idx]
    q_mass_balance = L * source_integral(t_out, windows) - inv_out
    return {
        "z": z,
        "t": t_out,
        "c": c_out,
        "flux": fluxes[idx],
        "cumulative_release": q_out,
        "cumulative_release_mass_balance": q_mass_balance,
        "inventory": inv_out,
        "n_cells": n_cells,
        "dt": dt,
    }


# Checks


def check_mass_balance(sol, times, n_quad=20001):
    """Return max relative |source - (inventory_quadrature + cumulative)| over times.

    The inventory is a Simpson quadrature of c(z, t) over z, not the series
    closed form, so this tests c against Q.
    """
    z = np.linspace(0.0, sol.L, n_quad)
    c = np.atleast_2d(sol.concentration(z, times))
    from scipy.integrate import simpson

    inv = simpson(c, x=z, axis=1)
    q = np.atleast_1d(sol.cumulative_release(times))
    src = np.atleast_1d(sol.source_integral(times))
    return float(np.max(np.abs(src - (inv + q)) / src))


def check_flux_integral(sol, times):
    """Return max relative |Q(t) - integral_0^t J dt'| over times, J integrated by quad."""
    edges = sorted({w.start for w in sol.windows} | {w.stop for w in sol.windows})
    worst = 0.0
    for t in np.atleast_1d(times):
        pts = [e for e in edges if 0.0 < e < t]
        val, _ = quad(lambda tt: float(sol.flux_top(tt)), 0.0, t, points=pts or None, limit=400)
        q = float(sol.cumulative_release(t))
        worst = max(worst, abs(val - q) / q)
    return float(worst)


def run_checks(D, L, windows, sample_times, n_modes=2000, fd_grids=((1000, 60.0), (2000, 30.0))):
    """Run both checks. Returns a dict of figures. Sample times are rounded to the coarsest dt."""
    sol = SlabSolution(D, L, windows, n_modes=n_modes)
    dt_round = max(g[1] for g in fd_grids)
    t_cmp = np.round(np.asarray(sample_times, dtype=float) / dt_round) * dt_round
    out = {
        "D_m2_s": D,
        "L_m": L,
        "n_modes": n_modes,
        "windows": [[w.start, w.stop, w.rate] for w in windows],
        "compare_times_s": t_cmp.tolist(),
        "mass_balance_max_rel_err": check_mass_balance(sol, t_cmp),
        "flux_integral_max_rel_err": check_flux_integral(sol, t_cmp),
        "fd": [],
    }
    q_an = np.atleast_1d(sol.cumulative_release(t_cmp))
    z_an = None
    for n_cells, dt in fd_grids:
        fd = finite_difference(D, L, windows, t_cmp, n_cells=n_cells, dt=dt)
        c_an = np.atleast_2d(sol.concentration(fd["z"], t_cmp))
        c_err = float(np.max(np.abs(fd["c"] - c_an)) / np.max(np.abs(c_an)))
        q_err = float(np.max(np.abs(fd["cumulative_release"] - q_an) / q_an))
        q_mb_err = float(np.max(np.abs(fd["cumulative_release_mass_balance"] - q_an) / q_an))
        out["fd"].append(
            {
                "n_cells": n_cells,
                "dt_s": dt,
                "concentration_max_rel_err": c_err,
                "cumulative_release_max_rel_err": q_err,
                "cumulative_release_mass_balance_max_rel_err": q_mb_err,
            }
        )
        z_an = fd["z"]
    out["fd_finest_cumulative_release_max_rel_err"] = out["fd"][-1]["cumulative_release_max_rel_err"]
    return out


def windows_from_case(case, which_tbr=None):
    """Build Window objects from case.yaml: rate = TBR * Gamma / V (1/m3/s)."""
    which = which_tbr or case["release"]["which_TBR_drives_release"]["value"]
    tbr = case["tbr"]["modelled_mean" if which == "modelled" else "measured"]["value"]
    gamma = case["irradiation"]["neutron_rate"]["value"]
    volume_m3 = case["geometry"]["salt_volume_cm3"]["value"] * 1e-6
    rate = tbr * gamma / volume_m3
    return [
        Window(w["start"]["value"], w["stop"]["value"], rate)
        for w in case["irradiation"]["windows_s"]
    ]


def _main(argv):
    """Run the checks on the run 1 case. Usage: python release/analytical.py [--check|--write] [case.yaml] [D_m2_s]."""
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Analytical solution self-checks")
    parser.add_argument("--check", action="store_true", help="Compare in-memory checks against committed expected results")
    parser.add_argument("--write", action="store_true", help="Write results/release/analytical_check.json")
    parser.add_argument("case_yaml", nargs="?", default=None, help="Path to case.yaml")
    parser.add_argument("D", nargs="?", type=float, default=3.9e-9, help="Diffusivity in m2/s")
    args = parser.parse_args(argv[1:])

    repo = Path(__file__).resolve().parents[1]
    case_path = Path(args.case_yaml) if args.case_yaml else repo / "case.yaml"
    case = yaml.safe_load(case_path.read_text())
    D = args.D
    L = case["release"]["depth_cm"]["value"] * 1e-2
    windows = windows_from_case(case)
    sample_times = np.asarray(case["sampling"]["IV"]["times_day"]["value"]) * 86400.0
    figures = run_checks(D, L, windows, sample_times)

    if figures["mass_balance_max_rel_err"] > 1e-6:
        raise ValueError(
            f"mass balance max rel err {figures['mass_balance_max_rel_err']:.2e} exceeds tolerance 1e-6"
        )
    if figures["fd_finest_cumulative_release_max_rel_err"] > 1e-6:
        raise ValueError(
            f"FD finest cumulative release err {figures['fd_finest_cumulative_release_max_rel_err']:.2e} exceeds tolerance 1e-6"
        )

    out_dir = repo / "results" / "release"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "analytical_check.json"
    payload = {"generated_by": "release/analytical.py, do not edit", **figures}

    if args.check:
        if not out_path.exists():
            print(f"analytical check: {out_path} is missing, run with --write first", file=sys.stderr)
            sys.exit(1)
        committed = json.loads(out_path.read_text())
        # Compare key figures
        diffs = []
        rtol = 1e-4
        for k in ("mass_balance_max_rel_err", "flux_integral_max_rel_err", "fd_finest_cumulative_release_max_rel_err"):
            v_fresh = figures[k]
            v_comm = committed[k]
            if not np.isclose(v_fresh, v_comm, rtol=rtol, atol=1e-8):
                diffs.append(f"{k}: fresh {v_fresh!r} vs committed {v_comm!r}")
        if diffs:
            print("analytical check: fresh run differs from committed expected results:", file=sys.stderr)
            for d in diffs:
                print(f"  {d}", file=sys.stderr)
            sys.exit(1)
        print("analytical check: matches committed expected results")
        print(f"D = {D:.3e} m2/s, L = {L:.5f} m, n_modes = {figures['n_modes']}")
        print(f"mass balance max rel err     : {figures['mass_balance_max_rel_err']:.2e}")
        print(f"flux integral vs Q max rel err: {figures['flux_integral_max_rel_err']:.2e}")
        print(f"FD finest cumulative release : {figures['fd_finest_cumulative_release_max_rel_err']:.2e}")
        return

    out_path.write_text(json.dumps(payload, indent=2))
    print(f"D = {D:.3e} m2/s, L = {L:.5f} m, n_modes = {figures['n_modes']}")
    print(f"mass balance max rel err     : {figures['mass_balance_max_rel_err']:.2e}")
    print(f"flux integral vs Q max rel err: {figures['flux_integral_max_rel_err']:.2e}")
    for fd in figures["fd"]:
        print(
            f"FD n_cells={fd['n_cells']} dt={fd['dt_s']}s: c {fd['concentration_max_rel_err']:.2e}, "
            f"Q {fd['cumulative_release_max_rel_err']:.2e}, "
            f"Q(mass balance) {fd['cumulative_release_mass_balance_max_rel_err']:.2e}"
        )
    print(f"wrote {out_path}")


if __name__ == "__main__":
    _main(sys.argv)
