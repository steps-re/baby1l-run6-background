"""Well-mixed 0D tritium release model (SPEC.md section 6).

Solves the ordinary differential equation

    V dc/dt = TBR * Gamma(t) - A_top * k_top * c - A_wall * k_wall * c - lambda * V * c

where
    c(t) is tritium concentration (particles/m3)
    V is salt volume (m3)
    TBR is tritium breeding ratio (particles/neutron)
    Gamma(t) is neutron rate (neutrons/s), piecewise constant in time
    A_top is top surface area (m2)
    k_top is top mass transfer coefficient (m/s)
    A_wall is crucible wall contact area (m2)
    k_wall is wall mass transfer coefficient (m/s)
    lambda is the radioactive decay constant of tritium (1/s)

The solution is piecewise analytical across irradiation windows. An independent
numerical solve with solve_ivp checks the algebra.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.integrate import solve_ivp

# Tritium half-life: 4500 days = 12.32 years (Lucas and Unterweger 2000, case.yaml)
HALF_LIFE_DAYS = 4500.0
HALF_LIFE_S = HALF_LIFE_DAYS * 86400.0
DEFAULT_LAMBDA = math.log(2.0) / HALF_LIFE_S


@dataclass(frozen=True)
class IrradiationWindow:
    """One irradiation window: neutron rate (n/s) on (start, stop]."""

    start: float
    stop: float
    neutron_rate: float

    def __post_init__(self):
        if self.stop <= self.start:
            raise ValueError(f"stop {self.stop} must exceed start {self.start}")


class R0Model:
    """Well-mixed 0D release model.

    Args:
        volume: salt volume V in m3
        A_top: top surface area in m2
        k_top: top mass transfer coefficient in m/s
        A_wall: wall contact area in m2
        k_wall: wall mass transfer coefficient in m/s
        TBR: tritium breeding ratio (particles/neutron)
        windows: list of IrradiationWindow
        decay: whether to include tritium radioactive decay
        decay_constant: custom decay constant in 1/s. Defaults to 12.32 y half-life.
    """

    def __init__(
        self,
        volume: float,
        A_top: float,
        k_top: float,
        A_wall: float = 0.0,
        k_wall: float = 0.0,
        TBR: float = 1.0,
        windows: Sequence[IrradiationWindow] | None = None,
        decay: bool = True,
        decay_constant: float | None = None,
    ):
        if volume <= 0:
            raise ValueError("volume must be positive")
        if A_top < 0 or k_top < 0 or A_wall < 0 or k_wall < 0:
            raise ValueError("areas and transfer coefficients must be non-negative")

        self.V = float(volume)
        self.A_top = float(A_top)
        self.k_top = float(k_top)
        self.A_wall = float(A_wall)
        self.k_wall = float(k_wall)
        self.TBR = float(TBR)
        self.windows = list(windows or [])
        if decay:
            self.decay_constant = DEFAULT_LAMBDA if decay_constant is None else float(decay_constant)
        else:
            self.decay_constant = 0.0

        # Loss rate constant alpha (1/s)
        self.alpha = (self.A_top * self.k_top + self.A_wall * self.k_wall) / self.V + self.decay_constant
        self.tau = 1.0 / self.alpha if self.alpha > 0 else float("inf")

    @classmethod
    def from_case(
        cls,
        case: dict,
        which_tbr: str | None = None,
        decay: bool = True,
    ) -> R0Model:
        """Construct an R0Model from a loaded case.yaml dictionary."""
        geom = case["geometry"]
        rel = case["release"]
        irrad = case["irradiation"]

        radius_m = (
            geom.get("lab_0d_model_radius_cm", geom.get("salt_radius_cm"))["value"] * 1e-2
        )
        height_m = (
            geom.get("lab_0d_model_height_cm", geom.get("salt_thickness_cm"))["value"] * 1e-2
        )
        l_wall_m = rel.get("lab_0d_wall_thickness_inch", {}).get("value", 0.06) * 0.0254
        a_top_m2 = math.pi * radius_m**2
        a_wall_m2 = 2.0 * math.pi * (radius_m + l_wall_m) * height_m + a_top_m2
        volume_m3 = a_top_m2 * height_m

        which = which_tbr or rel["which_TBR_drives_release"]["value"]
        tbr_val = case["tbr"]["modelled_mean" if which == "modelled" else "measured"]["value"]
        gamma_val = irrad["neutron_rate"]["value"]
        k_top_val = rel["k_top_m_s"]["value"]
        k_wall_val = rel["k_wall_m_s"]["value"]

        windows = [
            IrradiationWindow(w["start"]["value"], w["stop"]["value"], gamma_val)
            for w in irrad["windows_s"]
        ]

        half_life_y = rel.get("tritium_half_life_y", {}).get("value", 12.32)
        lam = math.log(2.0) / (half_life_y * 365.25 * 86400.0)

        return cls(
            volume=volume_m3,
            A_top=a_top_m2,
            k_top=k_top_val,
            A_wall=a_wall_m2,
            k_wall=k_wall_val,
            TBR=tbr_val,
            windows=windows,
            decay=decay,
            decay_constant=lam,
        )

    def _intervals(self, t_max: float) -> list[tuple[float, float, float]]:
        """Partition [0, t_max] into intervals with constant source rate (particles/s)."""
        edges = {0.0, float(t_max)}
        for w in self.windows:
            if w.start < t_max:
                edges.add(max(0.0, w.start))
            if w.stop < t_max:
                edges.add(max(0.0, w.stop))
        sorted_edges = sorted(edges)
        intervals = []
        for a, b in zip(sorted_edges[:-1], sorted_edges[1:]):
            mid = 0.5 * (a + b)
            s_rate = 0.0
            for w in self.windows:
                if w.start <= mid < w.stop:
                    s_rate += self.TBR * w.neutron_rate
            intervals.append((a, b, s_rate))
        return intervals

    def solve(self, times: Sequence[float]) -> dict[str, np.ndarray]:
        """Compute exact state at the given evaluation times (s).

        Returns dict containing:
            times: evaluation times (s)
            concentration: c(t) (particles/m3)
            inventory: V * c(t) (particles)
            top_flux: k_top * c(t) (particles/m2/s)
            cumulative_top_release: particles released from top
            cumulative_top_release_per_area: particles/m2 released from top
            cumulative_wall_release: particles lost to wall
            cumulative_decay: particles decayed
            cumulative_source: total particles created up to t
        """
        t_arr = np.asarray(times, dtype=float)
        if np.any(t_arr < 0):
            raise ValueError("evaluation times must be non-negative")
        if len(t_arr) == 0:
            return {}

        t_max = float(np.max(t_arr))
        intervals = self._intervals(t_max)

        # March through intervals and compute exact integrated values
        alpha = self.alpha
        inv_alpha = 1.0 / alpha if alpha > 0 else 0.0

        c_state = 0.0
        q_top_state = 0.0
        q_wall_state = 0.0
        q_decay_state = 0.0
        src_state = 0.0

        # We will record values at each target time
        c_out = np.zeros_like(t_arr)
        q_top_out = np.zeros_like(t_arr)
        q_wall_out = np.zeros_like(t_arr)
        q_decay_out = np.zeros_like(t_arr)
        src_out = np.zeros_like(t_arr)

        for a, b, s_val in intervals:
            dt_full = b - a
            s_vol = s_val / self.V

            mask = (t_arr >= a) & (t_arr <= b)
            if np.any(mask):
                dt_eval = t_arr[mask] - a
                if alpha > 0:
                    e_eval = np.exp(-alpha * dt_eval)
                    c_eval = c_state * e_eval + (s_vol * inv_alpha) * (1.0 - e_eval)
                    int_c_eval = c_state * inv_alpha * (1.0 - e_eval) + (s_vol * inv_alpha) * (
                        dt_eval - inv_alpha * (1.0 - e_eval)
                    )
                else:
                    c_eval = c_state + s_vol * dt_eval
                    int_c_eval = c_state * dt_eval + 0.5 * s_vol * dt_eval**2

                c_out[mask] = c_eval
                q_top_out[mask] = q_top_state + self.A_top * self.k_top * int_c_eval
                q_wall_out[mask] = q_wall_state + self.A_wall * self.k_wall * int_c_eval
                q_decay_out[mask] = q_decay_state + self.decay_constant * self.V * int_c_eval
                src_out[mask] = src_state + s_val * dt_eval

            # Advance state across the full interval [a, b]
            if alpha > 0:
                e_full = math.exp(-alpha * dt_full)
                int_c_full = c_state * inv_alpha * (1.0 - e_full) + (s_vol * inv_alpha) * (
                    dt_full - inv_alpha * (1.0 - e_full)
                )
                c_state = c_state * e_full + (s_vol * inv_alpha) * (1.0 - e_full)
            else:
                int_c_full = c_state * dt_full + 0.5 * s_vol * dt_full**2
                c_state = c_state + s_vol * dt_full

            q_top_state += self.A_top * self.k_top * int_c_full
            q_wall_state += self.A_wall * self.k_wall * int_c_full
            q_decay_state += self.decay_constant * self.V * int_c_full
            src_state += s_val * dt_full

        inv_out = self.V * c_out
        top_flux_out = self.k_top * c_out
        per_area_out = q_top_out / self.A_top if self.A_top > 0 else np.zeros_like(q_top_out)

        return {
            "times": t_arr,
            "concentration": c_out,
            "inventory": inv_out,
            "top_flux": top_flux_out,
            "cumulative_top_release": q_top_out,
            "cumulative_top_release_per_area": per_area_out,
            "cumulative_wall_release": q_wall_out,
            "cumulative_decay": q_decay_out,
            "cumulative_source": src_out,
        }

    def solve_ode_numerical(self, times: Sequence[float]) -> dict[str, np.ndarray]:
        """Numerical solution via solve_ivp for verification."""
        t_arr = np.asarray(times, dtype=float)
        t_max = float(np.max(t_arr))
        intervals = self._intervals(t_max)

        y0 = [0.0, 0.0, 0.0, 0.0]  # c, q_top, q_wall, q_decay
        t_all = [0.0]
        y_all = [y0]

        for a, b, s_val in intervals:
            s_vol = s_val / self.V

            def rhs(_t, y):
                c, _q_t, _q_w, _q_d = y
                dc_dt = s_vol - self.alpha * c
                dq_top = self.A_top * self.k_top * c
                dq_wall = self.A_wall * self.k_wall * c
                dq_dec = self.decay_constant * self.V * c
                return [dc_dt, dq_top, dq_wall, dq_dec]

            sub_t = t_arr[(t_arr >= a) & (t_arr <= b)]
            eval_pts = sorted(set(sub_t) | {b})
            sol = solve_ivp(rhs, (a, b), y0, method="DOP853", t_eval=eval_pts, rtol=1e-11, atol=1e-12)
            y0 = [sol.y[0, -1], sol.y[1, -1], sol.y[2, -1], sol.y[3, -1]]
            for i, tt in enumerate(sol.t):
                t_all.append(tt)
                y_all.append([sol.y[0, i], sol.y[1, i], sol.y[2, i], sol.y[3, i]])

        t_nodes = np.array(t_all)
        y_nodes = np.array(y_all)

        c_interp = np.interp(t_arr, t_nodes, y_nodes[:, 0])
        q_top_interp = np.interp(t_arr, t_nodes, y_nodes[:, 1])
        q_wall_interp = np.interp(t_arr, t_nodes, y_nodes[:, 2])
        q_dec_interp = np.interp(t_arr, t_nodes, y_nodes[:, 3])

        return {
            "times": t_arr,
            "concentration": c_interp,
            "inventory": self.V * c_interp,
            "cumulative_top_release": q_top_interp,
            "cumulative_wall_release": q_wall_interp,
            "cumulative_decay": q_dec_interp,
        }


def check_self_consistency():
    """Verify R0 exact solution against numerical solve and mass conservation."""
    import yaml

    repo = Path(__file__).resolve().parents[1]
    case_path = repo / "case.yaml"
    case = yaml.safe_load(case_path.read_text(encoding="utf-8"))

    sample_days = np.array(case["sampling"]["IV"]["times_day"]["value"])
    sample_s = sample_days * 86400.0

    # Test with decay enabled
    m_decay = R0Model.from_case(case, decay=True)
    res_exact = m_decay.solve(sample_s)
    res_num = m_decay.solve_ode_numerical(sample_s)

    # Mass balance: Source = Inventory + Q_top + Q_wall + Q_decay
    src = res_exact["cumulative_source"]
    acc = (
        res_exact["inventory"]
        + res_exact["cumulative_top_release"]
        + res_exact["cumulative_wall_release"]
        + res_exact["cumulative_decay"]
    )
    mb_err = float(np.max(np.abs(src - acc) / np.maximum(src, 1e-12)))

    # Difference between exact and numerical
    q_diff = float(
        np.max(
            np.abs(res_exact["cumulative_top_release"] - res_num["cumulative_top_release"])
            / np.maximum(res_exact["cumulative_top_release"], 1e-12)
        )
    )

    # Test with decay disabled
    m_nodecay = R0Model.from_case(case, decay=False)
    res_nodecay = m_nodecay.solve(sample_s)
    src_nd = res_nodecay["cumulative_source"]
    acc_nd = (
        res_nodecay["inventory"]
        + res_nodecay["cumulative_top_release"]
        + res_nodecay["cumulative_wall_release"]
    )
    mb_err_nd = float(np.max(np.abs(src_nd - acc_nd) / np.maximum(src_nd, 1e-12)))

    return {
        "mass_balance_rel_err_with_decay": mb_err,
        "exact_vs_numerical_q_top_rel_err": q_diff,
        "mass_balance_rel_err_without_decay": mb_err_nd,
    }


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="R0 model self-checks")
    parser.add_argument("--check", action="store_true", help="Assert self-checks pass")
    args = parser.parse_args()

    figures = check_self_consistency()
    print("R0 model self-checks:")
    for k, v in figures.items():
        print(f"  {k}: {v:.2e}")

    if args.check:
        if figures["mass_balance_rel_err_with_decay"] > 1e-10:
            print("R0 check failed: mass balance with decay exceeded 1e-10", file=sys.stderr)
            sys.exit(1)
        if figures["exact_vs_numerical_q_top_rel_err"] > 1e-8:
            print("R0 check failed: exact vs numerical exceeded 1e-8", file=sys.stderr)
            sys.exit(1)
        print("R0 check: all self-checks passed")
