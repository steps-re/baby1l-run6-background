#!/usr/bin/env python3
"""Is salt diffusivity identifiable from the BABY-1L release curves?

For each run window, fix D on a grid, fit everything else (k_top, k_wall,
TBR scale, noise floors) with the 1D finite-volume column (fv1d.py, checked
against FESTIM in festim_check.py), and record the profile likelihood of D.
The large-D end is the well-mixed 0D model.

Reference diffusivities (study brief values, not re-derived here):
  Calderoni et al. 2008, T in FLiBe: D = 9.3e-7 exp(-42 kJ/mol / RT) m2/s,
  evaluated at each run's logged minimum salt temperature;
  H2 in FLiBe near 900 K: 7e-10 m2/s (Nakamura, Fukada and Nishiumi, J. Plasma
  Fusion Res. SERIES 11, 25 (2015), eq. 10: 2.09e-8 exp(-25.2 kJ/mol / RT) gives
  7.2e-10 at 900 K, data taken at 773-873 K).

Also runs a synthetic control on run 1's design: data generated well mixed
must give a profile that is flat above some D (a lower bound only), and data
generated at a slow D must be told apart from the well-mixed case. Two slow
truths: D = 7e-10 m2/s, and the Calderoni value at run 1's salt temperature
(the power that matters for reading "Calderoni not rejected").

Usage (two parts):
  run/refit/.venv/bin/python studies/release_identifiability/diffusivity.py --part windows --workers 8
  run/refit/.venv/bin/python studies/release_identifiability/diffusivity.py --part controls --control-reps 4 --workers 8
High statistics belong on a many-core machine (`make reproduce-full`). Writes diffusivity.json with provenance.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

import multiprocessing as _mp

if _mp.parent_process() is not None:  # spawned worker: single-threaded BLAS, set before numpy loads
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "NUMEXPR_NUM_THREADS"):
        os.environ[_v] = "1"

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fitlib as F  # noqa: E402
import fv1d  # noqa: E402
from analysis import WINDOWS, ds_for, make_pool, provenance  # noqa: E402
from fv1d import d_calderoni_T  # noqa: E402

D_GRID = np.round(np.arange(-10.0, -4.99, 0.25), 3)  # log10 m2/s
D_H2 = 7.0e-10
N_CELLS = 60
CONTROL_REPS = 4


class DDataset(F.Dataset):
    D = 1e-3

    def predict(self, theta, par="A"):
        ktops, kw, s, _ = self.physical(theta, par)
        src = s * self.tbr_meas * self.gamma
        top, wall = fv1d.cumulative_release(self.D, ktops, kw, src, self.irr, self.switch,
                                            self.times.get("IV", []), self.times.get("OV", []), n=N_CELLS)
        out = []
        for st, cum in (("IV", top), ("OV", wall)):
            if st in self.streams:
                out.append(np.diff(np.concatenate([[0.0], cum])) * F.BQ_PER_PARTICLE)
        return np.concatenate(out)


def make(w, tbr="free"):
    base = ds_for(w, "freeTBR")
    d = DDataset(base.name, base.run, base.t_max_day, base.step, streams=base.streams, tbr=tbr)
    return d


def d_profile(d, y=None, th0=None):
    """Profile NLL over D_GRID plus the well-mixed end (D = 1e-3 m2/s).

    One optimiser setting everywhere (L-BFGS-B, no Powell polish). Two warm-
    started sweeps, large D to small and back, and the better fit is kept at
    every D, the well-mixed end included. So no point, the well-mixed end in
    particular, gets a more thorough search than the others.
    """
    polish = F.POLISH["on"]
    F.POLISH["on"] = False
    d.D = 1e-3
    first = F.fit(d, y=y) if th0 is None else F.fit(d, y=y, starts=[th0] + d.starts()[2:5])
    pts = [1e-3] + [10 ** g for g in sorted(D_GRID, reverse=True)]
    best = {pts[0]: first}
    for order in (pts, pts[::-1]):
        prev = best[order[0]][0] if order[0] in best else first[0]
        for D in order:
            d.D = D
            starts = [prev, first[0], d.starts()[2], d.starts()[4]]
            th, v = F.fit(d, y=y, starts=starts)
            if D not in best or v < best[D][1]:
                best[D] = (th, v)
            prev = best[D][0]
    F.POLISH["on"] = polish
    th_inf, v_inf = best[1e-3]
    prof = {float(np.round(math.log10(D), 3)): (best[D][1], best[D][0]) for D in pts[1:]}
    vmin = min(v for _, v in best.values())
    grid = sorted(prof)
    d2 = np.array([2 * (prof[g][0] - vmin) for g in grid])
    return grid, d2, 2 * (v_inf - vmin), prof, th_inf


def bounds_from(grid, d2, d2_inf):
    """95% region in log10 D (the well-mixed limit counts as the right end of
    the grid). If the 95% set is disjoint, the status and edges describe the
    whole set, and the grid points outside the piece around the best point are
    listed as islands."""
    g = np.array(grid)
    d = np.append(np.array(d2), d2_inf)  # index len(g) = well-mixed
    ok = d < F.CHI2_95
    i0 = int(np.argmin(d))
    lo_i = i0
    while lo_i > 0 and ok[lo_i - 1]:
        lo_i -= 1
    hi_i = i0
    while hi_i < len(d) - 1 and ok[hi_i + 1]:
        hi_i += 1
    lo = None if lo_i == 0 else float(np.interp(F.CHI2_95, [d[lo_i], d[lo_i - 1]], [g[lo_i], g[lo_i - 1]]))
    if hi_i == len(d) - 1:
        hi = None
    elif hi_i + 1 < len(g):
        hi = float(np.interp(F.CHI2_95, [d[hi_i], d[hi_i + 1]], [g[hi_i], g[hi_i + 1]]))
    else:
        hi = float(g[-1])
    islands = [float(g[k]) for k in range(len(g)) if ok[k] and not (lo_i <= k <= hi_i)]
    if islands:
        # a disjoint 95% set: report the extent of the whole set, not the piece around the best point
        ks = np.where(ok)[0]
        lo = None if ks[0] == 0 else float(np.interp(F.CHI2_95, [d[ks[0]], d[ks[0] - 1]], [g[ks[0]], g[ks[0] - 1]]))
        if ks[-1] == len(d) - 1:
            hi = None
        elif ks[-1] + 1 < len(g):
            hi = float(np.interp(F.CHI2_95, [d[ks[-1]], d[ks[-1] + 1]], [g[ks[-1]], g[ks[-1] + 1]]))
        else:
            hi = float(g[-1])
    if lo is None and hi is None:
        status = "not identifiable"
    elif hi is None:
        status = "lower bound only"
    elif lo is None:
        status = "upper bound only"
    else:
        status = "identifiable"
    return {"lo_log10": lo, "hi_log10": hi, "status": status, "islands_log10": islands}


def _window_task(args):
    wname, tbr = args
    F.POLISH["on"] = False  # same optimiser at every D, the well-mixed end included
    w = {x[0]: x for x in WINDOWS}[wname]
    d = make(w, tbr)
    T_K = d.meta["temperature_salt_C"] + 273.15
    grid, d2, d2_inf, prof, _ = d_profile(d)
    D_cal = d_calderoni_T(T_K)
    at = lambda D: float(np.interp(math.log10(D), grid, d2))  # noqa: E731
    kbest = grid[int(np.argmin(d2))]
    return wname, tbr, {
        "T_K": T_K, "D_calderoni_m2_s": D_cal,
        "two_delta_nll": d2.tolist(), "two_delta_nll_well_mixed": d2_inf,
        "at_calderoni": at(D_cal), "at_H2_7e-10": at(D_H2),
        "best_log10_D": float(kbest) if np.min(d2) <= d2_inf else None,
        "k_top_at_calderoni_m_s": 10 ** float(np.interp(math.log10(D_cal), grid,
                                                        [prof[g][1][0] for g in grid])),
        "region95": bounds_from(grid, d2, d2_inf),
        "fit_along_profile": {str(g): dict(zip(d.names_A, [float(x) for x in prof[g][1]])) for g in grid},
    }


def _control_task(args):
    """Synthetic data on run 1's design at a known D; profile must cover it."""
    D_true, seed = args
    F.POLISH["on"] = False
    d = make(WINDOWS[0])
    d.D = D_true
    th_t, _ = F.fit(d)  # parameters that best reproduce run 1 at this D
    y = np.random.default_rng(seed).multivariate_normal(d.predict(th_t), d.sigma(th_t))
    grid, d2, d2_inf, _, _ = d_profile(d, y=y, th0=th_t)
    if D_true >= 10 ** grid[-1]:
        inside = d2_inf < F.CHI2_95
    else:
        inside = float(np.interp(math.log10(D_true), grid, d2)) < F.CHI2_95
    return {"D_true": D_true, "region95": bounds_from(grid, d2, d2_inf), "truth_inside_95": bool(inside),
            "well_mixed_rejected": bool(d2_inf >= F.CHI2_95),
            "two_delta_nll": d2.tolist(), "two_delta_nll_well_mixed": d2_inf}


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["windows", "controls"], required=True,
                    help="run in two parts so neither passes ~5 min on a laptop")
    ap.add_argument("--control-reps", type=int, default=CONTROL_REPS,
                    help="synthetic replicates per truth (three truths)")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    t0 = time.time()
    out_path = HERE / "diffusivity.json"
    res = json.loads(out_path.read_text()) if out_path.exists() else {}
    res.update({"_generated": "Written by diffusivity.py. Do not edit by hand.",
                "grid_log10_D": D_GRID.tolist(), "n_cells": N_CELLS,
                "fv_large_D_vs_0D_max_rel_gap": fv1d.check_large_D(),
                "fv_grid_n60_vs_n240_at_7e-10": fv1d.check_grid(7e-10),
                "fv_grid_n60_vs_n240_at_1e-10": fv1d.check_grid(1e-10),
                "D_H2_ref_m2_s": D_H2})
    pool = make_pool(args.workers)
    part_args = {"workers": args.workers} | ({"control_reps": args.control_reps} if args.part == "controls" else {})
    res.setdefault("provenance", {})[args.part] = provenance(part_args)
    if args.part == "windows":
        res["runs"] = {"free": {}, "fixed": {}}
        tasks = [(w[0], tbr) for tbr in ("free", "fixed") for w in WINDOWS]
        for wname, tbr, r in pool.map(_window_task, tasks):
            res["runs"][tbr][wname] = r
            print(f"{wname:14s} TBR {tbr:5s} 2dNLL at Calderoni {r['at_calderoni']:.1f}, at 7e-10 "
                  f"{r['at_H2_7e-10']:.1f}, well-mixed {r['two_delta_nll_well_mixed']:.1f}, "
                  f"region {r['region95']}", flush=True)
        res["runtime_windows_s"] = time.time() - t0
    else:
        n = args.control_reps
        T1 = make(WINDOWS[0]).meta["temperature_salt_C"] + 273.15
        D_cal1 = d_calderoni_T(T1)
        truths = (("well_mixed_truth", 1e-3), ("D_7e-10_truth", D_H2), ("D_calderoni_run1_T_truth", D_cal1))
        # seeds: the first 2n keep the original two truths' streams, the Calderoni truth gets new ones
        seeds = [int(x) for x in np.random.SeedSequence(20260930).generate_state(2 * n)]
        seeds += [int(x) for x in np.random.SeedSequence(20260931).generate_state(n)]
        tasks = [(D, sd) for k, (_, D) in enumerate(truths) for sd in seeds[k * n:(k + 1) * n]]
        rows = list(pool.map(_control_task, tasks))
        ctrl = {"calderoni_T_K": T1}
        for label, D_true in truths:
            rr = [r for r in rows if r["D_true"] == D_true]
            ctrl[label] = {"D_true": D_true, "reps": len(rr),
                           "truth_inside_95_fraction": float(np.mean([r["truth_inside_95"] for r in rr])),
                           "well_mixed_rejected_fraction": float(np.mean([r["well_mixed_rejected"] for r in rr])),
                           "statuses": [r["region95"]["status"] for r in rr],
                           "lower_bounds_log10": [r["region95"]["lo_log10"] for r in rr],
                           "example": rr[0]}
            print("control", label, ctrl[label]["truth_inside_95_fraction"],
                  ctrl[label]["well_mixed_rejected_fraction"], ctrl[label]["statuses"], flush=True)
        res["synthetic_control_run1_design"] = ctrl
        res["runtime_controls_s"] = time.time() - t0
    pool.shutdown()
    out_path.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(f"done in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
