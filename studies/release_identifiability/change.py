#!/usr/bin/env python3
"""Run 2: does the release rate change during the run?

The BABY-1L paper (arXiv 2509.26174, sec. 3.3) suggests the salt freezing
mid-run as a possible explanation for run 2's slower release. The freeze time
is not in the run's files. This script asks a simpler question of the data: fit
the 0D model (TBR at the lab's value) with k_top allowed to take a second value
from one of the sample times CHANGE_DAYS on, and compare with one k_top for the
whole run. The change days were chosen after looking at the data, so the best
of them is calibrated by a parametric bootstrap from the one-k_top fit (the max
improvement over the same days in each simulated data set).

Writes change.json. Run with the study venv:
  .venv/bin/python change.py --reps 200 --workers 10
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

import fitlib as F  # noqa: E402
import model0d as M  # noqa: E402

OUT = HERE / "change.json"
RESULTS = HERE / "results.json"
DAY = 86400.0
# the 6th to 10th IV sample times (days ~9 to ~27); the second k_top applies from each
CHANGE_DAYS = tuple(F.INC["runs"]["2"]["streams"]["IV"]["times_day"][5:10])


def ds_one():
    return F.Dataset("run2", 2, None, False, streams=("IV", "OV"), tbr="lab", noise="fitted")


def ds_step(day):
    ds = F.Dataset("run2", 2, None, True, streams=("IV", "OV"), tbr="lab", noise="fitted")
    ds.switch = day * DAY
    return ds


def _fit(args):
    day, y = args
    ds = ds_one() if day is None else ds_step(day)
    th, v = F.fit(ds, y=None if y is None else np.array(y))
    return day, v, th.tolist()


def alpha_of(kt, kw):
    return (M.A_TOP * kt + M.A_WALL * kw) / M.V_M3


def summarise(day, v, th, v0):
    ds = ds_step(day)
    ktops, kw, _, f = ds.physical(np.array(th))
    prof = F.profile(ds, np.array(th), v, "log10_ktop_p2")
    c = prof["ci95"]
    return {"change_day": day, "two_delta_nll_gain": 2 * (v0 - v),
            "k_top_before_m_s": ktops[0], "k_top_after_m_s": ktops[1], "k_wall_m_s": kw,
            "tau_before_day": 1 / alpha_of(ktops[0], kw) / DAY, "tau_after_day": 1 / alpha_of(ktops[1], kw) / DAY,
            "k_top_after_ci95_m_s": [10 ** c["lo"], 10 ** c["hi"]] if c["status"] == "identifiable" else None,
            "k_top_after_status": c["status"],
            "floors_Bq": {k: float(x) for k, x in f.items() if np.isfinite(x)}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--seed", type=int, default=20261009)
    a = ap.parse_args()
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, mp_context=mp.get_context("spawn")) as pool:
        fits = list(pool.map(_fit, [(None, None)] + [(d, None) for d in CHANGE_DAYS]))
        _, v0, th0 = fits[0]
        per_day = [summarise(d, v, th, v0) for d, v, th in fits[1:]]
        best = max(per_day, key=lambda r: r["two_delta_nll_gain"])
        # null: data from the one-k_top fit, refitted with and without each change day
        ds0 = ds_one()
        rng = np.random.default_rng(a.seed)
        ys = [rng.multivariate_normal(ds0.predict(np.array(th0)), ds0.sigma(np.array(th0))).tolist()
              for _ in range(a.reps)]
        res = list(pool.map(_fit, [(d, y) for y in ys for d in (None, *CHANGE_DAYS)], chunksize=6))
    n = len(CHANGE_DAYS) + 1
    gains = [2 * (res[i * n][1] - min(r[1] for r in res[i * n + 1:(i + 1) * n])) for i in range(a.reps)]
    g = np.array(gains)
    run1 = json.loads(RESULTS.read_text())["fits"]["run1"]["labTBR"]
    out = {"_generated": "Written by change.py. Do not edit by hand.",
           "one_k_top": {"nll": v0, "k_top_m_s": ds0.physical(np.array(th0))[0][0]},
           "per_change_day": per_day, "best": best,
           "null": {"reps": a.reps, "seed": a.seed, "gain_quantiles": {q: float(np.quantile(g, q)) for q in (0.5, 0.95)},
                    "p_value_best": (1 + int(np.sum(g >= best["two_delta_nll_gain"]))) / (1 + a.reps)},
           "run1_labTBR_k_top_m_s": run1["physical"]["k_top_m_s"][0],
           "provenance": {"increments_sha256": hashlib.sha256((HERE / "increments.json").read_bytes()).hexdigest(),
                          "results_fits_run1_sha256": hashlib.sha256(json.dumps(
                              json.loads(RESULTS.read_text())["fits"]["run1"], sort_keys=True).encode()).hexdigest(),
                          "reps": a.reps, "seconds": time.time() - t0}}
    OUT.write_text(json.dumps(out, indent=1, default=float) + "\n")
    for r in per_day:
        print(f"change at day {r['change_day']}: gain {r['two_delta_nll_gain']:.2f}, tau {r['tau_before_day']:.1f} -> "
              f"{r['tau_after_day']:.1f} d, k_top {r['k_top_before_m_s']:.2e} -> {r['k_top_after_m_s']:.2e} "
              f"({r['k_top_after_status']})")
    print(f"best day {best['change_day']} gain {best['two_delta_nll_gain']:.2f}, null 95% {out['null']['gain_quantiles'][0.95]:.2f}, "
          f"p {out['null']['p_value_best']:.3f}; wrote {OUT} in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
