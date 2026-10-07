#!/usr/bin/env python3
"""Identifiability and transferability of the 0D release model on BABY-1L.

Steps (numbers land in results.json; README.md is generated from it):
  1. Fit each run window (IV + OV increments) with the production term fixed
     at the lab's measured TBR, then again with TBR free as a nuisance scale,
     then IV only with TBR free.
  2. Profile likelihood for k_top, k_wall, alpha and TBR scale; observed
     Fisher information at the optimum, compared with the profiles.
  3. Transferability: joint fit with one shared k_top across the He runs,
     likelihood ratio against separate fits, p-value from a parametric
     bootstrap under the shared-k null.
  4. Synthetic controls: data generated from known parameters with the same
     sampling times and noise; recovery, 95% profile-CI coverage, a flat
     direction that must stay flat, and false-positive rate and power of the
     transferability test.

Usage (three parts; run fits first, the other two read its fits):
  run/refit/.venv/bin/python studies/release_identifiability/analysis.py --part fits --workers 8
  run/refit/.venv/bin/python studies/release_identifiability/analysis.py --part controls --reps 10 --workers 8
  run/refit/.venv/bin/python studies/release_identifiability/analysis.py --part transfer --boot-reps 20 --workers 8
High statistics belong on a many-core machine (`make reproduce-full`). Every part records its provenance
(git sha, time, arguments) in results.json; controls and transfer also record a
hash of the fits they started from, and make_readme.py refuses mismatched parts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# Worker processes (spawned) must run single-threaded BLAS, or 32 workers each
# start 32 BLAS threads. This runs before numpy is imported in every worker.
BLAS_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
if mp.parent_process() is not None:
    for _v in BLAS_VARS:
        os.environ[_v] = "1"

import numpy as np  # noqa: E402
from scipy.stats import chi2 as chi2dist  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fitlib as F  # noqa: E402
import model0d as M  # noqa: E402

OUT = HERE / "results.json"
CASE_YAML = HERE.parents[1] / "case.yaml"


def make_pool(workers):
    """Process pool whose workers run single-threaded BLAS (env inherited at spawn)."""
    for v in BLAS_VARS:
        os.environ[v] = "1"
    return ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"))


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fits_hash(fits):
    """Hash of the stored fits exactly as results.json holds them (after a JSON round trip)."""
    canon = json.dumps(json.loads(json.dumps(fits, default=float)), sort_keys=True)
    return hashlib.sha256(canon.encode()).hexdigest()


def provenance(arguments):
    def git(*a):
        try:
            return subprocess.run(["git", *a], cwd=HERE, capture_output=True, text=True, check=True).stdout.strip()
        except Exception:  # noqa: BLE001
            return None
    dirty = git("status", "--porcelain", "--", "*.py")
    return {"git_sha": git("rev-parse", "HEAD"), "code_dirty": bool(dirty) if dirty is not None else None,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "arguments": arguments, "cpu_count": os.cpu_count(),
            "increments_sha256": file_sha256(HERE / "increments.json"),
            "quench_sha256": file_sha256(HERE / "quench.json") if (HERE / "quench.json").exists() else None}
SW3 = F.INC["runs"]["3"]["gas_switch_s"] / 86400.0
SW4 = F.INC["runs"]["4"]["gas_switch_s"] / 86400.0

# name, run, keep samples up to (day), k_top steps at the gas switch, gas label
WINDOWS = [
    ("run1", 1, None, False, "He"),
    ("run2", 2, None, False, "He"),
    ("run3_He", 3, SW3, False, "He (before the switch to 1000 ppm H2)"),
    ("run4_1000ppm", 4, SW4, False, "1000 ppm H2 (before the switch to 3.5% H2)"),
    ("run6", 6, None, True, "3.5% H2, then He (k_top steps at the switch)"),
]
SECONDARY = [
    ("run3_full_step", 3, None, True, "He then 1000 ppm H2, one k_top per gas"),
    ("run4_full_step", 4, None, True, "1000 ppm then 3.5% H2, one k_top per gas"),
]
VARIANTS = {
    "fixedTBR": dict(streams=("IV", "OV"), tbr="fixed"),
    "freeTBR": dict(streams=("IV", "OV"), tbr="free"),
    "labTBR": dict(streams=("IV", "OV"), tbr="lab"),
    "IVonly_freeTBR": dict(streams=("IV",), tbr="free"),
}


# windows whose runs have a quench-matched background (quench.py)
QUENCH_WINDOWS = [w for w in WINDOWS if w[1] in (3, 4)]


def ds_for(w, variant, noise="fitted", background="published"):
    name, run, tmax, step, _ = w
    return F.Dataset(name, run, tmax, step, noise=noise, background=background, **VARIANTS[variant])


def phys_summary(ds, th):
    ktops, kw, s, sig = ds.physical(th)
    alpha = (M.A_TOP * ktops[0] + M.A_WALL * kw) / M.V_M3
    return {"k_top_m_s": ktops, "k_wall_m_s": kw, "tbr_scale": s,
            "TBR_used": s * ds.tbr_meas, "alpha_per_s": alpha, "tau_day": 1 / alpha / 86400.0,
            "top_share_phi": M.A_TOP * ktops[0] / (M.A_TOP * ktops[0] + M.A_WALL * kw),
            "noise_floor_Bq": {k: v for k, v in sig.items() if np.isfinite(v)}}


def profile_names(ds):
    names = ["log10_ktop"] + (["log10_ktop_p2"] if ds.step else []) + ["log10_kwall", "log10_alpha"]
    if ds.tbr == "free":
        names.append("log10_tbr_scale")
    return names


def fit_block(ds, do_profiles=True):
    t0 = time.time()
    th, v = F.fit(ds)
    out = {"n_points": ds.n(), "nll": v, "theta": dict(zip(ds.names_A, th.tolist())),
           "physical": phys_summary(ds, th),
           "chi2_counting_only_at_mle": ds.chi2(th),
           "floors_at_lower_bound": [n.replace("log10_sig_", "") for n, x in zip(ds.names_A, th)
                                     if n.startswith("log10_sig_") and x - F.BOUNDS[n][0] < 1e-3]}
    if do_profiles:
        out["profiles"] = {n: F.profile(ds, th, v, n) for n in profile_names(ds)}
        out["fisher"] = F.fisher_summary(ds, th)
        out["wald_vs_profile"] = wald_vs_profile(out)
    out["seconds"] = time.time() - t0
    return th, v, out


def wald_vs_profile(block):
    fs = block["fisher"]
    res = {}
    for n, p in block["profiles"].items():
        if n not in fs["names"]:
            continue
        w = fs["wald95_halfwidth_log10"][fs["names"].index(n)]
        ci = p["ci95"]
        if ci["lo"] is not None and ci["hi"] is not None:
            prof_hw = 0.5 * (ci["hi"] - ci["lo"])
            res[n] = {"wald95_halfwidth": w, "profile95_halfwidth": prof_hw, "ratio": w / prof_hw}
        else:
            res[n] = {"wald95_halfwidth": w, "profile95_halfwidth": None, "ratio": None,
                      "note": f"profile {ci['status']}; Wald interval is not meaningful here"}
    return res


def _fit_task(args):
    wname, var = args
    background = "quench" if var.endswith("_quench") else "published"
    saved = F.BOUNDS["log10_tbr_scale"]
    try:
        if var.endswith("_wideTBR"):
            F.BOUNDS["log10_tbr_scale"] = (-1.0, 1.0)  # sensitivity: TBR within a factor of 10
        ds = ds_for(WIN[wname], var.replace("_quench", "").replace("_wideTBR", ""), background=background)
        th, v, blk = fit_block(ds)
    finally:
        F.BOUNDS["log10_tbr_scale"] = saved
    return wname, var, th.tolist(), blk


def gof_counting(w):
    """Counting-noise-only fit: is counting noise enough to explain the scatter?"""
    ds = ds_for(w, "freeTBR", noise="counting")
    th, v = F.fit(ds)
    c2 = ds.chi2(th)
    dof = ds.n() - len(ds.names_A)
    return {"chi2": c2, "dof": dof, "chi2_per_dof": c2 / dof, "p_value": float(chi2dist.sf(c2, dof)),
            "k_top_m_s": ds.physical(th)[0], "k_wall_m_s": ds.physical(th)[1]}


# ---------------------------------------------------------------- synthetic
def simulate(ds, th, rng):
    mu = ds.predict(th)
    S = ds.sigma(th)
    return rng.multivariate_normal(mu, S)


WIN = {w[0]: w for w in WINDOWS + SECONDARY}


def _cov_one(args):
    """One synthetic replicate: fit, then LR at the truth for each parameter.
    Same optimiser as the real fits: the full 8-point start grid plus Powell polish."""
    wname, var, truth, seed, names = args
    ds = ds_for(WIN[wname], var)
    rng = np.random.default_rng(seed)
    y = simulate(ds, truth, rng)
    th, v = F.fit(ds, y=y)
    thB, tB = ds.to_B(th), ds.to_B(truth)
    out = {}
    for n in names:
        par = "B" if n == "log10_alpha" else "A"
        tv = (truth if par == "A" else tB)[ds.names(par).index(n)]
        if n == "log10_kwall" and tv <= F.BOUNDS[n][0] + 1e-3:
            continue
        err = (th if par == "A" else thB)[ds.names(par).index(n)] - tv
        _, vc = F.constrained_fit(ds, par, n, tv, th if par == "A" else thB, y)
        out[n] = (float(err), bool(2 * (vc - min(v, vc)) < F.CHI2_95))
    floors = {n.replace("log10_sig_", ""): bool(x - F.BOUNDS[n][0] < 1e-3)
              for n, x in zip(ds.names_A, th) if n.startswith("log10_sig_")}
    return {"par": out, "floor_at_lower_bound": floors}


def coverage(pool, wname, var, truth, seeds, names=("log10_ktop", "log10_kwall", "log10_alpha")):
    raw = list(pool.map(_cov_one, [(wname, var, truth, sd, names) for sd in seeds]))
    rows = [r["par"] for r in raw]
    ds = ds_for(WIN[wname], var)
    res = {"truth_floors_Bq": {n.replace("log10_sig_", ""): 10 ** float(x) for n, x in zip(ds.names_A, truth)
                               if n.startswith("log10_sig_")},
           "fitted_floor_at_lower_bound_fraction": {s: float(np.mean([r["floor_at_lower_bound"][s] for r in raw]))
                                                    for s in raw[0]["floor_at_lower_bound"]}}
    for n in names:
        e = [r[n][0] for r in rows if n in r]
        i = [r[n][1] for r in rows if n in r]
        res[n] = {"reps": len(i), "coverage95": float(np.mean(i)) if i else None,
                  "median_error_log10": float(np.median(e)) if e else None,
                  "fraction_error_over_1_decade": float(np.mean(np.abs(e) > 1.0)) if e else None,
                  "spread_log10_p16_p84": [float(np.percentile(e, 16)), float(np.percentile(e, 84))] if e else None}
    return res


def _flat_one(args):
    """Generate IV+OV data from run 1's fitted truth, then fit IV only with TBR
    free. The k_top/k_wall split has no information there: the k_top profile
    must stay flat, while the same data with OV included must bound it."""
    thF, seed = args
    iv = ds_for(WINDOWS[0], "IVonly_freeTBR")
    full = ds_for(WINDOWS[0], "freeTBR")
    yF = simulate(full, thF, np.random.default_rng(seed))
    y_iv = yF[: iv.n()]
    th1, v1 = F.fit(iv, y=y_iv)
    p_iv = F.profile(iv, th1, v1, "log10_ktop", y=y_iv, half_width=2.0, n=21)
    th2, v2 = F.fit(full, y=yF)
    p_full = F.profile(full, th2, v2, "log10_ktop", y=yF, half_width=2.0, n=21)
    return {"IV_only": {"ci95": p_iv["ci95"], "max_2dnll_within_1_decade_of_truth": _max_within(p_iv, thF[0], 1.0)},
            "IV_plus_OV": {"ci95": p_full["ci95"]}}


def flat_direction_control(pool, thF, seeds):
    out = list(pool.map(_flat_one, [(thF, sd) for sd in seeds]))
    return {"truth_log10_ktop": float(thF[0]), "reps": out,
            "IV_only_not_identifiable_fraction": float(np.mean(
                [r["IV_only"]["ci95"]["status"] != "identifiable" for r in out])),
            "IV_plus_OV_identifiable_fraction": float(np.mean(
                [r["IV_plus_OV"]["ci95"]["status"] == "identifiable" for r in out]))}


def _max_within(p, center, half):
    g = np.array(p["grid"]); d = np.array(p["two_delta_nll"])  # noqa: E702
    m = np.abs(g - center) <= half
    return float(np.nanmax(d[m]))


# ---------------------------------------------------------------- transferability
def joint_nll(dss, seps, shared, values, ys=None):
    """Joint NLL with `shared` parameters fixed at `values` in every run and
    all other parameters free per run (a sum of per-run constrained fits)."""
    tot = 0.0
    for i, ds in enumerate(dss):
        th_i, v_i = seps[i]
        y = None if ys is None else ys[i]
        _, v = F.constrained_fit(ds, "A", shared, values, th_i, y)
        tot += v
    return tot


def lr_test(dss, shared, ys=None, starts=None):
    """LR of 'shared parameters equal across runs' against separate fits.

    The joint optimum is found over the shared parameters only (one or two
    dimensions), each evaluation a set of per-run constrained fits warm-started
    at that run's own optimum. Starts: the median of the per-run values, and
    each run's own value.
    """
    seps = [F.fit(ds, y=None if ys is None else ys[i], starts=None if starts is None else starts(ds))
            for i, ds in enumerate(dss)]
    per_run = np.array([[s[0][d.names_A.index(n)] for n in shared] for s, d in zip(seps, dss)])
    cands = [np.median(per_run, axis=0)] + list(per_run)
    fun = lambda v: joint_nll(dss, seps, shared, v, ys)  # noqa: E731
    best = None
    for x0 in cands:
        if len(shared) == 1:
            # the summed profile can be bimodal: scan a grid, then refine
            from scipy.optimize import minimize_scalar
            lo, hi = per_run.min() - 0.5, per_run.max() + 0.5
            g = np.linspace(lo, hi, 13)
            vals = [fun(np.array([x])) for x in g]
            k = int(np.argmin(vals))
            a_, b_ = g[max(k - 1, 0)], g[min(k + 1, len(g) - 1)]
            r = minimize_scalar(lambda x: fun(np.array([x])), bounds=(a_, b_), method="bounded",
                                options={"xatol": 1e-4})
            best = (np.array([r.x]), float(r.fun)) if r.fun < vals[k] else (np.array([g[k]]), float(vals[k]))
            break
        from scipy.optimize import minimize
        r = minimize(fun, x0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 400})
        if best is None or r.fun < best[1]:
            best = (r.x, float(r.fun))
    lr = 2 * (best[1] - sum(s[1] for s in seps))
    return max(lr, 0.0), best[0], seps


def _lr_one(args):
    names, shared, truths, seed = args
    dss = [ds_for(WIN[n], "freeTBR") for n in names]
    rng = np.random.default_rng(seed)
    ys = [simulate(ds, t, rng) for ds, t in zip(dss, truths)]
    return lr_test(dss, shared, ys)[0]


def null_truths(dss, shared, shared_vals, seps):
    """Per-run parameters under the null: shared values imposed, the rest refit."""
    out = []
    for ds, (th_i, _) in zip(dss, seps):
        th, _ = F.constrained_fit(ds, "A", shared, shared_vals, th_i)
        out.append(th)
    return out


def transferability(pool, group, shared, seeds, power_seeds, power_factors=(1.5, 2.0)):
    """Real and bootstrap LRs use the same optimiser settings as the real fits
    (the full 8-point start grid plus Powell polish), so the null distribution is
    calibrated for the statistic actually computed and no replicate gets a
    weaker search than the real data."""
    names = [w[0] for w in group]
    dss = [ds_for(w, "freeTBR") for w in group]
    lr, shared_vals, seps = lr_test(dss, shared)
    df = (len(dss) - 1) * len(shared)
    null_truth = null_truths(dss, shared, shared_vals, seps)
    null_lr = np.array(list(pool.map(_lr_one, [(names, shared, null_truth, sd) for sd in seeds])))
    crit = float(np.percentile(null_lr, 95))
    power = {}
    for fac in power_factors:
        alt = [t.copy() for t in null_truth]
        alt[1][dss[1].names_A.index(shared[0])] += math.log10(fac)
        alt_lr = np.array(list(pool.map(_lr_one, [(names, shared, alt, sd) for sd in power_seeds])))
        power[f"x{fac}"] = float(np.mean(alt_lr > crit))
    return {
        "runs": names, "shared": shared, "variant": "IV+OV, TBR free per run",
        "separate": {n: [float(s[0][d.names_A.index(n)]) for s, d in zip(seps, dss)] for n in shared},
        "shared_values": [float(x) for x in shared_vals],
        "LR": lr, "df": df,
        "separate_fit_theta": {n: [float(x) for x in s[0]] for n, s in zip(names, seps)}, "chi2_crit95": float(chi2dist.isf(0.05, df)),
        "p_asymptotic": float(chi2dist.sf(lr, df)),
        "null_LR_max": float(np.max(null_lr)),
        "null_degenerate": bool(np.median(null_lr) < 0.05),
        "p_bootstrap": float((1 + np.sum(null_lr >= lr)) / (1 + len(null_lr))),
        "bootstrap_reps": int(len(null_lr)),
        "null_LR_median": float(np.median(null_lr)),
        "null_false_positive_rate_at_asymptotic_crit": float(np.mean(null_lr > chi2dist.isf(0.05, df))),
        "bootstrap_crit95": crit,
        "power_if_second_run_differs_by": power,
        "power_parameter": shared[0],
        "power_reps": len(power_seeds),
    }


def ov_lag_test(w, floors, n_mc=20000, seed=11, background="published"):
    """Model-free test of the 0D model's OV shape.

    With constant k_top and k_wall, IV and OV are both integrals of the same
    c(t), so OV(t)/IV(t) is constant. Test at a common time t*: the first OV
    sample time (not the last) at which the IV has collected at least 80% of
    its window total. IV cumulative at t* is taken from the IV sample at the
    same time when one exists (within 0.05 d), else linearly interpolated.
    Compared: the share of each stream's total collected after t*.
    95% intervals two ways: counting noise only (model free), and counting noise
    plus the per-stream scatter floors fitted to this window (`floors`, Bq). The
    fits need those floors (counting-only chi2/dof is far above 1), so the
    floor-included interval is the one the verdict uses. Not applicable when
    k_top steps inside the window. Untestable when the interval is wider than 1
    or reaches outside [-1, 1], the range a share difference can take.
    """
    name, run, tmax, step, _ = w
    if step:
        return {"applicable": False, "reason": "k_top steps inside this window, so OV/IV need not be constant"}
    ds = ds_for(w, "fixedTBR", noise="counting", background=background)
    if "OV" not in ds.times or len(ds.times["OV"]) < 2:
        return {"applicable": False, "reason": "fewer than two OV samples"}
    n_iv = len(ds.times["IV"])
    tiv, tov = ds.times["IV"], ds.times["OV"]

    def iv_at(civ, t):
        k = np.where(np.abs(tiv - t) < 0.05 * 86400)[0]
        if len(k):
            return civ[k[0]]
        return float(np.interp(t, np.concatenate([[0.0], tiv]), np.concatenate([[0.0], civ])))

    civ0 = np.cumsum(ds.y[:n_iv])
    cov0 = np.cumsum(ds.y[n_iv:])
    js = [j for j in range(len(tov) - 1) if iv_at(civ0, tov[j]) >= 0.8 * civ0[-1]]
    if not js:
        return {"applicable": False, "reason": "no OV sample (other than the last) after IV reached 80%"}
    j = js[0]
    t_star = tov[j]

    def shares(y):
        civ = np.cumsum(y[:n_iv])
        cov_ = np.cumsum(y[n_iv:])
        return 1 - iv_at(civ, t_star) / civ[-1], 1 - cov_[j] / cov_[-1]

    iv0, ov0 = shares(ds.y)

    def interval(S):
        sims = np.random.default_rng(seed).multivariate_normal(ds.y, S, size=n_mc)
        diffs = np.array([b - a for a, b in (shares(y) for y in sims)])
        lo, hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
        # a share difference lies in [-1, 1], so an interval wider than 1 or reaching
        # outside that range is set by noise, not by the data
        noisy = hi - lo > 1.0 or lo < -1.0 or hi > 1.0
        v = "untestable" if noisy else ("consistent" if lo <= 0 <= hi else "OV late")
        return [lo, hi], v

    ci_c, v_c = interval(ds.C)
    S_f = ds.C + np.diag(np.array([floors[s] for s in ds.sid]) ** 2)
    ci_f, v_f = interval(S_f)
    return {"applicable": True, "t_star_day": t_star / 86400, "IV_share_after": iv0, "OV_share_after": ov0,
            "OV_minus_IV": ov0 - iv0, "floors_Bq": {k: float(floors[k]) for k in ("IV", "OV")},
            "ci95_counting": ci_c, "verdict_counting": v_c,
            "ci95_with_floors": ci_f, "verdict_with_floors": v_f,
            "verdict": v_f, "verdict_basis": "counting noise plus fitted scatter floors",
            "consistent_with_0D": v_f != "OV late"}


def pairwise_ratios(fits):
    """k_top ratio between every pair of primary windows, with 95% CI from the
    two profiles combined as independent (quadrature of half-widths in log10)."""
    out = {}
    names = list(fits)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pa, pb = fits[a]["profiles"]["log10_ktop"], fits[b]["profiles"]["log10_ktop"]
            if any(p["ci95"]["status"] != "identifiable" for p in (pa, pb)):
                continue
            ha = 0.5 * (pa["ci95"]["hi"] - pa["ci95"]["lo"])
            hb = 0.5 * (pb["ci95"]["hi"] - pb["ci95"]["lo"])
            d = pb["mle"] - pa["mle"]
            h = math.hypot(ha, hb)
            out[f"{b}/{a}"] = {"ratio": 10 ** d, "ci95": [10 ** (d - h), 10 ** (d + h)],
                               "consistent_with_1": (d - h) <= 0 <= (d + h)}
    return out


def run_transfer(res, pool, seeds, args, t_start):
    res["boot_reps"] = args.boot_reps
    res["power_reps"] = args.power_reps
    he = [WINDOWS[0], WINDOWS[1], WINDOWS[2]]
    tests = [("transferability_He_ktop", he, ["log10_ktop"]),
             ("transferability_run1_run3He_ktop", [WINDOWS[0], WINDOWS[2]], ["log10_ktop"]),
             ("transferability_run1_run2_ktop", he[:2], ["log10_ktop"]),
             ("transferability_He_kwall", he, ["log10_kwall"]),
             ("transferability_H2_run4_run6_ktop", [WINDOWS[3], WINDOWS[4]], ["log10_ktop"])]
    for key, grp, shared in tests:
        res[key] = transferability(pool, grp, shared, seeds(args.boot_reps), seeds(args.power_reps))
        t = res[key]
        print(f"{key}: LR={t['LR']:.2f} p_boot={t['p_bootstrap']:.3f} crit={t['bootstrap_crit95']:.2f} "
              f"power={t['power_if_second_run_differs_by']}", flush=True)
    print(f"[{time.time() - t_start:.0f}s] transferability done", flush=True)


def run_controls(res, pool, seeds, args, t_start):
    res["reps"] = args.reps
    res["flat_reps"] = args.flat_reps
    thetas = {}
    for w in WINDOWS:
        for var in ("fixedTBR", "freeTBR"):
            ds = ds_for(w, var)
            d = res["fits"][w[0]][var]["theta"]
            thetas[(w[0], var)] = np.array([d[n] for n in ds.names_A])
    ctrl = {}
    for w in WINDOWS:
        for var in ("fixedTBR", "freeTBR"):
            ctrl[f"{w[0]}/{var}"] = coverage(pool, w[0], var, thetas[(w[0], var)], seeds(args.reps))
            print(f"coverage {w[0]}/{var}: " + ", ".join(
                f"{n}={c['coverage95']}" for n, c in ctrl[f'{w[0]}/{var}'].items()
                if isinstance(c, dict) and "coverage95" in c), flush=True)
    res["synthetic_coverage"] = ctrl
    print(f"[{time.time() - t_start:.0f}s] coverage done", flush=True)
    # The OV floor collapses to its lower bound in some real fits. Test that
    # collapse: synthetic truth with a non-zero OV floor, the median of the
    # floors that did not collapse (TBR free fits).
    ov_name = "log10_sig_OV"
    live = [res["fits"][w[0]]["freeTBR"]["theta"][ov_name] for w in WINDOWS
            if "OV" not in res["fits"][w[0]]["freeTBR"]["floors_at_lower_bound"]]
    ov_true = float(np.median(live))
    ctrl2 = {"truth_log10_sig_OV": ov_true, "truth_sig_OV_Bq": 10 ** ov_true,
             "from_windows": [w[0] for w in WINDOWS
                              if "OV" not in res["fits"][w[0]]["freeTBR"]["floors_at_lower_bound"]],
             "cases": {}}
    for w in WINDOWS:
        for var in ("fixedTBR", "freeTBR"):
            if "OV" not in res["fits"][w[0]][var]["floors_at_lower_bound"]:
                continue
            ds = ds_for(w, var)
            th = thetas[(w[0], var)].copy()
            th[ds.names_A.index(ov_name)] = ov_true
            ctrl2["cases"][f"{w[0]}/{var}"] = coverage(pool, w[0], var, th, seeds(args.reps))
            print(f"coverage, non-zero OV floor, {w[0]}/{var}: " + ", ".join(
                f"{n}={c['coverage95']}" for n, c in ctrl2["cases"][f"{w[0]}/{var}"].items()
                if isinstance(c, dict) and "coverage95" in c), flush=True)
    res["synthetic_coverage_nonzero_OV_floor"] = ctrl2
    res["synthetic_flat_direction"] = flat_direction_control(
        pool, thetas[("run1", "freeTBR")], seeds(args.flat_reps))
    print("flat direction:", res["synthetic_flat_direction"]["IV_only_not_identifiable_fraction"],
          res["synthetic_flat_direction"]["IV_plus_OV_identifiable_fraction"], flush=True)


LAB_CHECK_VARIANTS = ("freeTBR", "labTBR")
LAB_CHECK_VARIANTS_QUENCH = ("freeTBR_quench", "labTBR_quench")


def _lab_point_one(args):
    """2*delta NLL at the lab's hand-set values: k_top alone (k_wall and the rest
    re-fitted), and k_top with k_wall together. A lab k_wall of 0 sits on the
    lower bound of log10 k_wall (1e-13 m/s, zero for these data)."""
    wname, var, theta, nll_min = args
    w = WIN[wname]
    background = "quench" if var.endswith("_quench") else "published"
    ds = ds_for(w, var.replace("_quench", ""), background=background)
    th = np.array(theta)
    m = ds.meta
    lkt = math.log10(m["lab_k_top_m_s"])
    lkw = math.log10(m["lab_k_wall_m_s"]) if m["lab_k_wall_m_s"] > 0 else F.BOUNDS["log10_kwall"][0]
    starts = (th, *ds.starts())
    _, v1 = F.constrained_fit(ds, "A", "log10_ktop", lkt, th, extra_starts=starts)
    _, v2 = F.constrained_fit(ds, "A", ["log10_ktop", "log10_kwall"], [lkt, lkw], th, extra_starts=starts)
    return wname, var, {"lab_k_top_m_s": m["lab_k_top_m_s"], "lab_k_wall_m_s": m["lab_k_wall_m_s"],
                        "two_delta_nll_at_lab_ktop": 2.0 * (v1 - nll_min),
                        "two_delta_nll_at_lab_ktop_and_kwall": 2.0 * (v2 - nll_min),
                        "lab_ktop_outside_95": 2.0 * (v1 - nll_min) > F.CHI2_95}


def lab_point_check(pool, fits):
    tasks = [(w[0], var, list(fits[w[0]][var]["theta"].values()), fits[w[0]][var]["nll"])
             for w in WINDOWS for var in LAB_CHECK_VARIANTS]
    tasks += [(w[0], var, list(fits[w[0]][var]["theta"].values()), fits[w[0]][var]["nll"])
              for w in QUENCH_WINDOWS for var in LAB_CHECK_VARIANTS_QUENCH]
    out = {}
    for wname, var, r in pool.map(_lab_point_one, tasks):
        out.setdefault(wname, {})[var] = r
    return out


def lab_style_cumfit(w):
    """Noise-model check: the lab tunes k_top and k_wall by eye against the
    cumulative IV and OV curves at TBR_used_in_model. Mimic that with an
    unweighted least-squares fit to the cumulative curves (no noise model)."""
    from scipy.optimize import minimize
    ds = ds_for(w, "labTBR")
    niv = len(ds.times["IV"])
    ycum = np.concatenate([np.cumsum(ds.y[:niv]), np.cumsum(ds.y[niv:])])
    fixed = {n: -1.0 for n in ds.names_A if n.startswith("log10_sig_")}

    def pred(z):
        d = {"log10_ktop": z[0], "log10_kwall": z[1], **fixed}
        if ds.step:
            d["log10_ktop_p2"] = z[0]
        p = ds.predict(np.array([d[n] for n in ds.names_A]))
        return np.concatenate([np.cumsum(p[:niv]), np.cumsum(p[niv:])])

    best = None
    for kt in (-8.0, -7.0, -6.0, -5.0):
        for kw in (-11.0, -9.0):
            r = minimize(lambda z: float(np.sum((pred(z) - ycum) ** 2)), [kt, kw], method="Powell",
                         bounds=[F.BOUNDS["log10_ktop"], F.BOUNDS["log10_kwall"]])
            if best is None or r.fun < best.fun:
                best = r
    return {"k_top_m_s": 10 ** best.x[0], "k_wall_m_s": 10 ** best.x[1], "sse": float(best.fun)}


def quench_summary():
    """The quench-matched totals from quench.json, pinned by its hash."""
    q = F.QUENCH
    if not q.get("runs"):
        sys.exit("quench.json missing: run quench.py first")
    keep = ("quench_set", "measured_TBR_published", "high_tSIE_vials", "high_tSIE_vials_by_position",
            "vials_beyond_scan_tSIE", "blank_minus_curve_at_blank_tSIE_Bq", "published", "anchored",
            "anchored_no_clamp", "literal", "anchored_flat_extrapolation_IV_plus_OV_Bq",
            "validation_max_abs_diff_Bq", "reason")
    return {"quench_json_sha256": file_sha256(HERE / "quench.json"), "scan": q["scan"], "n_mc": q["n_mc"],
            "high_tSIE": q["high_tSIE"],
            "runs": {r: {k: v for k, v in d.items() if k in keep} for r, d in q["runs"].items()}}

def salt_inventory():
    """Logged salt mass against the 1 L column the 0D and 1D models use (case.yaml)."""
    import yaml
    c = yaml.safe_load(CASE_YAML.read_text())
    m = c["material"]["salt_mass_charged_kg"]["value"]
    rho = c["material"]["density_g_cm3"]["value"]
    h_best = c["variants"]["best_estimate"]["salt_height_cm"]["value"] / 100.0
    v_logged = m / (rho * 1000.0)  # m3
    return {"salt_mass_logged_kg": m, "density_g_cm3": rho,
            "density_temperature_C": c["material"]["density_temperature_c"]["value"],
            "volume_logged_m3": v_logged, "volume_model_m3": M.V_M3,
            "height_model_m": M.H_M, "height_best_estimate_m": h_best,
            "height_best_estimate_locator": c["variants"]["best_estimate"]["salt_height_cm"]["locator"],
            "volume_ratio_logged_over_model": v_logged / M.V_M3,
            "height_ratio_best_over_model": h_best / M.H_M,
            "diffusivity_scale_height_squared": (h_best / M.H_M) ** 2}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10, help="synthetic replicates per coverage control")
    ap.add_argument("--flat-reps", type=int, default=None,
                    help="replicates for the flat-direction control (default: max(reps // 5, 8))")
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--boot-reps", type=int, default=20, help="bootstrap replicates per LR test")
    ap.add_argument("--power-reps", type=int, default=None,
                    help="replicates per power factor (default: max(boot_reps // 2, 10))")
    ap.add_argument("--part", choices=["fits", "controls", "transfer"], required=True,
                    help="fits: fits, profiles, Fisher, OV test. controls: synthetic coverage and "
                         "flat-direction controls (needs fits). transfer: likelihood-ratio tests with "
                         "bootstrap. Three parts keep each invocation under ~5 min on an 8-core laptop.")
    args = ap.parse_args()
    if args.flat_reps is None:
        args.flat_reps = max(args.reps // 5, 8)
    if args.power_reps is None:
        args.power_reps = max(args.boot_reps // 2, 10)
    t_start = time.time()

    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    part_args = {"fits": {"workers": args.workers, "seed": args.seed},
                 "controls": {"reps": args.reps, "flat_reps": args.flat_reps, "workers": args.workers,
                              "seed": args.seed},
                 "transfer": {"boot_reps": args.boot_reps, "power_reps": args.power_reps,
                              "workers": args.workers, "seed": args.seed}}[args.part]
    prov = provenance(part_args)
    if args.part in ("controls", "transfer"):
        if "fits" not in res:
            sys.exit("run --part fits first")
        if res.get("provenance", {}).get("fits", {}).get("increments_sha256") != prov["increments_sha256"]:
            sys.exit("increments.json changed since --part fits; rerun --part fits first")
        if res.get("provenance", {}).get("fits", {}).get("quench_sha256") != prov["quench_sha256"]:
            sys.exit("quench.json changed since --part fits; rerun --part fits first")
        prov["fits_sha256"] = fits_hash(res["fits"])
    else:
        # new fits invalidate controls and transfer computed from the old ones
        stale = [k for k in res if k.startswith(("synthetic_", "transferability_"))]
        for k in stale + ["reps", "boot_reps", "flat_reps", "power_reps"]:
            res.pop(k, None)
        for k in ("controls", "transfer"):
            res.get("provenance", {}).pop(k, None)
    res.setdefault("provenance", {})[args.part] = prov
    res.update({"_generated": "Written by analysis.py (three parts). Do not edit by hand.",
           "seed": args.seed,
           "geometry": {"A_top_m2": M.A_TOP, "A_wall_m2": M.A_WALL, "V_m3": M.V_M3, "height_m": M.H_M},
           "solver_vs_release_r0_max_rel_gap": M.check_against_r0(),
           "poisson_recount_check_run1": F.INC["poisson_recount_check_run1"],
           "windows": {w[0]: {"run": w[1], "t_max_day": w[2], "gas": w[4],
                              "lab_k_top_m_s": F.INC["runs"][str(w[1])]["lab_k_top_m_s"],
                              "lab_k_wall_m_s": F.INC["runs"][str(w[1])]["lab_k_wall_m_s"],
                              "measured_TBR": F.INC["runs"][str(w[1])]["measured_TBR"],
                              "neutron_rate_rel_err": F.INC["runs"][str(w[1])]["neutron_rate_err_n_s"]
                              / F.INC["runs"][str(w[1])]["neutron_rate_n_s"]}
                         for w in WINDOWS + SECONDARY}})

    pool = make_pool(args.workers)
    ss = np.random.SeedSequence(args.seed)

    def seeds(n):
        return [int(x) for x in ss.spawn(1)[0].generate_state(n)]

    if args.part in ("transfer", "controls"):
        (run_transfer if args.part == "transfer" else run_controls)(res, pool, seeds, args, t_start)
        pool.shutdown()
        res[f"runtime_{args.part}_s"] = time.time() - t_start
        OUT.write_text(json.dumps(res, indent=1, default=float) + "\n")
        print(f"wrote {OUT} in {res[f'runtime_{args.part}_s']:.0f}s")
        return 0
    fits, thetas = {w[0]: {} for w in WINDOWS}, {}
    tasks = ([(w[0], var) for w in WINDOWS for var in VARIANTS] + [(w[0], "freeTBR") for w in SECONDARY]
             + [(w[0], v + "_quench") for w in QUENCH_WINDOWS for v in ("freeTBR", "labTBR", "fixedTBR")]
             + [(w[0], "freeTBR_wideTBR") for w in WINDOWS] + [("run1", "IVonly_freeTBR_wideTBR")])
    for wname, var, th, blk in pool.map(_fit_task, tasks):
        fits.setdefault(wname, {})[var] = blk
        thetas[(wname, var)] = np.array(th)
        print(f"{wname:14s} {var:15s} n={blk['n_points']:2d} k_top={blk['physical']['k_top_m_s']} "
              f"({blk['seconds']:.1f}s)", flush=True)
    for w in WINDOWS:
        fits[w[0]]["counting_only_gof"] = gof_counting(w)
    res["fits"] = fits
    res["ov_lag_test"] = {w[0]: ov_lag_test(w, fits[w[0]]["freeTBR"]["physical"]["noise_floor_Bq"])
                          for w in WINDOWS}
    res["ov_lag_test_quench"] = {
        w[0]: ov_lag_test(w, fits[w[0]]["freeTBR_quench"]["physical"]["noise_floor_Bq"], background="quench")
        for w in QUENCH_WINDOWS}
    res["pairwise_ktop_ratios_freeTBR"] = pairwise_ratios({w[0]: fits[w[0]]["freeTBR"] for w in WINDOWS})
    res["pairwise_ktop_ratios_fixedTBR"] = pairwise_ratios({w[0]: fits[w[0]]["fixedTBR"] for w in WINDOWS})
    print(f"[{time.time() - t_start:.0f}s] fits done", flush=True)
    res["lab_point_check"] = lab_point_check(pool, fits)
    res["lab_style_cumfit"] = {w[0]: lab_style_cumfit(w) for w in WINDOWS}
    res["quench"] = quench_summary()
    res["pairwise_ktop_ratios_labTBR"] = pairwise_ratios({w[0]: fits[w[0]]["labTBR"] for w in WINDOWS})
    print(f"[{time.time() - t_start:.0f}s] lab-point checks done", flush=True)

    pool.shutdown()
    res["salt_inventory"] = salt_inventory()
    res["runtime_fits_s"] = time.time() - t_start
    prov["fits_sha256"] = fits_hash(res["fits"])
    OUT.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(f"wrote {OUT} in {res['runtime_fits_s']:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
