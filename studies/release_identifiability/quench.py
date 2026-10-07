#!/usr/bin/env python3
"""Apply the lab's run-6 quench-matched background to runs 3 and 4.

Up to run 4 the lab subtracts one blank vial (tSIE ~290) from every sample vial
in the same LSC file. A vial's background in Bq depends on its tSIE: the run-6
blank scan (SLR_BLANK_SCAN.csv, 2 to 11 mL water in 10 mL cocktail) goes from
0.31 Bq at tSIE 278 to 0.17 Bq at tSIE 489. From run 6 on, the lab subtracts
background_curve(tSIE), a linear interpolation of that scan (BABY-1L-run-6
analysis/tritium/tritium_model.py, build_background_curve_from_file). The
high-tSIE vials are the first bubbler vial of a sample (run 2's general.json:
"Bubbler water in vials 2-X-1 replenished periodically to compensate for
evaporation"), so a fixed blank over-subtracts exactly those vials.

This script executes the lab's scripts unchanged (through build_increments.py's
wrappers), takes background_curve from the lab's own run-6 namespace, and
recomputes every vial of runs 3 and 4 two ways:

  anchored  blank_Bq * curve(tSIE_vial) / curve(tSIE_blank): the curve's shape,
            the level of the blank counted with that vial (main result)
  literal   curve(tSIE_vial): the run-6 procedure as written

Both keep the lab's clamp of negative net vials to zero. Checks:
  - the rebuilt published numbers equal increments.json (1e-9 Bq)
  - on run 6, both variants reproduce processed_data.json exactly

Runs 1 and 2 are not corrected. They were counted with a different quench set
(3H-UG-LLT-NM, not -NM-C as in the scan and runs 3 to 6), whose efficiency at
high tSIE differs (about 0.47 against 0.34 at tSIE ~500), CPMA is exported as
whole counts, and their blanks sit 0.06 to 0.27 Bq above the scan curve. A
faithful correction needs a blank scan with that quench set.

Uncertainties: Monte Carlo over the counting noise of every sample vial, every
in-file blank (shared by all vials that used it) and every scan point.

Writes quench.json. Run with the study venv:
  .venv/bin/python quench.py
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_increments as bi  # noqa: E402  (patches libra_toolbox on import)
from fetch_lsc import folder  # noqa: E402

lsc = bi.lsc
OUT = HERE / "quench.json"
OUT_VARIANTS = HERE / "quench_variants.json"
CORRECTED = (3, 4)
VALIDATE = 6
ALL = (1, 2, 3, 4, 6)
N_MC = 2000
SEED = 20261007
HIGH_TSIE = 330.0  # a vial is "high tSIE" above this (blanks sit at ~280 to 295)

_prev_sub = lsc.LSCSample.substract_background


def _sub(self, b):
    self._bg_sample = b
    _prev_sub(self, b)


lsc.LSCSample.substract_background = _sub


def _row(reader, label):
    return reader.data.iloc[bi._label_map(reader)[label]]


def _sigma(bq, cpm, minutes):
    """Poisson sigma of an LSC activity, as in build_increments._row_sigmas."""
    return abs(bq) / math.sqrt(max(cpm, 1e-9) * minutes)


def collect(run):
    """Per vial: gross Bq, the background the lab subtracted, both tSIE values."""
    ns = bi.exec_lab_script(run)
    streams = ns.get("gas_streams") or {"IV": ns["IV_stream"], "OV": ns["OV_stream"]}
    vials = []
    for s in ("IV", "OV"):
        for j, smp in enumerate(streams[s].samples):
            for k, v in enumerate(smp.samples):
                r = _row(v._reader, v._label)
                own = [key for key in v.terms if key[3] == v._label]
                d = {"stream": s, "sample": j, "vial": k + 1, "gross": float(v.gross),
                     "net_lab": float(v.activity.magnitude), "tsie": float(r["tSIE"]),
                     "sig_gross": math.sqrt(bi.ROW_VAR[own[0]]),
                     "quench": str(v._reader.quench_set)}
                b = getattr(v, "_bg_sample", None)
                if b is not None and hasattr(b, "_reader"):
                    rb = _row(b._reader, b._label)
                    d.update(bg=float(b.activity.magnitude), bg_tsie=float(rb["tSIE"]),
                             bg_key=f"{b._reader.file_path}|{b._label}",
                             sig_bg=math.sqrt(sum(c * c * bi.ROW_VAR[key] for key, c in b.terms.items())))
                elif b is not None:  # a typed scalar background (run 1 samples 1-2)
                    d.update(bg=float(b.activity.magnitude), bg_tsie=None, bg_key="scalar", sig_bg=0.0)
                else:  # run 6: the lab subtracted background_curve(tSIE) outside the toolbox
                    d.update(bg=float(ns["background_curve"](d["tsie"])), bg_tsie=None, bg_key="curve", sig_bg=0.0)
                vials.append(d)
    return ns, vials


def scan_points(ns6):
    """The run-6 blank scan exactly as the lab reads it, plus counting sigmas."""
    info = ns6["general_data"]["tritium_blank_set"]
    reader = ns6["background_reader"]
    pts = []
    for label in info["blanks"]:
        r = ns6["get_row_by_label"](reader, label)
        bq = float(lsc.LSCSample.from_file(reader, label).activity.magnitude)
        pts.append((float(r["tSIE"]), bq, _sigma(bq, float(r["CPMA"]), float(r["Count Time"]))))
    return np.array(sorted(pts))


def curve_from(points_t, points_bq):
    from scipy.interpolate import interp1d
    return interp1d(points_t, points_bq, bounds_error=False, fill_value="extrapolate")


def backgrounds(vials, curve, mode, bg=None):
    bg = np.array([v["bg"] for v in vials]) if bg is None else bg
    tv = np.array([v["tsie"] for v in vials])
    if mode == "published":
        return bg
    if mode == "literal":
        return curve(tv)
    tb = np.array([v["bg_tsie"] if v["bg_tsie"] is not None else np.nan for v in vials])
    anchored = bg * curve(tv) / curve(tb)
    if mode == "anchored_high_tSIE":  # only vials above HIGH_TSIE change
        return np.where(tv > HIGH_TSIE, anchored, bg)
    return anchored


def per_sample(vials, net):
    """Sum vials into (stream -> [sample totals]), plus soluble (vials 1-2) sums."""
    tot, sol = {"IV": {}, "OV": {}}, {"IV": {}, "OV": {}}
    for v, x in zip(vials, net):
        tot[v["stream"]][v["sample"]] = tot[v["stream"]].get(v["sample"], 0.0) + x
        if v["vial"] <= 2:
            sol[v["stream"]][v["sample"]] = sol[v["stream"]].get(v["sample"], 0.0) + x
    as_list = lambda d: [d[k] for k in sorted(d)]  # noqa: E731
    return ({s: as_list(tot[s]) for s in tot}, {s: as_list(sol[s]) for s in sol})


def summary(vials, net):
    tot, sol = per_sample(vials, net)
    out = {}
    for s in ("IV", "OV"):
        t, so = sum(tot[s]), sum(sol[s])
        out[s] = {"total_Bq": t, "soluble_share": so / t if t > 0 else float("nan")}
    out["IV_plus_OV_Bq"] = out["IV"]["total_Bq"] + out["OV"]["total_Bq"]
    return out


def evaluate(vials, curve, mode, clamp=True, gross=None, bg=None):
    g = np.array([v["gross"] for v in vials]) if gross is None else gross
    net = g - backgrounds(vials, curve, mode, bg)
    return np.maximum(net, 0.0) if clamp else net


def diagnostics(vials, net):
    """Low-activity vials (gross < 0.6 Bq) should not trend with tSIE if the
    background is right. Their mean net at high tSIE cannot be negative."""
    g = np.array([v["gross"] for v in vials])
    tv = np.array([v["tsie"] for v in vials])
    low, hi = g < 0.6, tv > HIGH_TSIE
    A = np.vstack([np.ones(low.sum()), tv[low] - 300.0]).T
    slope = float(np.linalg.lstsq(A, net[low], rcond=None)[0][1]) * 100.0
    return {"negative_vials": int(np.sum(net < 0)), "vials": len(net),
            "low_activity_vials": int(low.sum()),
            "slope_Bq_per_100_tSIE": slope,
            "mean_net_low_activity_high_tSIE_Bq": float(np.mean(net[low & hi])) if (low & hi).any() else None,
            "mean_net_low_activity_low_tSIE_Bq": float(np.mean(net[low & ~hi]))}


# Shapes of the background curve below tSIE ~300, where the scan has two points
# (0.306 Bq at 278, 0.266 at 294). Each builds a callable from scan points, sigmas.
FLAT_BELOW = 300.0


def _interp(t, b, s):
    return curve_from(t, b)


def _flat_below(t, b, s):
    c = curve_from(t, b)
    return lambda x: c(np.maximum(np.asarray(x, float), FLAT_BELOW))


def _poly(deg):
    def build(t, b, s):
        p = np.polyfit(t, b, deg, w=1.0 / s)
        return lambda x: np.polyval(p, np.asarray(x, float))
    return build


SHAPES = {"flat_below_300": _flat_below, "linear": _poly(1), "quadratic": _poly(2), "interp": _interp}


def monte_carlo(vials, pts, mode, rng, n=N_MC, shape=None):
    """Spread of the stream totals from counting noise (vials, blanks, scan)."""
    g0 = np.array([v["gross"] for v in vials])
    sg = np.array([v["sig_gross"] for v in vials])
    keys = sorted({v["bg_key"] for v in vials})
    kidx = np.array([keys.index(v["bg_key"]) for v in vials])
    bg0 = np.array([v["bg"] for v in vials])
    sb = np.array([v["sig_bg"] for v in vials])
    draws = []
    for _ in range(n):
        g = g0 + rng.standard_normal(len(g0)) * sg
        z = rng.standard_normal(len(keys))[kidx]
        bg = bg0 + z * sb
        bq = pts[:, 1] + rng.standard_normal(len(pts)) * pts[:, 2]
        c = curve_from(pts[:, 0], bq) if shape is None else shape(pts[:, 0], bq, pts[:, 2])
        s = summary(vials, evaluate(vials, c, mode, True, g, bg))
        draws.append([s["IV"]["total_Bq"], s["OV"]["total_Bq"], s["IV_plus_OV_Bq"],
                      s["IV"]["soluble_share"], s["OV"]["soluble_share"]])
    sd = np.std(np.array(draws), axis=0, ddof=1)
    return {"IV_total_Bq": sd[0], "OV_total_Bq": sd[1], "IV_plus_OV_Bq": sd[2],
            "IV_soluble_share": sd[3], "OV_soluble_share": sd[4]}


def main() -> int:
    rng = np.random.default_rng(SEED)
    inc = json.loads((HERE / "increments.json").read_text())
    data = {r: collect(r) for r in ALL}
    ns6 = data[VALIDATE][0]
    pts = scan_points(ns6)
    curve = curve_from(pts[:, 0], pts[:, 1])
    lab_curve = ns6["background_curve"]
    tt = np.linspace(pts[0, 0], 700.0, 50)
    if np.max(np.abs(curve(tt) - lab_curve(tt))) > 1e-12:
        raise RuntimeError("scan curve differs from the lab's background_curve")
    scan_file = folder(VALIDATE) / "data" / "tritium_detection" / ns6["general_data"]["tritium_blank_set"]["filename"]

    out = {"_generated": "Written by quench.py. Do not edit by hand.",
           "scan": {"file": f"BABY-1L-run-6/data/tritium_detection/{scan_file.name}",
                    "sha256": hashlib.sha256(scan_file.read_bytes()).hexdigest(),
                    "points_tSIE_Bq_sigma": pts.tolist(),
                    "quench_set": next(v["quench"] for v in data[VALIDATE][1])},
           "increments_sha256": hashlib.sha256((HERE / "increments.json").read_bytes()).hexdigest(),
           "n_mc": N_MC, "seed": SEED, "high_tSIE": HIGH_TSIE, "runs": {}}

    for run in ALL:
        ns, vials = data[run]
        proc = json.loads((folder(run) / "data" / "processed_data.json").read_text())
        pub = evaluate(vials, curve, "published")
        if np.max(np.abs(pub - np.array([v["net_lab"] for v in vials]))) > 1e-9:
            raise RuntimeError(f"run {run}: published nets not reproduced")
        tot_pub, _ = per_sample(vials, pub)
        for s in ("IV", "OV"):
            d = np.max(np.abs(np.array(tot_pub[s]) - np.array(inc["runs"][str(run)]["streams"][s]["inc_total_Bq"])))
            if d > 1e-9:
                raise RuntimeError(f"run {run} {s}: per-sample totals differ from increments.json by {d}")
        tv = np.array([v["tsie"] for v in vials])
        hi = tv > HIGH_TSIE
        r = {"quench_set": sorted({v["quench"] for v in vials}),
             "measured_TBR_published": proc["measured_TBR"]["value"],
             "high_tSIE_vials": int(hi.sum()),
             "high_tSIE_vials_by_position": {f"{s}{k}": int(np.sum(hi & (np.array([v['stream'] for v in vials]) == s)
                                                                 & (np.array([v['vial'] for v in vials]) == k)))
                                             for s in ("IV", "OV") for k in (1, 2, 3, 4)},
             "vials_beyond_scan_tSIE": int(np.sum(tv > pts[-1, 0])),
             "blank_minus_curve_at_blank_tSIE_Bq": None,
             "published": {**summary(vials, pub), "sd": monte_carlo(vials, pts, "published", rng),
                           "diagnostics": diagnostics(vials, evaluate(vials, curve, "published", clamp=False))}}
        tb = [v for v in vials if v["bg_tsie"] is not None]
        if tb:
            dd = [v["bg"] - float(curve(v["bg_tsie"])) for v in tb]
            r["blank_minus_curve_at_blank_tSIE_Bq"] = {"median": float(np.median(dd)), "min": float(min(dd)),
                                                      "max": float(max(dd))}
        if run == VALIDATE:
            # the lab's own procedure: the published nets ARE curve(tSIE) subtracted and clamped,
            # so this only checks the plumbing (close to tautological). The real validations
            # are the published-net check above and build_increments' match to processed_data.json.
            lit = evaluate(vials, curve, "literal")
            r["validation_max_abs_diff_Bq"] = float(np.max(np.abs(lit - pub)))
            if r["validation_max_abs_diff_Bq"] > 1e-12:
                raise RuntimeError("run 6: the run-6 procedure does not reproduce the lab's nets")
            r["anchored_delta_Bq"] = {s: [0.0] * len(tot_pub[s]) for s in ("IV", "OV")}
        elif run in CORRECTED:
            for mode, clamp, key in (("anchored", True, "anchored"), ("anchored", False, "anchored_no_clamp"),
                                     ("literal", True, "literal")):
                net = evaluate(vials, curve, mode, clamp)
                r[key] = {**summary(vials, net),
                          "diagnostics": diagnostics(vials, evaluate(vials, curve, mode, clamp=False))}
                if clamp:
                    r[key]["sd"] = monte_carlo(vials, pts, mode, rng)
                r[key]["measured_TBR"] = r["measured_TBR_published"] * r[key]["IV_plus_OV_Bq"] / r["published"]["IV_plus_OV_Bq"]
            # flat extrapolation beyond the last scan point instead of linear
            flat = curve_from(np.append(pts[:, 0], 5000.0), np.append(pts[:, 1], pts[-1, 1]))
            r["anchored_flat_extrapolation_IV_plus_OV_Bq"] = summary(vials, evaluate(vials, flat, "anchored"))["IV_plus_OV_Bq"]
            tot_anc, _ = per_sample(vials, evaluate(vials, curve, "anchored"))
            r["anchored_delta_Bq"] = {s: (np.array(tot_anc[s]) - np.array(tot_pub[s])).tolist() for s in ("IV", "OV")}
        else:
            r["anchored_delta_Bq"] = None
            r["reason"] = ("counted with quench set 3H-UG-LLT-NM, not the -NM-C set of the run-6 blank scan; "
                           "CPMA exported as whole counts; in-file blanks sit above the scan curve. "
                           "A correction needs a blank scan with this quench set.")
        out["runs"][str(run)] = r
        print(f"run {run}: done", flush=True)

    OUT.write_text(json.dumps(out, indent=1, default=float) + "\n")
    write_variants(data, pts, out)
    for run in (3, 4):
        r = out["runs"][str(run)]
        p, a = r["published"], r["anchored"]
        print(f"run {run}: IV {p['IV']['total_Bq']:.2f} -> {a['IV']['total_Bq']:.2f}+-{a['sd']['IV_total_Bq']:.2f}  "
              f"OV {p['OV']['total_Bq']:.2f} -> {a['OV']['total_Bq']:.2f}+-{a['sd']['OV_total_Bq']:.2f}  "
              f"TBR {r['measured_TBR_published']:.3e} -> {a['measured_TBR']:.3e}")
    print(f"wrote {OUT}")
    return 0


def write_variants(data, pts, out):
    """Which vials carry the change, and how it depends on the curve's shape below
    tSIE 300, where the scan has two points. Two vial sets (only the vials above
    HIGH_TSIE, or all vials) times four shapes. A separate file and RNG, so
    quench.json (and the fits built on it) do not move."""
    rng = np.random.default_rng(SEED + 1)
    res = {"_generated": "Written by quench.py. Do not edit by hand.", "high_tSIE": HIGH_TSIE,
           "flat_below": FLAT_BELOW,
           "shapes": {"flat_below_300": f"lab's curve held flat below tSIE {FLAT_BELOW:.0f}",
                      "linear": "weighted straight-line fit to the scan",
                      "quadratic": "weighted quadratic fit to the scan",
                      "interp": "lab's curve (linear interpolation of the scan)"},
           "vial_sets": {"anchored_high_tSIE": f"only vials above tSIE {HIGH_TSIE:.0f}", "anchored": "all vials"},
           "n_mc": N_MC, "runs": {}}
    for run in CORRECTED:
        _, vials = data[run]
        pub_tot = out["runs"][str(run)]["published"]["IV_plus_OV_Bq"]
        tbr0 = out["runs"][str(run)]["measured_TBR_published"]
        tv = np.array([v["tsie"] for v in vials])
        hi = tv > HIGH_TSIE
        pos = np.array([f"{v['stream']}{v['vial']}" for v in vials])
        r = {"high_tSIE_vials": int(hi.sum()),
             "high_tSIE_first_bubbler_vials": int(np.sum(hi & np.isin(pos, ["IV1", "OV1"]))),
             "other_vial_tSIE_range": [float(tv[~hi].min()), float(np.median(tv[~hi])), float(tv[~hi].max())],
             "blank_tSIE_range": [float(min(v["bg_tsie"] for v in vials)), float(max(v["bg_tsie"] for v in vials))],
             "variants": {}}
        pub = evaluate(vials, None, "published")
        for mode in ("anchored_high_tSIE", "anchored"):
            for name, shape in SHAPES.items():
                c = shape(pts[:, 0], pts[:, 1], pts[:, 2])
                net = evaluate(vials, c, mode)
                sm = summary(vials, net)
                r["variants"][f"{mode}/{name}"] = {
                    **sm, "measured_TBR": tbr0 * sm["IV_plus_OV_Bq"] / pub_tot,
                    "delta_Bq": float((net - pub).sum()),
                    "diagnostics": diagnostics(vials, evaluate(vials, c, mode, clamp=False)),
                    "sd": monte_carlo(vials, pts, mode, rng, shape=shape)}
        res["runs"][str(run)] = r
        print(f"run {run}: TBR {tbr0:.3e} -> " + ", ".join(f"{k} {x['measured_TBR'] / tbr0 - 1:+.1%}"
                                                          for k, x in r["variants"].items()))
    OUT_VARIANTS.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(f"wrote {OUT_VARIANTS}")


if __name__ == "__main__":
    sys.exit(main())
