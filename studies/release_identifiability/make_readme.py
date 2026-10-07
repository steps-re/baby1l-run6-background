#!/usr/bin/env python3
"""Write README.md (one screen) and APPENDIX.md (everything else) from the
study's JSON outputs. No result is typed by hand; typed qualitative sentences
are guarded and the build fails if the data stop supporting them.

Inputs: increments.json, quench.json (via results.json), quench_variants.json, change.json,
results.json, diffusivity.json, festim_check.json.
Usage: python3 make_readme.py [--draft]
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
INC = json.loads((HERE / "increments.json").read_text())
R = json.loads((HERE / "results.json").read_text())
DF = json.loads((HERE / "diffusivity.json").read_text())
FC = json.loads((HERE / "festim_check.json").read_text())
QV_ = json.loads((HERE / "quench_variants.json").read_text())
CH = json.loads((HERE / "change.json").read_text())

W = ["run1", "run2", "run3_He", "run4_1000ppm", "run6"]
NAME = {"run1": "Run 1", "run2": "Run 2", "run3_He": "Run 3 (He part)", "run4_1000ppm": "Run 4 (1000 ppm part)",
        "run6": "Run 6", "run3_full_step": "Run 3 (whole run)", "run4_full_step": "Run 4 (whole run)"}
GAS = {"run1": "He", "run2": "He", "run3_He": "He", "run4_1000ppm": "1000 ppm H2 in He", "run6": "3.5% H2, then He"}


CHI2_95 = 3.841458820694124


def fail(msg):
    raise SystemExit(f"refusing to build README.md: {msg}")


def fits_hash(fits):
    """Same hash as analysis.fits_hash (sorted-key JSON of the stored fits)."""
    return hashlib.sha256(json.dumps(fits, sort_keys=True).encode()).hexdigest()


def check_provenance():
    """Every part must exist, come from the current increments.json, and controls
    and transfer must have been computed from the fits stored now."""
    inc_sha = hashlib.sha256((HERE / "increments.json").read_bytes()).hexdigest()
    pr, pd = R.get("provenance", {}), DF.get("provenance", {})
    for src, part, prov in (("results.json", k, pr.get(k)) for k in ("fits", "controls", "transfer")):
        if not prov:
            fail(f"{src} has no provenance for part '{part}'; rerun it")
        if prov["increments_sha256"] != inc_sha:
            fail(f"{src} part '{part}' was computed from a different increments.json")
    for part in ("windows", "controls"):
        if not pd.get(part):
            fail(f"diffusivity.json has no provenance for part '{part}'; rerun it")
        if pd[part]["increments_sha256"] != inc_sha:
            fail(f"diffusivity.json part '{part}' was computed from a different increments.json")
    h = fits_hash(R["fits"])
    for part in ("fits", "controls", "transfer"):
        if pr[part].get("fits_sha256") != h:
            fail(f"results.json part '{part}' was computed from different fits than the ones stored; "
                 "rerun controls and transfer after fits")
    if R["reps"] != pr["controls"]["arguments"]["reps"] or R["boot_reps"] != pr["transfer"]["arguments"]["boot_reps"]:
        fail("replicate counts in results.json disagree with their provenance")
    for k, c in R["synthetic_coverage"].items():
        if c["log10_ktop"]["reps"] != R["reps"]:
            fail(f"coverage {k} has {c['log10_ktop']['reps']} replicates, provenance says {R['reps']}")
    for k in (k for k in R if k.startswith("transferability_")):
        t = R[k]
        if t["bootstrap_reps"] != R["boot_reps"]:
            fail(f"{k} has {t['bootstrap_reps']} bootstrap replicates, provenance says {R['boot_reps']}")
        for w, th in t["separate_fit_theta"].items():
            stored = list(R["fits"][w]["freeTBR"]["theta"].values())
            gap = max(abs(a_ - b_) for a_, b_ in zip(th, stored))
            if gap > 1e-3:
                fail(f"{k}: the separate fit of {w} differs from the stored fit by {gap:.3g} in log10")
    q_sha = hashlib.sha256((HERE / "quench.json").read_bytes()).hexdigest()
    if pr["fits"].get("quench_sha256") != q_sha or R["quench"]["quench_json_sha256"] != q_sha:
        fail("results.json was computed from a different quench.json; rerun analysis.py --part fits")
    chp = CH["provenance"]
    if chp["increments_sha256"] != inc_sha:
        fail("change.json was computed from a different increments.json; rerun change.py")
    run1_sha = hashlib.sha256(json.dumps(R["fits"]["run1"], sort_keys=True).encode()).hexdigest()
    if chp["results_fits_run1_sha256"] != run1_sha:
        fail("change.json was computed against different run 1 fits; rerun change.py")
    for lab, c in DF["synthetic_control_run1_design"].items():
        if isinstance(c, dict) and c["reps"] != pd["controls"]["arguments"]["control_reps"]:
            fail(f"diffusivity control {lab} replicate count disagrees with provenance")


def names_list(ws, empty="none"):
    """'Run 1', 'Run 1 and Run 2', 'Run 1, Run 2 and Run 4'."""
    n = [NAME[w] for w in ws]
    if not n:
        return empty
    return n[0] if len(n) == 1 else ", ".join(n[:-1]) + " and " + n[-1]


def e(x, d=2):
    """Compact scientific notation, e.g. 1.03e-7."""
    if x is None:
        return "open"
    m, ex = f"{x:.{d}e}".split("e")
    return f"{m}e{int(ex)}"


def pct(x, d=0):
    return f"{100 * x:.{d}f}%"


def ci_lin(block, name, d=2, plain=False):
    """'value (lo to hi)' in linear units from a log10 profile. Open sides are named."""
    p = block["profiles"][name]
    c = p["ci95"]
    lo = None if (c["lo"] is None or c["status"] in ("upper bound only", "not identifiable")) else 10 ** c["lo"]
    hi = None if (c["hi"] is None or c["status"] in ("lower bound only", "not identifiable")) else 10 ** c["hi"]
    v = 10 ** p["mle"]
    if plain:
        e_ = lambda x, _d: "open" if x is None else f"{x:.2f}"  # noqa: E731
    else:
        e_ = e
    if lo is None and hi is None:
        return f"{e_(v, d)} (no bound from data)"
    if lo is None:
        return f"{e_(v, d)} (below {e_(hi, d)})"
    if hi is None:
        return f"{e_(v, d)} (above {e_(lo, d)})"
    return f"{e_(v, d)} ({e_(lo, d)} to {e_(hi, d)})"


def status(block, name):
    return block["profiles"][name]["ci95"]["status"]


def tau_ci(block):
    p = block["profiles"]["log10_alpha"]
    c = p["ci95"]
    t = lambda la: 1 / 10 ** la / 86400  # noqa: E731
    if c["status"] != "identifiable":
        return f"{t(p['mle']):.2f} ({c['status']})"
    return f"{t(p['mle']):.2f} ({t(c['hi']):.2f} to {t(c['lo']):.2f})"


def fits(w, v="freeTBR"):
    return R["fits"][w][v]


MIN_REPS = {"controls reps": 1000, "transfer boot_reps": 500, "diffusivity control_reps": 100,
            "change null reps": 200}
DRAFT = False


def check_replicates(draft):
    """The claims (p-values, coverage, power) need these replicate counts. Below
    them the README only builds with --draft, which stamps it DRAFT."""
    have = {"controls reps": R["reps"], "transfer boot_reps": R["boot_reps"],
            "diffusivity control_reps": DF["provenance"]["controls"]["arguments"]["control_reps"],
            "change null reps": CH["null"]["reps"]}
    low = {k: (have[k], m) for k, m in MIN_REPS.items() if have[k] < m}
    if low and not draft:
        fail("replicate counts below what the claims need: "
             + ", ".join(f"{k} = {h} (need {m})" for k, (h, m) in low.items())
             + ". Rerun with `make reproduce-full` or pass --draft.")
    return bool(low)


def main():
    global DRAFT
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--draft", action="store_true", help="build below the minimum replicate counts, stamped DRAFT")
    args = ap.parse_args()
    check_provenance()
    DRAFT = check_replicates(args.draft)
    L = []
    a = L.append
    f1, f2, f3, f4 = (fits(w) for w in W[:4])
    tr = R["transferability_He_ktop"]
    r12 = R["pairwise_ktop_ratios_freeTBR"]["run2/run1"]
    r13 = R["pairwise_ktop_ratios_freeTBR"]["run3_He/run1"]
    r41 = R["pairwise_ktop_ratios_freeTBR"]["run4_1000ppm/run1"]
    r61 = R["pairwise_ktop_ratios_freeTBR"]["run6/run1"]
    t13 = R["transferability_run1_run3He_ktop"]
    tH2 = R["transferability_H2_run4_run6_ktop"]
    tkw = R["transferability_He_kwall"]
    gof = {w: R["fits"][w]["counting_only_gof"] for w in W}
    lag = R["ov_lag_test"]
    # verdicts: "verdict_with_floors" (the one the README concludes from) and "verdict_counting"
    lagf_bad = [w for w in W if lag[w].get("verdict_with_floors") == "OV late"]
    lagf_unt = [w for w in W if lag[w].get("verdict_with_floors") == "untestable"]
    lagf_ok = [w for w in W if lag[w].get("verdict_with_floors") == "consistent"]
    lagc_bad = [w for w in W if lag[w].get("verdict_counting") == "OV late"]
    lagq = R["ov_lag_test_quench"]
    lag_na = [w for w in W if not lag[w].get("applicable")]
    d_free, d_fix = DF["runs"]["free"], DF["runs"]["fixed"]
    cal_rej = [w for w in W if d_free[w]["at_calderoni"] >= CHI2_95]
    cal_ok = [w for w in W if w not in cal_rej]
    cal_str = ", ".join(f"{d_free[w]['at_calderoni']:.1f}" for w in cal_rej)
    h2_rej = [w for w in W if d_free[w]["at_H2_7e-10"] >= CHI2_95]
    h2_ok = [w for w in W if w not in h2_rej]
    # "the best fit mixes faster than Calderoni" must hold in every window where Calderoni is rejected
    for w in cal_rej:
        b = d_free[w]["best_log10_D"]
        if not (b is None or b > math.log10(d_free[w]["D_calderoni_m2_s"])):
            fail(f"{w}: Calderoni rejected but the best-fit D is below it")
    cal_islands_below = [w for w in cal_rej if any(x < math.log10(d_free[w]["D_calderoni_m2_s"])
                                                   for x in d_free[w]["region95"]["islands_log10"])]
    iv_only = {w: R["fits"][w]["IVonly_freeTBR"] for w in W}
    cov = R["synthetic_coverage"]
    cov_k = [c["log10_ktop"]["coverage95"] for c in cov.values()]
    cov_a = [c["log10_alpha"]["coverage95"] for c in cov.values()]
    n_k = sum(c["log10_ktop"]["reps"] for c in cov.values())
    pooled_k = sum(c["log10_ktop"]["coverage95"] * c["log10_ktop"]["reps"] for c in cov.values()) / n_k
    pooled_se = math.sqrt(0.95 * 0.05 / n_k)
    narrow = pooled_k < 0.95 - 2 * pooled_se
    frag6 = max(cov[f"run6/{v}"]["log10_ktop"]["fraction_error_over_1_decade"] for v in ("fixedTBR", "freeTBR"))
    run6_fragile = frag6 >= 0.1
    flat = R["synthetic_flat_direction"]
    dctl = DF["synthetic_control_run1_design"]
    dcal7, dcalc = dctl["D_7e-10_truth"], dctl["D_calderoni_run1_T_truth"]
    s3 = fits("run3_He")["physical"]["tbr_scale"]
    s3c = fits("run3_He")["profiles"]["log10_tbr_scale"]["ci95"]
    kt3fix = 10 ** R["fits"]["run3_He"]["fixedTBR"]["profiles"]["log10_ktop"]["mle"]
    kt3free = 10 ** fits("run3_He")["profiles"]["log10_ktop"]["mle"]

    # Qualitative sentences below are typed. Fail loudly if the data stop supporting them.
    pins_now = [(tb, w) for tb, dd in (("free", d_free), ("fixed", d_fix)) for w in W
                if dd[w]["region95"]["status"] == "identifiable"]
    if pins_now != [("fixed", "run3_He")]:
        fail(f"answer 6 says only run 3's He part with TBR fixed pins D; data say {pins_now}")
    if tr["null_degenerate"] or t13["null_degenerate"]:
        fail("answer 3 relies on usable He k_top tests")
    if any(R["fits"][w]["IVonly_freeTBR"]["profiles"]["log10_ktop"]["ci95"]["status"] == "identifiable" for w in W):
        fail("answer 1 says IV-only k_top is never identifiable")
    if lag_na != ["run6"]:
        fail(f"answer 2 says the OV test does not apply to run 6 only; not applicable in {lag_na}")

    def sig(t):
        return t["p_bootstrap"] < 0.05

    def p_text(t):
        n = t["bootstrap_reps"]
        floor = t["p_bootstrap"] <= 1 / (n + 1) + 1e-12
        return (f"bootstrap p = {t['p_bootstrap']:.3f} "
                + (f"(the smallest possible with {n} replicates, whose largest null LR was {t['null_LR_max']:.1f})"
                   if floor else f"({n} replicates)"))

    # D1: TBR range sensitivity and "TBR lands at 1" are computed, not typed
    wide_shift = {}
    for w in W:
        base, wide = fits(w), R["fits"][w]["freeTBR_wideTBR"]
        d = 0.0
        for n in ("log10_ktop", "log10_kwall", "log10_alpha"):
            pb, pw = base["profiles"][n], wide["profiles"][n]
            d = max(d, abs(pb["mle"] - pw["mle"]))
            for side in ("lo", "hi"):
                xb, xw = pb["ci95"][side], pw["ci95"][side]
                if (xb is None) != (xw is None):
                    d = float("inf")
                elif xb is not None:
                    d = max(d, abs(xb - xw))
            if pb["ci95"]["status"] != pw["ci95"]["status"]:
                d = float("inf")
        wide_shift[w] = d
    wide_same = [w for w in W if wide_shift[w] < 0.01]
    wide_moved = [w for w in W if w not in wide_same]
    tbr_at_1 = []
    for w in W:
        c = fits(w)["profiles"]["log10_tbr_scale"]["ci95"]
        if c["status"] == "identifiable" and c["lo"] <= 0 <= c["hi"]:
            tbr_at_1.append(w)
    tbr_not_1 = [w for w in W if w not in tbr_at_1]
    floor_bound = {w: fits(w)["floors_at_lower_bound"] for w in W}
    ov_collapsed = [w for w in W if "OV" in floor_bound[w]]
    iv_collapsed = [w for w in W if "IV" in floor_bound[w]]
    c_ov = R["synthetic_coverage_nonzero_OV_floor"]
    si = R["salt_inventory"]

    # ================= README.md: one screen. Every number computed, typed claims guarded.
    lp, cf = R["lab_point_check"], R["lab_style_cumfit"]
    lab_kt = {w: R["windows"][w]["lab_k_top_m_s"] for w in W}
    lab_kw = {w: R["windows"][w]["lab_k_wall_m_s"] for w in W}
    lab21 = lab_kt["run2"] / lab_kt["run1"]
    rL21 = R["pairwise_ktop_ratios_labTBR"]["run2/run1"]
    QS = R["quench"]["runs"]
    q3, q4 = QS["3"], QS["4"]
    QV = QV_["runs"]
    scan = R["quench"]["scan"]["points_tSIE_Bq_sigma"]
    hiT = R["quench"]["high_tSIE"]
    shapes = ("flat_below_300", "linear", "quadratic", "interp")
    SHAPE_LABEL = {"flat_below_300": f"held flat below tSIE {QV_['flat_below']:.0f}",
                   "linear": "straight-line fit", "quadratic": "quadratic fit", "interp": "as interpolated (run 6)"}

    def dtbr(r, key):
        return QV[r]["variants"][key]["measured_TBR"] / QS[r]["measured_TBR_published"] - 1

    hi_rng = {r: (min(dtbr(r, f"anchored_high_tSIE/{s_}") for s_ in shapes),
                  max(dtbr(r, f"anchored_high_tSIE/{s_}") for s_ in shapes)) for r in ("3", "4")}
    all_rng = {r: (min(dtbr(r, f"anchored/{s_}") for s_ in shapes),
                   max(dtbr(r, f"anchored/{s_}") for s_ in shapes)) for r in ("3", "4")}
    # guards for section 1
    for r, q in (("3", q3), ("4", q4)):
        v = QV[r]
        if not v["high_tSIE_first_bubbler_vials"] >= v["high_tSIE_vials"] - 1:
            fail(f"section 1 says the high-tSIE vials are the first bubbler vial (run {r})")
        if not q["published"]["diagnostics"]["mean_net_low_activity_high_tSIE_Bq"] < 0:
            fail(f"section 1 says run {r}'s high-tSIE low-activity vials average below zero as published")
        if not all(QV[r]["variants"][f"anchored_high_tSIE/{s_}"]["delta_Bq"] > 0 for s_ in shapes):
            fail(f"section 1 says every curve shape adds activity in the high-tSIE vials of run {r}")
        if not (hi_rng[r][1] <= all_rng[r][1] + 1e-12):
            fail("section 1 says correcting all vials is the upper end")
    for r in ("1", "2"):
        if not QS[r]["published"]["diagnostics"]["mean_net_low_activity_high_tSIE_Bq"] < 0:
            fail(f"section 1 says run {r} shows the same signature")
        if not all("NM-C" not in x for x in QS[r]["quench_set"]):
            fail(f"section 1 says run {r} used a different quench set")
    # run 2 note
    geo = R["geometry"]

    def alpha(kt, kw):
        return (geo["A_top_m2"] * kt + geo["A_wall_m2"] * kw) / geo["V_m3"]
    f1L, f2L = R["fits"]["run1"]["labTBR"]["physical"], R["fits"]["run2"]["labTBR"]["physical"]
    a_lab = alpha(lab_kt["run2"], lab_kw["run2"]) / alpha(lab_kt["run1"], lab_kw["run1"])
    a_fit = f2L["alpha_per_s"] / f1L["alpha_per_s"]
    k2_vs_lab1 = f2L["k_top_m_s"][0] / lab_kt["run1"]
    best = CH["best"]
    cnull = CH["null"]
    floor2 = f2L["noise_floor_Bq"]["IV"]
    floors_he = {w: R["fits"][w]["labTBR"]["physical"]["noise_floor_Bq"]["IV"] for w in ("run1", "run3_He")}
    if not floor2 > max(floors_he.values()):
        fail("run 2 note says run 2 has the largest IV scatter floor of the He runs")
    if not cnull["p_value_best"] < 0.05:
        fail("run 2 note says a change in release rate is supported")
    if not best["tau_after_day"] < best["tau_before_day"]:
        fail("run 2 note says release speeds up after the change")
    kafter_lab1 = best["k_top_after_m_s"] / lab_kt["run1"]
    ci_after = [x / lab_kt["run1"] for x in best["k_top_after_ci95_m_s"]]

    def pm(v, sd, d=2):
        return f"{v:.{d}f} ± {sd:.{d}f}"

    a("<!-- GENERATED FILE. Built by make_readme.py from the JSON files in this folder. Do not edit by hand. -->")
    a("")
    a("# BABY-1L: the run-6 tSIE background, applied to runs 3 and 4")
    a("")
    if DRAFT:
        a("> **DRAFT.** Built from low replicate counts (see APPENDIX.md, Reproduce). P-values, coverage and power "
          "are not final.")
        a("")
    a("The per-vial data are rebuilt from the public raw LSC files with each run's own `tritium_model.py` and match "
      "`processed_data.json` exactly. Every number below is computed by the scripts in this folder. Fits, controls "
      "and caveats are in APPENDIX.md.")
    a("")
    v3, v4 = QV["3"], QV["4"]
    a(f"Up to run 4 each vial has its file's blank subtracted, counted at tSIE {v3['blank_tSIE_range'][0]:.0f} to "
      f"{v4['blank_tSIE_range'][1]:.0f}. In runs 3 and 4, {v3['high_tSIE_vials'] + v4['high_tSIE_vials']} vials sit "
      f"above tSIE {hiT:.0f}, and {v3['high_tSIE_first_bubbler_vials'] + v4['high_tSIE_first_bubbler_vials']} of "
      "them are the first bubbler vial of a sample, whose water evaporates between samples. The run-6 blank scan "
      f"(`SLR_BLANK_SCAN.csv`) puts the background at those tSIE values well below the blank: {scan[0][1]:.2f} Bq "
      f"at tSIE {scan[0][0]:.0f}, {scan[3][1]:.2f} at {scan[3][0]:.0f}, {scan[-1][1]:.2f} at {scan[-1][0]:.0f}. From "
      "run 6 on, `tritium_model.py` subtracts `background_curve(tSIE)`. Applied to runs 3 and 4, with the curve "
      "scaled so it equals each file's own blank at the blank's tSIE:")
    a("")
    a("| Vials corrected | Curve below tSIE " + f"{QV_['flat_below']:.0f}" + " | Run 3 IV + OV (Bq) | Run 3 measured TBR | Run 3 OV soluble | "
      "Run 4 IV + OV (Bq) | Run 4 measured TBR | Run 4 OV soluble |")
    a("| --- | --- | --- | --- | --- | --- | --- | --- |")
    p3, p4 = q3["published"], q4["published"]
    a(f"| none (as published) | | {pm(p3['IV_plus_OV_Bq'], p3['sd']['IV_plus_OV_Bq'])} | "
      f"{q3['measured_TBR_published']:.2e} | {pct(p3['OV']['soluble_share'])} | "
      f"{pm(p4['IV_plus_OV_Bq'], p4['sd']['IV_plus_OV_Bq'])} | {q4['measured_TBR_published']:.2e} | "
      f"{pct(p4['OV']['soluble_share'])} |")
    for mode, lab_ in (("anchored_high_tSIE", f"above tSIE {hiT:.0f} only"), ("anchored", "all")):
        for s_ in shapes:
            x3, x4 = QV["3"]["variants"][f"{mode}/{s_}"], QV["4"]["variants"][f"{mode}/{s_}"]
            a(f"| {lab_} | {SHAPE_LABEL[s_]} | {pm(x3['IV_plus_OV_Bq'], x3['sd']['IV_plus_OV_Bq'])} | "
              f"{x3['measured_TBR']:.2e} ({dtbr('3', f'{mode}/{s_}'):+.1%}) | {pct(x3['OV']['soluble_share'])} | "
              f"{pm(x4['IV_plus_OV_Bq'], x4['sd']['IV_plus_OV_Bq'])} | {x4['measured_TBR']:.2e} "
              f"({dtbr('4', f'{mode}/{s_}'):+.1%}) | {pct(x4['OV']['soluble_share'])} |")
    a("")
    a(f"Correcting only the vials above tSIE {hiT:.0f} raises the measured TBR by {hi_rng['3'][0]:.1%} to "
      f"{hi_rng['3'][1]:.1%} in run 3 and {hi_rng['4'][0]:.1%} to {hi_rng['4'][1]:.1%} in run 4 across the four "
      "curve shapes. That is the part that holds up. Correcting every vial adds more, up to "
      f"{all_rng['3'][1]:.1%} and {all_rng['4'][1]:.1%}, but that depends on the curve between tSIE "
      f"{scan[0][0]:.0f} and {scan[1][0]:.0f}, where the scan's first two points drop by "
      f"{scan[0][1] - scan[1][1]:.2f} Bq (about {(scan[0][1] - scan[1][1]) / math.hypot(scan[0][2], scan[1][2]):.0f} "
      "sigma) and the ordinary vials sit just above their blanks. Read the all-vials rows as the upper end of a "
      "range. ± is counting noise (vials, blanks, scan points, "
      f"{QV_['n_mc']} draws) and does not include the choice of curve shape. With the published background, "
      f"low-activity vials above tSIE {hiT:.0f} average "
      f"{q3['published']['diagnostics']['mean_net_low_activity_high_tSIE_Bq']:+.3f} Bq (run 3) and "
      f"{q4['published']['diagnostics']['mean_net_low_activity_high_tSIE_Bq']:+.3f} Bq (run 4), below zero.")
    a("")
    a("Runs 1 and 2 show the same sign (low-activity vials above tSIE "
      f"{hiT:.0f} average {QS['1']['published']['diagnostics']['mean_net_low_activity_high_tSIE_Bq']:+.3f} and "
      f"{QS['2']['published']['diagnostics']['mean_net_low_activity_high_tSIE_Bq']:+.3f} Bq). They were counted "
      f"with the {QS['1']['quench_set'][0].split(': ')[-1]} quench set, not the "
      f"{R['quench']['scan']['quench_set'].split(': ')[-1]} set of the scan, so the run-6 curve cannot be carried "
      "over. A blank scan with that quench set would settle them.")
    a("")
    a("## A note on run 2")
    a("")
    a(f"The run-2 script sets k_top to {lab21:.1f} times run 1's value and k_wall to 0 (`analysis/tritium/"
      "tritium_model.py`, lines 223 to 224 at v0.5). Fitted with one k_top at the lab's TBR, run 2's k_top is "
      f"{rL21['ratio']:.2f} of run 1's fitted value and {k2_vs_lab1:.2f} of run 1's script value. On the total "
      f"release rate (top plus wall), the scripts imply {a_lab:.2f} and the fits give {a_fit:.2f}. But one k_top "
      f"does not describe run 2: its IV scatter floor is {floor2:.2f} Bq, against "
      f"{min(floors_he.values()):.2f} to {max(floors_he.values()):.2f} in the other He runs. Letting k_top change "
      f"once fits much better with the change at the day-{best['change_day']:.0f} sample (2ΔNLL "
      f"{best['two_delta_nll_gain']:.0f}, p = {cnull['p_value_best']:.3f}, the smallest possible with "
      f"{cnull['reps']} simulations. The {len(CH['per_change_day'])} days tried were picked after looking at the data, and the best is calibrated by simulation). The "
      f"release time goes from {best['tau_before_day']:.1f} to {best['tau_after_day']:.1f} days, and the later k_top "
      f"is {kafter_lab1:.2f} of run 1's script value (95%: {ci_after[0]:.2f} to {ci_after[1]:.2f}), "
      + ("consistent with" if ci_after[0] <= lab21 <= ci_after[1] else "just above" if lab21 < ci_after[0] else "just below")
      + f" the {lab21:.1f} in the script. The paper suggests the salt freezing mid-run as a possible explanation for "
      "run 2's slower release (sec. 3.3). A partly frozen salt could look like that. A full stop in release fits "
      "worse than one constant k_top, and the run's files do not say when the salt froze, so this is not tested. "
      "Details: APPENDIX.md, Run 2.")
    a("")
    a("## What this does not show")
    a("")
    a("- The 0D model lumps salt mixing, the gas-side film and surface kinetics into k_top. The intervals are "
      f"conditional on it, and synthetic tests put the nominal 95% intervals for k_top at {pct(pooled_k)} "
      "(APPENDIX.md, Controls).")
    dk = {w: abs(R["fits"][w]["labTBR_quench"]["physical"]["k_top_m_s"][0]
                 / R["fits"][w]["labTBR"]["physical"]["k_top_m_s"][0] - 1) for w in ("run3_He", "run4_1000ppm")}
    a(f"- Refitting runs 3 and 4 with every vial corrected (the upper end) moves k_top by {pct(dk['run3_He'])} and "
      f"{pct(dk['run4_1000ppm'])} (lab TBR). The background changes the collected totals and chemical form, not "
      "the release rates.")
    a("- These curves do not measure the salt diffusivity, and nothing here tests the paper's diffusion-limited "
      "reading (APPENDIX.md, Salt diffusivity).")
    a("- Runs 7 to 10 and the FLiBe-1L runs are not used. Paper: arXiv 2509.26174.")
    a("")
    (HERE / "README.md").write_text("\n".join(L))
    # ================= APPENDIX.md
    L.clear()
    a("<!-- GENERATED FILE. Built by make_readme.py. Do not edit by hand. -->")
    a("")
    a("# Appendix: fits, controls and reproduction")
    a("")
    if DRAFT:
        a("> **DRAFT.** Built from low replicate counts. P-values, coverage and power are not final.")
        a("")
    a("README.md has the two main results. This file has everything behind them and the rest of the "
      "identifiability study: which release-model parameters the BABY-1L data can pin down.")
    a("")
    a("## Answers, parameter by parameter")
    a("")
    alpha_ok = [w for w in W if status(fits(w), "log10_alpha") == "identifiable"]
    kw_two_sided = [w for w in W if status(fits(w), "log10_kwall") == "identifiable"]
    cn_iv_all = [R["fits"][w]["IVonly_freeTBR"]["fisher"]["physics_block"] for w in W]
    cn_iv_fin = [b["condition_number"] for b in cn_iv_all if b and math.isfinite(b["condition_number"])]
    cn_f = [R["fits"][w]["freeTBR"]["fisher"]["physics_block"]["condition_number"] for w in W[:4]]
    if any(status(fits(w), "log10_ktop") != "identifiable" for w in W):
        fail("answer 1 says IV + OV gives k_top a two-sided interval in every window")
    ivk = iv_only["run1"]["profiles"]["log10_ktop"]["ci95"]
    ivk_wide = R["fits"]["run1"]["IVonly_freeTBR_wideTBR"]["profiles"]["log10_ktop"]["ci95"]
    a("1. **The total release rate is identifiable. The split between top and wall is not, from the inner "
      "vessel alone.** The inner-vessel (IV) curve pins alpha = (A_top k_top + A_wall k_wall)/V, the release "
      f"time constant, in {names_list(alpha_ok)} (IV + OV, TBR free{', run 6 only softly, see Controls' if run6_fragile else ''}). Run 1: tau = {tau_ci(f1)} "
      "days, 95% profile interval. With IV data only and TBR free, k_top gets only an upper bound: the value it "
      f"takes if all loss goes through the top ({e(10 ** ivk['hi'])} m/s in run 1). Below that the profile is "
      "flat, and its lower edge is wherever the allowed TBR range ends "
      f"({e(10 ** ivk['lo'])} m/s for a factor of 3.16, {e(10 ** ivk_wide['lo'])} m/s for a factor of 10). "
      "That holds with TBR free. With TBR fixed, the IV total itself sets the top share, so the split then rests "
      "on the TBR assumption. Adding the OV stream gives k_top a two-sided interval in every window (condition "
      f"number of the physics block {min(cn_f):.0f} to {max(cn_f):.0f} in runs 1 to 4, against "
      f"{e(min(cn_iv_fin), 1) if cn_iv_fin else 'infinite'} to infinite with IV only).")
    head2 = ("k_top gets a two-sided interval only when the outer-vessel (OV) stream is fitted too, and the "
             "test of that part of the model is weak.")
    body2 = (f"On counting noise alone, OV tritium looks late (still arriving after the IV has stopped) in "
             f"{names_list(lagc_bad)}. But counting noise alone does not describe these data (answer 7). Once each "
             "window's fitted scatter floors are added, "
             + (f"{names_list(lagf_bad)} still {'shows' if len(lagf_bad) == 1 else 'show'} a late OV"
                + (f", with the OV floor at its lower bound in {names_list([w for w in lagf_bad if w in ov_collapsed])}"
                   if any(w in ov_collapsed for w in lagf_bad) else "") + ". " if lagf_bad else "no window shows a late OV. ")
             + (f"{names_list(lagf_unt)} {'is' if len(lagf_unt) == 1 else 'are'} untestable. " if lagf_unt else "")
             + (f"{names_list(lagf_ok)} {'is' if len(lagf_ok) == 1 else 'are'} consistent. " if lagf_ok else "")
             + "So a late OV, which the 0D model cannot produce, is suggested but not established.")
    a(f"2. **{head2}** The OV total sets the wall share. {body2} The test does not apply to run 6 (its k_top "
      f"steps). k_wall has a two-sided interval in {names_list(kw_two_sided)}. Read k_wall as an OV amplitude, "
      "not a wall property.")
    head3 = ("Run 2 released about half as fast as run 1, as the paper reports." if sig(tr) and r12["ci95"][1] < 0.7
             else "The He-run k_top comparison is not conclusive at this bootstrap size.")
    agree13 = not sig(t13)
    a(f"3. **{head3}** The paper links run 2's slower release to a heater malfunction that froze the salt mid-run "
      f"(arXiv 2509.26174, sec. 3.3). The fits put run 2's k_top at {r12['ratio']:.2f} times run 1's "
      f"(95%: {r12['ci95'][0]:.2f} to {r12['ci95'][1]:.2f}, TBR free). The run-2 script sets {lab21:.2f}. "
      "README.md (a note on run 2) and Run 2 below look at a change in rate during the run. A shared k_top for runs 1, 2 and 3 (He part) "
      f"{'is rejected' if sig(tr) else 'is not rejected'}: likelihood ratio {tr['LR']:.1f}, {p_text(tr)}, asymptotic "
      f"p = {tr['p_asymptotic']:.0e}. Runs 1 and 3 (He part) {'agree' if agree13 else 'differ'} on the lab's numbers "
      f"(ratio {r13['ratio']:.2f}, 95%: {r13['ci95'][0]:.2f} to {r13['ci95'][1]:.2f}, LR {t13['LR']:.1f}, "
      f"p = {t13['p_bootstrap']:.2f}).")
    dT = INC["runs"]["4"]["temperature_salt_C"] - INC["runs"]["1"]["temperature_salt_C"]
    lab41 = lab_kt["run4_1000ppm"] / lab_kt["run1"]
    ds6 = fits("run6")
    n6, p6 = ds6["n_points"], len(ds6["theta"])
    if not n6 < 2 * p6:
        fail("answer 4 says run 6 has few points per free parameter")
    a(f"4. **The H2 effect the paper reports (sec. 3.4) shows up the same way.** The `processed_data.json` values "
      f"already put run 4 (1000 ppm H2) at {lab41:.1f} times run 1. The fits give {r41['ratio']:.0f} times "
      f"(95%: {r41['ci95'][0]:.0f} to {r41['ci95'][1]:.0f}), so this adds an interval and nothing else. There is "
      f"one run per H2 condition, and gas is confounded with run order, irradiation schedule and a {dT:.0f} C "
      f"difference in logged salt temperature. Run 6 has {n6} data points for {p6} free parameters (k_top before "
      "and after its switch, k_wall, TBR scale and two scatter floors), so its fitted values are not used as "
      "results"
      + (". The run 4 against run 6 test is not usable (its bootstrap null is degenerate)." if tH2["null_degenerate"]
         else f". A shared k_top for runs 4 and 6 gives LR {tH2['LR']:.1f}, {p_text(tH2)}, so the two H2 runs "
              f"{'differ' if sig(tH2) else 'are not told apart'} at this bootstrap size."))
    f3L = R["fits"]["run3_He"]["labTBR"]
    kt3lab = f3L["physical"]["k_top_m_s"][0]
    c3L = f3L["profiles"]["log10_ktop"]["ci95"]
    pen3 = lp["run3_He"]["labTBR"]["two_delta_nll_at_lab_ktop"]
    tbr3_used, tbr3_meas = INC["runs"]["3"]["TBR_used_in_model"], INC["runs"]["3"]["measured_TBR"]
    tbr3_free = s3 * tbr3_meas
    if not pen3 < CHI2_95:
        fail("answer 5 says the lab's run 3 k_top sits inside the interval at the lab's TBR")
    if not abs(tbr3_free / tbr3_used - 1) < 0.15:
        fail("answer 5 says the free-TBR fit lands near the lab's modelled TBR for run 3")
    a(f"5. **At the lab's own production assumption, run 3's fit agrees with `processed_data.json`.** For run 3 "
      f"the lab's model uses the modelled TBR ({tbr3_used:.2e}), not the measured {tbr3_meas:.2e}, which the paper "
      f"puts down to extra release (sec. 3.2 and 3.4). At that TBR the He part gives k_top = {e(kt3lab)} m/s (95%: "
      f"{e(10 ** c3L['lo'])} to {e(10 ** c3L['hi'])}), and the `processed_data.json` value {e(lab_kt['run3_He'])} "
      f"costs 2ΔNLL {pen3:.1f}. Letting TBR float lands at {s3:.2f} of the measured value (95%: "
      f"{10 ** s3c['lo']:.2f} to {10 ** s3c['hi']:.2f}), a TBR of {tbr3_free:.2e}, within "
      f"{pct(abs(tbr3_free / tbr3_used - 1))} of the modelled one. Fixing TBR at the measured value instead would "
      f"give k_top = {e(kt3fix)} m/s, {kt3free / kt3fix:.1f} times lower than the free fit, which is why the lab's "
      "choice matters for this run.")
    a("6. **These curves do not measure the salt diffusivity.** The 1D column used here puts all of the wall "
      "loss on its bottom face, uses the 1 L model height, and has low power: on run 1's design, synthetic data "
      f"made at the Calderoni value were told apart from a well-mixed salt in {pct(dcalc['well_mixed_rejected_fraction'])} "
      f"of {dcalc['reps']} replicates. No window pins D unless TBR is fixed (run 3's He part), and that interval "
      "comes from the TBR assumption. Nothing here tests the paper's Sherwood analysis (sec. 3.3) or its "
      "diffusion-limited reading. The per-window profiles are under Salt diffusivity.")
    a(f"7. **Counting noise is not what limits these fits. Model error is.** With LSC counting noise alone, "
      f"chi2 per degree of freedom runs from {min(g['chi2_per_dof'] for g in gof.values()):.0f} to "
      f"{max(g['chi2_per_dof'] for g in gof.values()):.0f}. Every interval here includes a fitted per-stream "
      f"scatter floor on top of counting noise. In synthetic tests ({R['reps']} replicates per window and TBR "
      f"mode) the nominal 95% intervals covered {pct(min(cov_k))} to {pct(max(cov_k))} of the time for k_top, "
      f"{pct(pooled_k, 1)} pooled (binomial SE {100 * pooled_se:.1f} points)"
      + (f", so read them as roughly {100 * pooled_k:.0f}% intervals." if narrow else
         ", consistent with 95% at this replicate count."))
    a("")
    a("Every number above comes from the scripts in this folder, listed at the end.")
    a("")

    a("## Fitted values per run window")
    a("")
    a("IV + OV fitted, TBR free as a nuisance scale. 95% profile-likelihood intervals. \"open\" means the data "
      "bound that side only up to a physical or prior limit. The lab's values are from each run's "
      "`processed_data.json` at the pinned ref. Run 6 is shown for completeness only (answer 4).")
    a("")
    a("| Window | Gas | tau (d) | k_top (m/s) | lab k_top | k_wall (m/s) | lab k_wall | TBR / measured | IV floor (Bq) | OV floor (Bq) |")
    a("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for w in W:
        b = fits(w)
        ph = b["physical"]
        fl = ph["noise_floor_Bq"]
        a(f"| {NAME[w]} | {GAS[w]} | {tau_ci(b)} | {ci_lin(b, 'log10_ktop')} | {e(R['windows'][w]['lab_k_top_m_s'])} | "
          f"{ci_lin(b, 'log10_kwall')} | {e(R['windows'][w]['lab_k_wall_m_s']) if R['windows'][w]['lab_k_wall_m_s'] else '0'} | {ci_lin(b, 'log10_tbr_scale', plain=True)} | "
          f"{fl.get('IV', 0):.2f}{' (at bound)' if 'IV' in b['floors_at_lower_bound'] else ''} | "
          f"{fl.get('OV', 0):.2f}{' (at bound)' if 'OV' in b['floors_at_lower_bound'] else ''} |")
    a("")
    a(f"\"(at bound)\": the fitted floor sits on its lower limit, {10 ** -4:.0e} Bq, so it is effectively zero. "
      f"That happens for the OV floor in {names_list(ov_collapsed)}"
      + (f" and the IV floor in {names_list(iv_collapsed)}" if iv_collapsed else "")
      + " (TBR free). Controls below test whether that collapse is real.")
    a("")
    a("With TBR fixed at the value the lab's model uses (`TBR_used_in_model`: measured for runs 1, 2, 4 and 6, "
      "modelled for run 3), and 2ΔNLL at the lab's k_top with k_wall re-fitted:")
    a("")
    a("| Window | tau (d) | k_top (m/s) | k_wall (m/s) | 2ΔNLL at lab k_top | lab-style cumulative fit k_top |")
    a("| --- | --- | --- | --- | --- | --- |")
    for w in W:
        b = fits(w, "labTBR")
        a(f"| {NAME[w]} | {tau_ci(b)} | {ci_lin(b, 'log10_ktop')} | {ci_lin(b, 'log10_kwall')} | "
          f"{lp[w]['labTBR']['two_delta_nll_at_lab_ktop']:.1f} | {e(cf[w]['k_top_m_s'])} |")
    a("")
    a("The lab-style column is an unweighted least-squares fit to the cumulative IV and OV curves at the lab's "
      "TBR, with no noise model, as a check that the comparison does not hinge on the noise model.")
    a("")
    a("With TBR fixed at the lab's measured value instead:")
    a("")
    a("| Window | tau (d) | k_top (m/s) | k_wall (m/s) |")
    a("| --- | --- | --- | --- |")
    for w in W:
        b = fits(w, "fixedTBR")
        a(f"| {NAME[w]} | {tau_ci(b)} | {ci_lin(b, 'log10_ktop')} | {ci_lin(b, 'log10_kwall')} |")
    a("")
    a("")
    a("Runs 3 and 4 with the run-6 background (README section 1), at the lab's TBR and with TBR free, and 2ΔNLL "
      "at the lab's k_top:")
    a("")
    a("| Window | TBR | k_top (m/s) | k_wall (m/s) | 2ΔNLL at lab k_top |")
    a("| --- | --- | --- | --- | --- |")
    for w in ("run3_He", "run4_1000ppm"):
        for v, lab_ in (("labTBR_quench", "lab"), ("freeTBR_quench", "free")):
            b = R["fits"][w][v]
            a(f"| {NAME[w]} | {lab_} | {ci_lin(b, 'log10_ktop')} | {ci_lin(b, 'log10_kwall')} | "
              f"{lp[w][v]['two_delta_nll_at_lab_ktop']:.1f} |")
    a("")
    b6 = fits("run6")
    a(f"Run 6's k_top after its switch to He: {status(b6, 'log10_ktop_p2')}. Only one IV sample falls after the "
      "switch. Run 6's measured TBR carries tritium that left before the switch, so the He value cannot be "
      "tested against runs 1 to 3.")
    a("")
    for w in ("run3_full_step", "run4_full_step"):
        b = R["fits"][w]["freeTBR"]
        a(f"{NAME[w]} with one k_top per gas: k_top before the switch is {status(b, 'log10_ktop')}, after the switch "
          f"{status(b, 'log10_ktop_p2')}. The scatter floor rises to {b['physical']['noise_floor_Bq']['IV']:.2f} Bq (IV) "
          f"and {b['physical']['noise_floor_Bq']['OV']:.2f} Bq (OV). A single salt reservoir cannot produce the "
          "late releases after the gas switch.")
    a("")
    a("![k per run](fig_k_per_run.png)")
    a("")

    a("## Identifiability, parameter by parameter")
    a("")
    a("Status from the profile likelihood (95%). Edges are found by root-finding, not grid interpolation. A side "
      "is open when the TBR scale is held at its allowed range there (factor 3.16, a prior, so no information from "
      "the data), or when the profile dips back below the 95% line further out. The physical limit k_wall = 0 "
      "does count: an upper bound on k_top from \"all loss through the top\" is real information.")
    a("")
    a("| Window | Data used | k_top | k_wall | alpha (total rate) | TBR scale |")
    a("| --- | --- | --- | --- | --- | --- |")
    for w in W:
        for v, lab in (("fixedTBR", "IV+OV, TBR fixed"), ("freeTBR", "IV+OV, TBR free"), ("IVonly_freeTBR", "IV only, TBR free")):
            b = R["fits"][w][v]
            tb = status(b, "log10_tbr_scale") if "log10_tbr_scale" in b["profiles"] else "fixed"
            a(f"| {NAME[w]} | {lab} | {status(b, 'log10_ktop')} | {status(b, 'log10_kwall')} | "
              f"{status(b, 'log10_alpha')} | {tb} |")
    a("")
    wv = [(w, v, n, x) for w in W for v in ("fixedTBR", "freeTBR") for n, x in R["fits"][w][v]["wald_vs_profile"].items()
          if x["ratio"] is not None and n == "log10_ktop"]
    rat = [x["ratio"] for *_, x in wv]
    cn_iv = [R["fits"][w]["IVonly_freeTBR"]["fisher"]["physics_block"]["condition_number"] for w in W
             if R["fits"][w]["IVonly_freeTBR"]["fisher"]["physics_block"]]
    cn_full = [R["fits"][w]["freeTBR"]["fisher"]["physics_block"]["condition_number"] for w in W[:4]]
    a(f"Fisher information against the profiles: where k_top is identifiable, the Wald 95% half-width is "
      f"{min(rat):.2f} to {max(rat):.2f} times the profile half-width. Fisher is a fair shortcut there. It is not "
      f"one for k_wall, whose profiles are lopsided. The Fisher matrix does flag the flat direction. Its condition "
      f"number for (k_top, k_wall, TBR) is {min(cn_iv):.0e} to {max(cn_iv):.0e} with IV only, against "
      f"{min(cn_full):.0f} to {max(cn_full):.0f} with IV + OV (runs 1 to 4).")
    a("")
    a("![profiles](fig_profiles.png)")
    a("")

    a("## Transferability tests")
    a("")
    a("Likelihood ratio (LR) of \"one shared value\" against \"one value per run\". Every other parameter, "
      "including TBR scale and scatter floors, stays per run. The p-value comes from a parametric bootstrap under "
      f"the shared-value fit ({R['boot_reps']} replicates, so the smallest reportable p is "
      f"{1 / (R['boot_reps'] + 1):.3f}). Power is the rejection rate when the second run's value is really 1.5 or "
      "2 times the first.")
    a("")
    a("| Runs | Shared | LR | p (bootstrap) | p (chi2) | null LR median / 95% / max | power x1.5 | power x2 | usable |")
    a("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for k in ("transferability_He_ktop", "transferability_run1_run2_ktop", "transferability_run1_run3He_ktop",
              "transferability_He_kwall", "transferability_H2_run4_run6_ktop"):
        t = R[k]
        pw = t["power_if_second_run_differs_by"]
        a(f"| {', '.join(NAME[x] for x in t['runs'])} | {t['shared'][0].replace('log10_', '')} | {t['LR']:.1f} | "
          f"{t['p_bootstrap']:.3f} | {t['p_asymptotic']:.1e} | {t['null_LR_median']:.1f} / "
          f"{t['bootstrap_crit95']:.1f} / {t['null_LR_max']:.1f} | {pct(pw['x1.5'])} | {pct(pw['x2.0'])} | "
          f"{'no, null degenerate' if t['null_degenerate'] else 'yes'} |")
    a("")
    tkeys = [k for k in R if k.startswith("transferability_")]
    degen = [k for k in tkeys if R[k]["null_degenerate"]]
    a("A test is marked not usable when half or more of its null replicates give an LR near zero. The fits then "
      "sit on a flat or bounded direction, and the bootstrap cannot calibrate the statistic. "
      + ("That applies to " + "; ".join(", ".join(NAME[x] for x in R[k]["runs"]) + " (" +
                                         R[k]["shared"][0].replace("log10_", "") + ")" for k in degen) + ". "
         if degen else "No test here is in that state. ")
      + "Real and bootstrap LRs use the same optimiser as the real fits (the 8-point start grid plus Powell "
      "polish). The separate fits inside each test reproduce the stored per-run fits, which the README build "
      f"checks. Power uses {tr['power_reps']} replicates per factor (binomial SE up to "
      f"{100 * math.sqrt(0.25 / tr['power_reps']):.0f} points). The run 1 against run 3 test has "
      f"{pct(t13['power_if_second_run_differs_by']['x2.0'])} power at a factor of 2 and "
      f"{'finds no difference' if not sig(t13) else 'finds a difference'}. But run 3's "
      "value moves by more than that factor with the TBR assumption (answer 5).")
    a("")

    a("## Is the OV stream consistent with the 0D model?")
    a("")
    a("With constant k_top and k_wall, IV and OV are both integrals of the same salt concentration, so "
      "OV(t)/IV(t) is constant. Test at a common time t*: the first OV sample time (other than the last) by "
      "which the IV has collected 80% of its total. Compared: each stream's share collected after t*. 95% "
      "intervals two ways: counting noise only, and counting noise plus the window's fitted per-stream scatter "
      "floors (TBR free fit). The fits need those floors, so the floor-included verdict is the one concluded from. "
      "\"Untestable\" means the interval is wider than 1 or reaches outside [-1, 1], the range a share difference "
      "can take.")
    a("")
    a("| Window | Data | t* (day) | IV share after | OV share after | OV minus IV | 95%, counting only | "
      "Verdict, counting only | 95%, with floors | Verdict, with floors |")
    a("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for w in W:
        for lab, src in (("as published", lag), ("run-6 background", lagq)):
            if w not in src:
                continue
            x = src[w]
            if not x.get("applicable"):
                a(f"| {NAME[w]} | {lab} | | | | | | | | not applicable: {x['reason']} |")
                continue
            a(f"| {NAME[w]} | {lab} | {x['t_star_day']:.1f} | {pct(x['IV_share_after'])} | "
              f"{pct(x['OV_share_after'])} | {x['OV_minus_IV']:+.2f} | {x['ci95_counting'][0]:+.2f} to "
              f"{x['ci95_counting'][1]:+.2f} | {x['verdict_counting']} | {x['ci95_with_floors'][0]:+.2f} to "
              f"{x['ci95_with_floors'][1]:+.2f} | {x['verdict_with_floors']} |")
    a("")
    a("A delayed OV is what permeation through the crucible wall or holdup in the outer vessel would look like. "
      "The 0D model has neither. "
      + ("With floors included it shows in " + names_list(lagf_bad) + ", " +
         ("a lead, not a finding, given the collapsed OV floor there." if any(w in ov_collapsed for w in lagf_bad)
          else "but in no other window, so it is a lead, not a finding.") if lagf_bad else
         "With floors included no window shows it, so this is a lead, not a finding."))
    a("")

    a("## Chemical form by stream")
    a("")
    a("Share of each stream's whole-run total that was water soluble in the LSC split (roughly HTO/TF, against "
      "HT after the furnace). The 0D model ignores this.")
    a("")
    a("| Run | Gas | IV soluble share | OV soluble share | IV / OV soluble share, run-6 background |")
    a("| --- | --- | --- | --- | --- |")
    for run, w in (("1", "run1"), ("2", "run2"), ("3", "run3_He"), ("4", "run4_1000ppm"), ("6", "run6")):
        st = INC["runs"][run]["streams"]
        sh = {k: sum(st[k]["inc_soluble_Bq"]) / sum(st[k]["inc_total_Bq"]) for k in ("IV", "OV")}
        gas = INC["runs"][run]["gas_initial"] + (f", then {INC['runs'][run]['gas_switched_to']}"
                                                 if INC["runs"][run]["gas_switched_to"] else "")
        qa = R["quench"]["runs"][run].get("anchored")
        qcol = (f"{pct(qa['IV']['soluble_share'])} / {pct(qa['OV']['soluble_share'])}" if qa else
                "same (the lab's own curve)" if run == "6" else "not corrected (other quench set)")
        a(f"| {run} | {gas} | {pct(sh['IV'])} | {pct(sh['OV'])} | {qcol} |")
    a("")

    a("## Run-6 background: method and checks")
    a("")
    a("`quench.py` executes each run's `tritium_model.py` unchanged and takes `background_curve` from run 6's own "
      "namespace (`build_background_curve_from_file`, a linear interpolation of the 7 blanks in "
      f"`SLR_BLANK_SCAN.csv`, 2 to 11 mL water, quench set {R['quench']['scan']['quench_set'].split(': ')[-1]}). For a "
      "vial at tSIE t counted in a file whose blank sits at tSIE t_b and reads B, the background is "
      "B x curve(t) / curve(t_b): the scan's shape, the level of the blank counted with that vial. \"As written\" "
      "is curve(t) alone, as run 6 does it. Both keep the toolbox's clamp of negative net vials to zero. After the "
      "correction few vials go negative, so the clamp barely matters (run 3 IV + OV "
      f"{QS['3']['anchored']['IV_plus_OV_Bq']:.2f} Bq clamped, {QS['3']['anchored_no_clamp']['IV_plus_OV_Bq']:.2f} not, "
      "every vial corrected with the curve as interpolated).")
    a("")
    a("Checks, and what they do and do not show. (1) For every run, `quench.py` recomputes the lab's published "
      "net activity of every vial from the gross count and the background the lab subtracted, to 1e-9 Bq, and its "
      "per-sample totals equal `increments.json`, which `build_increments.py` matches to each "
      "`processed_data.json` to 1e-9 Bq. Those are the real validations. (2) The curve read from the scan file "
      "equals the lab's own `background_curve` object to 1e-12 Bq. (3) On run 6, subtracting that curve "
      "reproduces the published nets exactly. That last check is close to tautological, because run 6's "
      "published nets are that same curve subtracted. It confirms the plumbing, not the method.")
    a("")
    a("The curve between the blanks (tSIE "
      f"{QV['3']['blank_tSIE_range'][0]:.0f} to {QV['4']['blank_tSIE_range'][1]:.0f}) and ordinary vials (median "
      f"{QV['3']['other_vial_tSIE_range'][1]:.0f} in run 3) rests on the scan's two lowest points. How much each "
      "curve shape adds, split by where it comes from:")
    a("")
    a("| Run | Curve below tSIE " + f"{QV_['flat_below']:.0f}" + " | Added, vials above tSIE " + f"{hiT:.0f}" + " (Bq) | Added, all vials (Bq) | Negative vials, all corrected | Mean net, low-activity vials above tSIE " + f"{hiT:.0f}" + " (Bq) |")
    a("| --- | --- | --- | --- | --- | --- |")
    for r in ("3", "4"):
        for s_ in shapes:
            hv, av = QV[r]["variants"][f"anchored_high_tSIE/{s_}"], QV[r]["variants"][f"anchored/{s_}"]
            a(f"| {r} | {SHAPE_LABEL[s_]} | {hv['delta_Bq']:+.2f} | {av['delta_Bq']:+.2f} | "
              f"{av['diagnostics']['negative_vials']} | {av['diagnostics']['mean_net_low_activity_high_tSIE_Bq']:+.3f} |")
    a("")
    a("| Run | Quench set | Vials above tSIE " + f"{hiT:.0f}" + " | Blank minus curve (Bq, median) | Negative vials, published / run-6 | "
      "Mean net, low-activity vials above / below tSIE " + f"{hiT:.0f}" + ", published (Bq) | same, run-6 background (Bq) | Measured TBR, curve as written |")
    a("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in ("1", "2", "3", "4", "6"):
        q = QS[r]
        dp = q["published"]["diagnostics"]
        dq = q["anchored"]["diagnostics"] if "anchored" in q else None
        bm = q["blank_minus_curve_at_blank_tSIE_Bq"]
        bm_s = "" if bm is None else f"{bm['median']:+.3f}"
        neg = f"{dp['negative_vials']} / {dq['negative_vials']}" if dq else f"{dp['negative_vials']} / -"
        pub = (f"{dp['mean_net_low_activity_high_tSIE_Bq']:+.3f} / "
               f"{dp['mean_net_low_activity_low_tSIE_Bq']:+.3f}")
        cor = (f"{dq['mean_net_low_activity_high_tSIE_Bq']:+.3f} / {dq['mean_net_low_activity_low_tSIE_Bq']:+.3f}"
               if dq else "-")
        lit = f"{q['literal']['measured_TBR']:.2e}" if "literal" in q else "-"
        a(f"| {r} | {q['quench_set'][0].split(': ')[-1]} | {q['high_tSIE_vials']} of {dp['vials']} | {bm_s} | "
          f"{neg} | {pub} | {cor} | {lit} |")
    a("")
    a("Low-activity means gross below 0.6 Bq. Those vials still carry real tritium, so their mean need not be zero, "
      "but it cannot be negative. The first bubbler vial is the one whose water evaporates between samples (run "
      "2's `general.json`: \"Bubbler water in vials 2-X-1 replenished periodically to compensate for "
      "evaporation\"). Less water means less quench and higher tSIE, which is what the blank scan varies. "
      f"{QS['3']['vials_beyond_scan_tSIE']} run 3 vials lie beyond the scan's highest tSIE and use the linear "
      "extrapolation run 6 uses. Holding the curve flat there instead changes run 3's IV + OV total by "
      f"{QS['3']['anchored_flat_extrapolation_IV_plus_OV_Bq'] - QS['3']['anchored']['IV_plus_OV_Bq']:+.2f} Bq. "
      "Measured TBR scales with IV + OV, as in each run's `tritium_model.py` (T_produced / T_consumed). The ± "
      "values do not include the change of the scan's level between its count and runs 3 and 4, which the "
      "\"as written\" column brackets.")
    a("")
    a("## Run 2")
    a("")
    a("One k_top for the whole run against one that changes at a sample time (lab TBR, k_wall and floors "
      "re-fitted). The change days were picked after looking at the data, so the best one is calibrated by "
      f"simulation: {cnull['reps']} data sets from the one-k_top fit, each refitted at every change day, best gain "
      f"kept. Their 95th percentile is {cnull['gain_quantiles']['0.95']:.1f}.")
    a("")
    a("| Change at day | 2ΔNLL gain | tau before (d) | tau after (d) | k_top after (m/s) | IV floor (Bq) |")
    a("| --- | --- | --- | --- | --- | --- |")
    for r_ in CH["per_change_day"]:
        ci_ = r_["k_top_after_ci95_m_s"]
        a(f"| {r_['change_day']:.1f} | {r_['two_delta_nll_gain']:.1f} | {r_['tau_before_day']:.1f} | "
          f"{r_['tau_after_day']:.1f} | {e(r_['k_top_after_m_s'])}"
          + (f" ({e(ci_[0])} to {e(ci_[1])})" if ci_ else f" ({r_['k_top_after_status']})")
          + f" | {r_['floors_Bq']['IV']:.2f} |")
    a("")
    a(f"One k_top: {e(CH['one_k_top']['k_top_m_s'])} m/s, IV floor {floor2:.2f} Bq. The gain peaks sharply at one "
      "sample, so the change is located to within one sampling interval, not dated. The paper's freeze is not "
      "in the run's files. Ratios in README.md use run 1's lab-TBR fit "
      f"({e(f1L['k_top_m_s'][0])} m/s) or the run-1 script value ({e(lab_kt['run1'])} m/s).")
    a("")
    a("## Salt diffusivity")
    a("")
    a(f"1D salt column (`fv1d.py`), height {R['geometry']['height_m'] * 100:.2f} cm, Robin release at the top, "
      "wall loss on the bottom face, exact in time. D fixed on a grid from 1e-10 to 1e-5 m2/s, all else refit. "
      "The right end of each curve is the well-mixed 0D model. Checks: the solver matches the 0D model at large D "
      f"to {DF['fv_large_D_vs_0D_max_rel_gap']:.0e}. It matches FESTIM (a 1D FESTIM model with a Robin top, 50 cells, not included in the share bundle) to "
      f"{pct(FC['cases'][0]['max_rel_gap_fv_vs_festim'], 2)} at D = {e(FC['cases'][0]['D_m2_s'], 1)} and "
      f"{pct(FC['cases'][1]['max_rel_gap_fv_vs_festim'], 2)} at D = {e(FC['cases'][1]['D_m2_s'], 1)} m2/s. At those "
      f"two D values a stagnant salt changes the run 1 curve by {pct(FC['cases'][0]['max_rel_gap_0d_vs_festim'])} "
      f"and {pct(FC['cases'][1]['max_rel_gap_0d_vs_festim'])} against the well-mixed model, so the test has "
      f"something to see. Grid: 60 against 240 cells differ by {pct(DF['fv_grid_n60_vs_n240_at_7e-10'], 2)} at "
      "7e-10 m2/s.")
    a("")
    a(f"Column height caveat: the column here is the lab's 1 L model ({100 * si['height_model_m']:.2f} cm). The "
      f"logged salt charge ({si['salt_mass_logged_kg']:.2f} kg at {si['density_g_cm3']:.4f} g/cm3) fills about "
      f"{100 * si['height_best_estimate_m']:.2f} cm. Diffusive times scale with height squared, so every D "
      f"threshold in this section would move up by a factor of about {si['diffusivity_scale_height_squared']:.2f} "
      "with the logged charge (see Data and noise model).")
    a("")
    a("Reference values: Calderoni's T in FLiBe, D = 9.3e-7 exp(-42 kJ/mol / RT) at each run's logged salt "
      "temperature (Calderoni et al. 2008, doi:10.1016/j.fusengdes.2008.05.016), and 7e-10 m2/s for H2 near 900 K "
      "(Nakamura, Fukada and Nishiumi 2015, H2 in Flibe, as compiled in h-transport-materials: 7.1e-10 at 900 K, "
      "a short extrapolation of a fit measured at 773 to 873 K). Both are FLiBe numbers. We found no measured value for "
      "ClLiF.")
    a("")
    a("| Window | TBR | 95% region for D (m2/s) | 2dNLL at Calderoni | 2dNLL at 7e-10 | 2dNLL well mixed |")
    a("| --- | --- | --- | --- | --- | --- |")
    for tb, dd in (("free", d_free), ("fixed", d_fix)):
        for w in W:
            x = dd[w]
            g = x["region95"]
            reg = (f"above {e(10 ** g['lo_log10'], 1)}" if g["status"] == "lower bound only" else
                   f"{e(10 ** g['lo_log10'], 1)} to {e(10 ** g['hi_log10'], 1)}" if g["status"] == "identifiable" else
                   g["status"])
            if g["islands_log10"]:
                reg += f", plus isolated grid points at log10 D = {', '.join(str(v) for v in g['islands_log10'])}"
            a(f"| {NAME[w]} | {tb} | {reg} | {x['at_calderoni']:.1f} | {x['at_H2_7e-10']:.1f} | "
              f"{x['two_delta_nll_well_mixed']:.1f} |")
    a("")
    fp4 = d_free["run4_1000ppm"]["fit_along_profile"]
    isl = d_free["run4_1000ppm"]["region95"]["islands_log10"]
    if isl:
        fi = fp4[str(isl[0])]
        a(f"Run 4's isolated points near log10 D = {isl[0]} are degenerate fits. k_top sits at its upper limit "
          f"(a perfectly absorbing top), TBR is {10 ** fi['log10_tbr_scale']:.1f} times the measured value, and the IV "
          f"scatter floor grows to {10 ** fi['log10_sig_IV']:.1f} Bq against "
          f"{fits('run4_1000ppm')['physical']['noise_floor_Bq']['IV']:.2f} Bq well mixed. They vanish with TBR fixed.")
    pins = [(tb, w) for tb, dd in (("free", d_free), ("fixed", d_fix)) for w in W
            if dd[w]["region95"]["status"] == "identifiable"]
    a("Windows where the 95% set for D is closed on both sides: "
      + (", ".join(f"{NAME[w]} with TBR {tb}" for tb, w in pins) or "none") + ".")
    if ("fixed", "run3_He") in pins:
        a("For run 3's He part that happens only with TBR fixed. The slow tail that a fixed TBR forces is being "
          "explained by slow diffusion. That is the TBR assumption talking, not a measurement of D.")
    a("The 1D rung puts all the wall loss on the bottom face. A real crucible loses through its side wall too, "
      "which a 2D model would carry. The minimum salt temperature is the lab's logged `temperature_salt`, which "
      "the lab notes is a minimum.")
    a("")

    a("## Controls: does the method recover known truth?")
    a("")
    a("Synthetic data use each window's real sampling times, irradiation schedule, counting covariance and "
      "fitted scatter floors, with the fitted parameters as truth. Fits start from a fixed grid, never from "
      "the truth.")
    a("")
    cov_w = [c["log10_kwall"]["coverage95"] for c in cov.values()]
    npts = [fits(w)["n_points"] for w in W]
    a(f"- Synthetic fits use the same optimiser as the real fits: the 8-point start grid plus Powell polish.")
    a(f"- Coverage of nominal 95% profile intervals ({R['reps']} replicates per window and TBR mode, binomial SE "
      f"about {100 * math.sqrt(0.95 * 0.05 / R['reps']):.1f} points each): k_top {pct(min(cov_k))} to "
      f"{pct(max(cov_k))} ({pct(pooled_k, 1)} pooled), k_wall {pct(min(cov_w))} to {pct(max(cov_w))}, alpha "
      f"{pct(min(cov_a))} to {pct(max(cov_a))}. "
      + (f"Intervals run narrow. With {min(npts)} to {max(npts)} points per window, maximum likelihood "
         f"underestimates the scatter floors. Read the intervals as roughly {100 * pooled_k:.0f}%. "
         if narrow else "Pooled k_top coverage is consistent with 95% at this replicate count. ")
      + "Coverage is checked at the truth, so it does not test how the interval edges were located.")
    ovc = c_ov["cases"]
    ovc_k = [c["log10_ktop"]["coverage95"] for c in ovc.values()]
    ovc_coll = [c["fitted_floor_at_lower_bound_fraction"]["OV"] for c in ovc.values()]
    base_coll = [cov[k]["fitted_floor_at_lower_bound_fraction"]["OV"] for k in ovc]
    a(f"- OV floor collapse. With the collapsed OV floor as truth (as fitted), the refit OV floor lands on its "
      f"bound in {pct(min(base_coll))} to {pct(max(base_coll))} of replicates for "
      f"{', '.join(NAME[k.split('/')[0]] + (' TBR fixed' if k.endswith('fixedTBR') else ' TBR free') for k in ovc)}. With a non-zero OV floor as truth instead ({c_ov['truth_sig_OV_Bq']:.2f} Bq, the median "
      f"of the non-collapsed floors in {names_list(c_ov['from_windows'])}), it collapses in {pct(min(ovc_coll))} to "
      f"{pct(max(ovc_coll))}, and k_top coverage is {pct(min(ovc_k))} to {pct(max(ovc_k))} ({R['reps']} replicates "
      "each). "
      + ("So a collapse in the real fits is evidence that the OV scatter there really is small, not an optimiser "
         "artefact." if max(ovc_coll) <= 0.2 else
         "So the fit collapses the OV floor even when it is not zero; read the OV floors as unreliable."))
    c14 = {k: c for k, c in cov.items() if not k.startswith("run6")}
    r6 = INC["runs"]["6"]
    t6_first_h = r6["streams"]["IV"]["times_day"][0] * 24
    t6_irr_h = sum(b - a_ for a_, b in r6["irradiations_s"]) / 3600
    t6_half_h = math.log(2) * fits("run6")["physical"]["tau_day"] * 24
    a(f"- Median error of fitted log10 k_top, runs 1 to 4: "
      f"{min(c['log10_ktop']['median_error_log10'] for c in c14.values()):+.3f} to "
      f"{max(c['log10_ktop']['median_error_log10'] for c in c14.values()):+.3f} decades. No material bias there.")
    a(f"- Run 6's first sample comes {t6_first_h:.1f} hours after the start of a {t6_irr_h:.1f} hour irradiation, "
      f"and its fitted release half-time is {t6_half_h:.1f} hours, so its sampling barely resolves the rate. In "
      f"{pct(cov['run6/fixedTBR']['log10_ktop']['fraction_error_over_1_decade'])} (TBR fixed) and "
      f"{pct(cov['run6/freeTBR']['log10_ktop']['fraction_error_over_1_decade'])} (TBR free) of synthetic "
      "replicates the fitted k_top was off by more than a decade. "
      + ("Treat run 6's k_top and alpha as soft." if run6_fragile else
         "With the full optimiser that failure did not show up at this replicate count. An earlier two-start fit "
         "without polish did show it, which points to the optimiser rather than to run 6."))
    a(f"- Flat direction stays flat: IV-only fits of synthetic run 1 data gave k_top no lower limit from the "
      f"data in {pct(flat['IV_only_not_identifiable_fraction'])} of {len(flat['reps'])} replicates. The same data "
      f"with OV gave a two-sided interval in {pct(flat['IV_plus_OV_identifiable_fraction'])}.")
    a(f"- Transfer test, both sides: under a true shared k_top the bootstrap 95% LR for the three He runs is "
      f"{tr['bootstrap_crit95']:.1f} (chi2 says {tr['chi2_crit95']:.1f}). With run 2's k_top truly 2 times higher, it rejects in "
      f"{pct(tr['power_if_second_run_differs_by']['x2.0'])} of replicates.")
    a(f"- Diffusivity, both sides (run 1 design, {dctl['well_mixed_truth']['reps']} replicates per truth): with "
      f"well-mixed truth the method returned a lower bound only in "
      f"{sum(s == 'lower bound only' for s in dctl['well_mixed_truth']['statuses'])} of "
      f"{dctl['well_mixed_truth']['reps']} and kept the truth inside 95% in "
      f"{pct(dctl['well_mixed_truth']['truth_inside_95_fraction'])}. With D = 7e-10 truth it kept the truth "
      f"in {pct(dcal7['truth_inside_95_fraction'])} and rejected well mixed in "
      f"{pct(dcal7['well_mixed_rejected_fraction'])}. With the Calderoni value at run 1's temperature as truth "
      f"({e(dcalc['D_true'], 1)} m2/s) it kept the truth in {pct(dcalc['truth_inside_95_fraction'])} and rejected "
      f"well mixed in {pct(dcalc['well_mixed_rejected_fraction'])}. The diffusivity fits keep their own single "
      "optimiser setting (L-BFGS-B without polish) at every D, the well-mixed end included.")
    a("")

    a("## Data and noise model")
    a("")
    rc = R["poisson_recount_check_run1"]
    a("- Runs 1 to 4 at the tags `case.yaml` (bundle root) uses (v0.6, v0.5, v0.2, v0.1). Run 6 at commit `82c5af9` (no tag "
      "exists). The raw LSC exports were fetched by `fetch_lsc.py` and pinned by sha256 in `lsc_manifest.json`. "
      "The three non-LSC files per run match `upstream/MANIFEST.json` (bundle root).")
    a("- `build_increments.py` runs each lab's own `tritium_model.py` unchanged against those files, with the "
      "toolbox wrapped so every vial carries its Poisson error (Bq / sqrt(CPMA x count time)). Shared blanks "
      "give correlated errors. The rebuilt cumulative curves equal `processed_data.json` to 1e-9 Bq in every run "
      "and stream.")
    a(f"- The Poisson formula checks out on run 1's six-fold recounts: the scatter of {rc['n_groups']} vials "
      f"counted 6 times is {rc['pooled_ratio']:.2f} times the predicted sigma ({rc['pooled_dof']} dof).")
    a("- Fits use per-sample increments, not the cumulative curve, so errors do not pile up along the curve.")
    a("- Scatter floor: counting noise plus a fitted additive term per stream (Bq). Two alternatives were "
      "rejected. A term proportional to the data let run 4's fit ignore its big early samples. A term proportional "
      "to the model stuck in local minima.")
    a("- Background: the main fits use the lab's numbers (one blank per file, negative vials set to zero by the "
      "toolbox). The run-6 background variant for runs 3 and 4 is described above. It changes the per-sample "
      "increments. The counting covariance is kept from the published background, which slightly understates the "
      "blank's weight in the corrected vials.")
    a("- The 0D model and geometry are the lab's: radius 7 cm, V = 1 L, wall 0.06 in, decay on. The solver matches "
      f"`release/r0.py` (bundle root) to {R['solver_vs_release_r0_max_rel_gap']:.0e}.")
    a("- Windows: run 3 is cut at its switch to H2 (day "
      f"{R['windows']['run3_He']['t_max_day']:.2f}) and run 4 at its switch to 3.5% H2 (day "
      f"{R['windows']['run4_1000ppm']['t_max_day']:.2f}). Run 6 is fitted whole, with k_top stepping at its switch "
      "to He.")
    a("- TBR range sensitivity: widening the allowed TBR range from a factor of 3.16 to 10 "
      + (f"leaves the IV + OV k_top, k_wall and alpha estimates and intervals unchanged (within 0.01 in log10) in "
         f"{names_list(wide_same)}" if wide_same else "changes every window")
      + (f", and moves them in {names_list(wide_moved)} (largest shift "
         f"{max(wide_shift[w] for w in wide_moved):.2g} in log10, or a status change)" if wide_moved else "")
      + ". See `freeTBR_wideTBR` in results.json. It also moves the prior-set lower edge of IV-only k_top.")
    a("- TBR as a nuisance: a scale on the lab's measured TBR, allowed within a factor of 3.16. Its 95% interval "
      f"includes 1 in {names_list(tbr_at_1)}"
      + (f" and excludes 1 in {names_list(tbr_not_1)}" if tbr_not_1 else "")
      + ". Landing at 1 is expected, because the measured TBR is computed from these same collected totals. "
      "It is not independent evidence.")
    geo = R["geometry"]  # side-wall area grows with height; the bottom (= top area) does not
    a_wall_best = geo["A_top_m2"] + (geo["A_wall_m2"] - geo["A_top_m2"]) * si["height_ratio_best_over_model"]
    kw_factor = si["volume_ratio_logged_over_model"] / (a_wall_best / geo["A_wall_m2"])
    a(f"- Salt volume caveat. The models use the lab's 1 L column ({100 * si['height_model_m']:.2f} cm at 7 cm "
      f"radius). The logged charge is {si['salt_mass_logged_kg']:.2f} kg, which at the lab's density of "
      f"{si['density_g_cm3']:.4f} g/cm3 ({si['density_temperature_C']} C) is {1000 * si['volume_logged_m3']:.3f} L, "
      f"a column about {100 * si['height_best_estimate_m']:.2f} cm tall (`case.yaml` best estimate, bundle root). The fits "
      "measure release rates per unit volume (alpha), so with the logged charge every absolute k_top here would "
      f"be {100 * (si['volume_ratio_logged_over_model'] - 1):.0f}% higher (volume over top area, which does "
      f"not change) and k_wall about {100 * (kw_factor - 1):.0f}% higher (the wall area grows with the column). "
      f"Diffusivity thresholds scale with height squared, a factor of {si['diffusivity_scale_height_squared']:.2f}. "
      "Ratios between runs, alpha and tau, and every identifiability status are unaffected.")
    a("- `requirements.txt` pins `libra-toolbox` at `4c7ef88`: the repo's pinned `cb41e93` plus the merged PR #94, "
      "which only lets the package import without OpenMC. Nothing under `libra_toolbox/tritium/` differs.")
    a("")

    a("## What this does not show")
    a("")
    a("- It does not say which physics sets k_top. The 0D model lumps salt mixing, the gas-side film and "
      "surface kinetics into one number.")
    a("- The intervals are conditional on the model. "
      + ("The OV test shows the model is wrong for OV timing, so k_wall intervals mean little, "
         if lagf_bad else
         "The OV test cannot confirm the model's OV timing once scatter floors are included, so k_wall intervals "
         "mean little, ")
      + "and k_top intervals are best read as \"what this model needs\".")
    a("- Runs 7 to 10 and the FLiBe-1L runs are not used.")
    a("")

    a("## Reproduce")
    a("")
    pr, pd = R["provenance"], DF["provenance"]

    def argstr(d):
        return " ".join(f"--{k.replace('_', '-')} {v}" for k, v in d.items() if k != "seed")

    def when(pv, rt):
        return (f"# {pv['timestamp_utc'][:16]}, {rt / 60:.1f} min, {pv['cpu_count']} cores, "
                f"code {(pv['git_sha'] or 'unknown')[:7]}{' + local edits' if pv['code_dirty'] else ''}")
    py = "python studies/release_identifiability/"
    a("From this folder, with Python 3.12 or later:")
    a("")
    a("```bash")
    a("make venv             # .venv from requirements.txt (libra-toolbox pinned to 4c7ef88)")
    a("make data             # re-download the raw LSC files (sha256-pinned), rebuild increments.json and quench.json, check both are identical")
    a("make reproduce-fast   # every stage at low replicate counts, APPENDIX.md and README.md stamped DRAFT")
    a("make reproduce-full   # the replicate counts the claims need, hours of CPU, use a many-core machine")
    a("```")
    a("")
    a("Measured times are in the `Makefile`. Both reproduce targets overwrite the committed results. "
      "Re-extracting the bundle restores them. The FESTIM cross-check needs the repository's FESTIM model, which "
      "is not in the bundle. Without it the check is skipped with a message and the committed "
      "`festim_check.json` is used.")
    a("")
    a("The commands and replicate counts that produced the committed results, read from the provenance "
      "records:")
    a("")
    a("```bash")
    a("python fetch_lsc.py  # raw LSC files, sha256-pinned")
    a("python build_increments.py")
    a("python quench.py")
    a(f"python analysis.py --part fits {argstr(pr['fits']['arguments'])}  {when(pr['fits'], R['runtime_fits_s'])}")
    a(f"python analysis.py --part controls {argstr(pr['controls']['arguments'])}  "
      f"{when(pr['controls'], R['runtime_controls_s'])}")
    a(f"python analysis.py --part transfer {argstr(pr['transfer']['arguments'])}  "
      f"{when(pr['transfer'], R['runtime_transfer_s'])}")
    a(f"python change.py --reps {CH['provenance']['reps']}  # {CH['provenance']['seconds'] / 60:.1f} min")
    a("python festim_check.py  # optional, needs FESTIM and the repository's FESTIM model")
    a(f"python diffusivity.py --part windows {argstr(pd['windows']['arguments'])}  "
      f"{when(pd['windows'], DF['runtime_windows_s'])}")
    a(f"python diffusivity.py --part controls {argstr(pd['controls']['arguments'])}  "
      f"{when(pd['controls'], DF['runtime_controls_s'])}")
    a("python make_figures.py")
    a("python make_readme.py")
    a("```")
    a("")
    a("`make_readme.py` refuses to build if any part is missing its provenance, was computed from a different "
      "`increments.json` or `quench.json`, if controls or transfer were computed from different fits than the ones "
      "stored, or if replicate counts are below what the claims need (unless `--draft`).")
    a("")
    a("## Files")
    a("")
    for fn, what in (("fetch_lsc.py", "download raw lab files at pinned refs"),
                     ("lsc_manifest.json", "sha256 pins for those files"),
                     ("build_increments.py", "run the lab's scripts, attach counting covariance"),
                     ("increments.json", "per-sample increments and covariance, runs 1-4 and 6"),
                     ("quench.py / quench.json / quench_variants.json", "the run-6 tSIE background applied to runs 3 and 4, curve-shape variants"),
                     ("model0d.py", "0D model with an optional gas step"),
                     ("fitlib.py", "likelihood, fits, profiles, Fisher"),
                     ("analysis.py", "fits, identifiability, transferability, synthetic controls"),
                     ("results.json", "all 0D numbers"),
                     ("change.py / change.json", "run 2 with a change in k_top during the run"),
                     ("fv1d.py", "1D finite-volume salt column"),
                     ("festim_check.py / festim_check.json", "FESTIM cross-check of fv1d"),
                     ("diffusivity.py / diffusivity.json", "profile of D and its controls"),
                     ("make_figures.py", "the two figures"),
                     ("make_readme.py", "writes README.md and this file"),
                     ("requirements.txt / Makefile", "pinned environment and the reproduce targets")):
        a(f"- `{fn}`: {what}")
    a("")
    (HERE / "APPENDIX.md").write_text("\n".join(L))
    print("wrote README.md and APPENDIX.md")


if __name__ == "__main__":
    main()
