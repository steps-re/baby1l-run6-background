#!/usr/bin/env python3
"""Rebuild each run's release increments, with a counting covariance, by running
the lab's own tritium_model.py against the raw LSC files.

The lab script is executed unchanged. Three toolbox methods are wrapped so that
every vial activity also carries a linear combination of raw LSC rows:

  LSCFileReader.read_file   per row, gross Poisson sigma = Bq / sqrt(CPMA * T)
  LSCSample.from_file       the vial starts as +1 * (its row)
  LSCSample.substract_background   subtracts the blank's row combination

A background the lab subtracted as a bare number (run 1 samples 1-2 use a typed
0.320 Bq; run 6 interpolates blanks against tSIE) has no row. It gets one shared
pseudo-row per run whose variance is the median blank variance of that run.

The lab script stops with FileNotFoundError at the OpenMC statepoint, after the
gas streams are built. That is where this script reads them.

Covariance of the per-sample increments is A diag(var) A^T over all raw rows,
so blanks shared between samples and streams are correlated correctly.

Check: the cumulative sums rebuilt here must equal processed_data.json at the
same pinned ref to 1e-9 Bq, or the run is rejected.

A second check tests the Poisson sigma itself: run 1 recounted five OV vials
and one blank six times each (OV-1-1_with_avg.csv). The scatter of the six
recounts is compared with the predicted counting sigma.

Run with the refit venv (libra-toolbox at cb41e93 + one import-only commit):
  run/refit/.venv/bin/python studies/release_identifiability/build_increments.py
Writes studies/release_identifiability/increments.json.
"""

from __future__ import annotations

import json
import math
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fetch_lsc import RUNS, folder  # noqa: E402

import libra_toolbox.tritium.lsc_measurements as lsc  # noqa: E402
from libra_toolbox.tritium.model import ureg  # noqa: E402

OUT = HERE / "increments.json"

ROW_VAR: dict = {}  # row id -> variance (Bq^2); NaN = fill later
ROW_KIND: dict = {}  # row id -> "vial" | "blank" | "scalar_bg"
CURRENT_RUN = [None]

_orig_read = lsc.LSCFileReader.read_file
_orig_from_file = lsc.LSCSample.from_file
_orig_sub = lsc.LSCSample.substract_background
_orig_init = lsc.LSCSample.__init__


def _row_sigmas(df: pd.DataFrame) -> list[float]:
    """Gross counting sigma in Bq per data row, from CPMA and count time.

    Rows with no count time are the counter's average of the rows above them
    (run 1's six-fold recount). Their effective count time is the sum.
    """
    sig = []
    acc_t = 0.0
    for _, row in df.iterrows():
        bq = pd.to_numeric(row.get("Bq:1"), errors="coerce")
        cpm = pd.to_numeric(row.get("CPMA"), errors="coerce")
        t = pd.to_numeric(row.get("Count Time"), errors="coerce")
        if pd.isna(t):
            t_eff = acc_t
            acc_t = 0.0
        else:
            t_eff = float(t)
            acc_t += float(t)
        if pd.isna(bq) or pd.isna(cpm) or cpm <= 0 or t_eff <= 0:
            sig.append(float("nan"))
            if pd.isna(bq):
                acc_t = 0.0
        else:
            sig.append(float(bq) / math.sqrt(float(cpm) * t_eff))
    return sig


def read_file(self):
    _orig_read(self)
    self._sigmas = _row_sigmas(self.data)


def _label_map(reader) -> dict:
    """Same zip semantics as get_bq1_values_with_labels: later label wins."""
    out = {}
    for i, lab in enumerate(reader.vial_labels):
        out[lab] = i
    return out


def from_file(reader, vial_name):
    s = _orig_from_file(reader, vial_name)
    i = _label_map(reader)[vial_name]
    rid = (CURRENT_RUN[0], os.path.basename(reader.file_path), i, vial_name)
    ROW_VAR[rid] = reader._sigmas[i] ** 2
    ROW_KIND[rid] = "blank" if ("BL" in str(vial_name) and "IV" not in str(vial_name)
                                 and "OV" not in str(vial_name)) else "vial"
    s.terms = {rid: 1.0}
    s.gross = float(s.activity.magnitude)
    s._reader, s._label = reader, vial_name
    return s


def init(self, activity, name):
    _orig_init(self, activity, name)
    self.terms = None
    self.gross = float(activity.magnitude)
    self.clamped = False


def substract_background(self, background_sample):
    before = float(self.activity.magnitude)
    bg = float(background_sample.activity.magnitude)
    _orig_sub(self, background_sample)
    bterms = getattr(background_sample, "terms", None)
    if not bterms:
        rid = (CURRENT_RUN[0], "scalar_bg", -1, "scalar_bg")
        ROW_VAR.setdefault(rid, float("nan"))
        ROW_KIND[rid] = "scalar_bg"
        bterms = {rid: 1.0}
    terms = dict(self.terms)
    for k, v in bterms.items():
        terms[k] = terms.get(k, 0.0) - v
    self.terms = terms
    self.clamped = before - bg < 0
    self.clamp_deficit = max(0.0, bg - before)  # Bq added by the lab's clamp to zero


lsc.LSCFileReader.read_file = read_file
lsc.LSCSample.from_file = staticmethod(from_file)
lsc.LSCSample.__init__ = init
lsc.LSCSample.substract_background = substract_background


def exec_lab_script(run: int) -> dict:
    tdir = folder(run) / "analysis" / "tritium"
    src = (tdir / "tritium_model.py").read_text()
    ns: dict = {"__name__": "lab_tritium_model"}
    cwd = os.getcwd()
    os.chdir(tdir)
    CURRENT_RUN[0] = run
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exec(compile(src, str(tdir / "tritium_model.py"), "exec"), ns)
        raise RuntimeError(f"run {run}: lab script ran past the statepoint; expected it to stop there")
    except FileNotFoundError as e:
        if "statepoint" not in str(e):
            raise
    finally:
        os.chdir(cwd)
    return ns


def fix_scalar_bg(run: int, ns: dict):
    """Samples whose background the lab subtracted as a number outside the
    toolbox (run 6's tSIE curve) have no background term yet. Give them the
    run's shared scalar-background pseudo-row."""
    rid = (run, "scalar_bg", -1, "scalar_bg")
    streams = ns.get("gas_streams") or {"IV": ns["IV_stream"], "OV": ns["OV_stream"]}
    for st in streams.values():
        for smp in st.samples:
            for v in smp.samples:
                if len(v.terms) == 1 and v.background_substracted:
                    ROW_VAR.setdefault(rid, float("nan"))
                    ROW_KIND[rid] = "scalar_bg"
                    v.terms = {**v.terms, rid: -1.0}
                    v.clamped = float(v.activity.magnitude) == 0.0 and v.gross > 0
                    v.clamp_deficit = 0.0
                    if v.clamped:
                        # the lab subtracted background_curve(tSIE), see its create_sample
                        row = ns["get_row_by_label"](v._reader, v._label)
                        bg = float(ns["background_curve"](float(row["tSIE"])))
                        v.clamp_deficit = max(0.0, bg - v.gross)
    blanks = [ROW_VAR[k] for k in ROW_VAR if k[0] == run and ROW_KIND[k] == "blank"
              and not math.isnan(ROW_VAR[k])]
    if not blanks:  # run 6 blanks come from the blank-set file, read via from_file too
        blanks = [ROW_VAR[k] for k in ROW_VAR if k[0] == run and "BL" in str(k[3]).upper()
                  and not math.isnan(ROW_VAR[k])]
    med = float(np.median(blanks))
    for k in ROW_VAR:
        if k[0] == run and ROW_KIND[k] == "scalar_bg":
            ROW_VAR[k] = med
    return streams, med, len(blanks)


def seconds_rel(general: dict, start):
    from datetime import datetime
    sw = general.get("cover_gas", {}).get("switched_to")
    if isinstance(sw, dict) and sw.get("gas_switch_time"):
        t = datetime.strptime(sw["gas_switch_time"], "%m/%d/%Y %H:%M")
        return (t - start).total_seconds(), sw["type"]
    return None, None


def build_run(run: int) -> dict:
    ns = exec_lab_script(run)
    streams, med_blank_var, n_blanks = fix_scalar_bg(run, ns)
    general = json.loads((folder(run) / "data" / "general.json").read_text())
    processed = json.loads((folder(run) / "data" / "processed_data.json").read_text())
    start = ns["start_time"]

    # the lab builds irradiations and cumulative curves in list order, make sure
    # samples are time-ordered so increments are well defined
    rows, meta = [], []
    for sname in ("IV", "OV"):
        st = streams[sname]
        t = [x.total_seconds() for x in st.relative_times]
        if any(b < a for a, b in zip(t, t[1:])):
            raise RuntimeError(f"run {run} {sname}: samples not time ordered")
        for j, smp in enumerate(st.samples):
            for form, vials in (("soluble", smp.samples[:2]), ("insoluble", smp.samples[2:])):
                val = sum(float(v.activity.magnitude) for v in vials)
                terms: dict = {}
                for v in vials:
                    for k, c in v.terms.items():
                        terms[k] = terms.get(k, 0.0) + c
                rows.append((sname, form, j, t[j], val, terms,
                             sum(bool(getattr(v, "clamped", False)) for v in vials),
                             sum(getattr(v, "clamp_deficit", 0.0) for v in vials)))
    keys = sorted({k for r in rows for k in r[5]}, key=str)
    kidx = {k: i for i, k in enumerate(keys)}
    var = np.array([ROW_VAR[k] for k in keys])
    if np.any(~np.isfinite(var)):
        raise RuntimeError(f"run {run}: rows without a counting variance: "
                           f"{[k for k in keys if not np.isfinite(ROW_VAR[k])]}")
    A = np.zeros((len(rows), len(keys)))
    for i, r in enumerate(rows):
        for k, c in r[5].items():
            A[i, kidx[k]] += c
    cov_form = A @ np.diag(var) @ A.T

    # totals = soluble + insoluble per (stream, sample)
    out_streams = {}
    tot_rows, tot_A = [], []
    for sname in ("IV", "OV"):
        idx_s = [i for i, r in enumerate(rows) if r[0] == sname and r[1] == "soluble"]
        idx_i = [i for i, r in enumerate(rows) if r[0] == sname and r[1] == "insoluble"]
        times = [rows[i][3] / 86400.0 for i in idx_s]
        sol = [rows[i][4] for i in idx_s]
        ins = [rows[i][4] for i in idx_i]
        tot = [a + b for a, b in zip(sol, ins)]
        for a, b in zip(idx_s, idx_i):
            tot_A.append(A[a] + A[b])
            tot_rows.append((sname, rows[a][2]))
        cum_lab = processed["cumulative_tritium_release"][sname]
        diffs = {
            "total": float(np.max(np.abs(np.cumsum(tot) - np.array(cum_lab["total"]["value"])))),
            "soluble": float(np.max(np.abs(np.cumsum(sol) - np.array(cum_lab["soluble"]["value"])))),
            "insoluble": float(np.max(np.abs(np.cumsum(ins) - np.array(cum_lab["insoluble"]["value"])))),
            "times_day": float(np.max(np.abs(np.array(times) - np.array(cum_lab["sampling_times"]["value"])))),
        }
        if max(diffs.values()) > 1e-9:
            raise RuntimeError(f"run {run} {sname}: rebuilt curve differs from processed_data.json: {diffs}")
        out_streams[sname] = {
            "times_day": times,
            "inc_total_Bq": tot, "inc_soluble_Bq": sol, "inc_insoluble_Bq": ins,
            "clamped_vials": [rows[a][6] + rows[b][6] for a, b in zip(idx_s, idx_i)],
            "clamp_added_Bq": [rows[a][7] + rows[b][7] for a, b in zip(idx_s, idx_i)],
            "clamp_added_soluble_Bq": [rows[a][7] for a in idx_s],
            "clamp_added_insoluble_Bq": [rows[b][7] for b in idx_i],
            "max_abs_diff_vs_processed_data": diffs,
        }
    TA = np.array(tot_A)
    cov_tot = TA @ np.diag(var) @ TA.T
    sw_s, sw_type = seconds_rel(general, start)
    return {
        "run": run,
        "ref": f"{RUNS[run]['repo']}@{RUNS[run]['ref']} ({RUNS[run]['commit'][:12]})",
        "gas_initial": general["cover_gas"]["type"],
        "gas_switch_s": sw_s,
        "gas_switched_to": sw_type,
        "temperature_salt_C": general["general_data"].get("temperature_salt", {}).get("value")
        if isinstance(general["general_data"].get("temperature_salt"), dict) else None,
        "irradiations_s": [[i["start_time"]["value"], i["stop_time"]["value"]] for i in processed["irradiations"]],
        "neutron_rate_n_s": processed["neutron_rate_used_in_model"]["value"],
        "neutron_rate_err_n_s": processed["neutron_rate_used_in_model"].get("error"),
        "measured_TBR": processed["measured_TBR"]["value"],
        "TBR_used_in_model": processed["TBR_used_in_model"]["value"],
        "lab_k_top_m_s": processed["k_top"]["value"],
        "lab_k_wall_m_s": processed["k_wall"]["value"],
        "modelled_radius_cm": processed["modelled_baby_radius"]["value"],
        "streams": out_streams,
        "total_vector_order": [f"{s}[{j}]" for s, j in tot_rows],
        "cov_total_Bq2": cov_tot.tolist(),
        "n_raw_rows": len(keys),
        "scalar_background_variance_Bq2": med_blank_var,
        "n_blanks_for_scalar_background": n_blanks,
    }


def recount_check() -> dict:
    """Run 1 recounted OV vials six times: compare scatter with Poisson sigma."""
    p = folder(1) / "data" / "tritium_detection" / "OV-1-1_with_avg.csv"
    r = lsc.LSCFileReader(str(p), vial_labels=[None] * 1)  # labels unused here
    _orig_read(r)
    df = r.data
    sig = _row_sigmas(df)
    groups, cur = [], []
    for i, row in df.iterrows():
        t = pd.to_numeric(row.get("Count Time"), errors="coerce")
        bq = pd.to_numeric(row.get("Bq:1"), errors="coerce")
        if pd.isna(bq):
            cur = []
            continue
        if pd.isna(t):
            if cur:
                groups.append(cur)
            cur = []
            continue
        cur.append((float(bq), sig[i]))
    ratios = []
    for g in groups:
        vals = np.array([a for a, _ in g])
        pred = float(np.mean([b for _, b in g]))
        ratios.append(float(np.std(vals, ddof=1) / pred))
    # pooled: sum of squared deviations over pooled predicted variance
    ss = sum(float(np.sum((np.array([a for a, _ in g]) - np.mean([a for a, _ in g])) ** 2)) for g in groups)
    pv = sum((len(g) - 1) * float(np.mean([b for _, b in g])) ** 2 for g in groups)
    dof = sum(len(g) - 1 for g in groups)
    return {"n_groups": len(groups), "n_recounts_each": [len(g) for g in groups],
            "sd_over_predicted_sigma": ratios,
            "pooled_ratio": math.sqrt(ss / pv), "pooled_dof": dof}


def main() -> int:
    from libra_toolbox.tritium.model import quantity_to_activity
    bq_per_particle = float(quantity_to_activity(1 * ureg.particle).to(ureg.Bq).magnitude)
    out = {"_generated": "Written by build_increments.py from the lab's raw LSC files. Do not edit.",
           "Bq_per_particle": bq_per_particle,
           "Bq_per_particle_origin": "libra_toolbox.tritium.model.quantity_to_activity (the lab's conversion)",
           "runs": {}}
    for run in RUNS:
        out["runs"][str(run)] = build_run(run)
        r = out["runs"][str(run)]
        print(f"run {run}: IV {len(r['streams']['IV']['times_day'])} pts, "
              f"OV {len(r['streams']['OV']['times_day'])} pts, {r['n_raw_rows']} raw rows, "
              f"rebuilt = processed_data.json")
    out["poisson_recount_check_run1"] = recount_check()
    print("recount check:", out["poisson_recount_check_run1"])
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
