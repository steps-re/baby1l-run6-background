#!/usr/bin/env python3
"""Cross-check the 1D finite-volume solver (fv1d.py) against FESTIM.

Uses the baby-benchmark repository's FESTIM rung, run/festim/r1_robin.py
(zero-flux bottom, Robin top; not part of the share bundle, so in the bundle
this check is skipped and the committed festim_check.json is used), at two finite diffusivities that bracket the study range:
the H2-in-FLiBe value near 900 K (7e-10 m2/s) and the Calderoni T value at
the run 1 salt temperature. Run 1 irradiation schedule, k_top = 1e-7 m/s,
k_wall = 0 (FESTIM's rung has no bottom flux).

Run in the FESTIM env:
  UCX_TLS=tcp,self,sm conda run -n festim-dev python studies/release_identifiability/festim_check.py
Writes festim_check.json next to this file.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "run" / "festim"))

import fv1d  # noqa: E402
from fv1d import d_calderoni_T  # noqa: E402
import model0d as M0  # noqa: E402
try:
    from r1_robin import FESTIM_COMMIT, run_robin_simulation  # noqa: E402
except ImportError as exc:  # FESTIM (dolfinx) not installed
    print(f"festim_check skipped: FESTIM is not importable here ({exc}). This step is an optional "
          "cross-check of fv1d.py; the committed festim_check.json is used as is. Run it in a FESTIM "
          "environment to regenerate it.")
    sys.exit(0)
from release.analytical import Window  # noqa: E402

def main() -> int:
    inc = json.loads((HERE / "increments.json").read_text())
    r1 = inc["runs"]["1"]
    T_K = r1["temperature_salt_C"] + 273.15
    irr = [tuple(x) for x in r1["irradiations_s"]]
    src = r1["measured_TBR"] * r1["neutron_rate_n_s"]
    t = np.array(r1["streams"]["IV"]["times_day"]) * 86400.0
    k_top = 1.0e-7
    out = {"festim_commit_of_rung": FESTIM_COMMIT, "k_top_m_s": k_top, "T_K": T_K, "cases": []}
    for D in (7.0e-10, d_calderoni_T(T_K)):
        windows = [Window(a, b, src / M0.V_M3) for a, b in irr]
        t0 = time.time()
        fe = run_robin_simulation(M0.H_M, D, k_top, windows, t, n_cells=50)
        secs = time.time() - t0
        q_fe = np.asarray(fe["cumulative_release_at_samples"]) * M0.A_TOP
        q_fv, _ = fv1d.cumulative_release(D, [k_top], 0.0, src, irr, None, t, [], n=60, decay=False)
        q_0d, _ = M0.cumulative_release([k_top], 0.0, src, irr, None, t, [], decay=False)
        out["cases"].append({
            "D_m2_s": D, "festim_seconds": secs, "festim_steps": fe["n_steps"],
            "festim_cum_top": q_fe.tolist(), "fv_cum_top": q_fv.tolist(), "zeroD_cum_top": q_0d.tolist(),
            "max_rel_gap_fv_vs_festim": float(np.max(np.abs(q_fv / q_fe - 1))),
            "max_rel_gap_0d_vs_festim": float(np.max(np.abs(q_0d / q_fe - 1))),
        })
        print(f"D={D:.3g}: FESTIM {secs:.0f}s, FV vs FESTIM max gap "
              f"{out['cases'][-1]['max_rel_gap_fv_vs_festim']:.2e}, 0D vs FESTIM "
              f"{out['cases'][-1]['max_rel_gap_0d_vs_festim']:.2e}", flush=True)
    (HERE / "festim_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
