#!/usr/bin/env python3
"""Figures from results.json and diffusivity.json (no numbers typed here).

  fig_profiles.png   profile likelihoods: k_top, k_wall, salt diffusivity D
  fig_k_per_run.png  fitted k_top and k_wall per run window, 95% profile CIs

Usage: run/refit/.venv/bin/python studies/release_identifiability/make_figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
R = json.loads((HERE / "results.json").read_text())
DF = json.loads((HERE / "diffusivity.json").read_text())

# validated categorical order (dataviz reference palette, light mode), fixed per window
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
MARKERS = ["o", "s", "^", "D", "v"]
WINDOWS = ["run1", "run2", "run3_He", "run4_1000ppm", "run6"]
LABEL = {"run1": "Run 1 (He)", "run2": "Run 2 (He)", "run3_He": "Run 3 He part",
         "run4_1000ppm": "Run 4 (1000 ppm H2)", "run6": "Run 6 (3.5% H2)"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CHI2_95 = 3.841458820694124


def style(ax):
    ax.set_facecolor("#fcfcfb")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.axhline(CHI2_95, color=MUTED, lw=1, ls=(0, (4, 3)))


def fig_profiles():
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), facecolor="#fcfcfb")
    for ax, pname, title in ((axs[0], "log10_ktop", "k_top (m/s)"), (axs[1], "log10_kwall", "k_wall (m/s)")):
        style(ax)
        for i, w in enumerate(WINDOWS):
            p = R["fits"][w]["freeTBR"]["profiles"][pname]
            ax.plot(p["grid"], np.minimum(p["two_delta_nll"], 30), color=COLORS[i], lw=2,
                    marker=MARKERS[i], ms=4, markevery=4, label=LABEL[w])
        if pname == "log10_ktop":
            p = R["fits"]["run1"]["IVonly_freeTBR"]["profiles"][pname]
            ax.plot(p["grid"], np.minimum(p["two_delta_nll"], 30), color=COLORS[0], lw=2, ls=":",
                    label="Run 1, IV only (flat)")
        ax.set_ylim(0, 30)
        ax.set_xlabel(f"log10 {title}", color=INK, fontsize=9)
        ax.set_title(f"Profile of {title}, TBR free", color=INK, fontsize=10, loc="left")
    axs[0].set_ylabel("2 x (NLL - min)   (dashed: 95%)", color=INK, fontsize=9)
    axs[0].legend(fontsize=7, frameon=False, loc="upper left")
    ax = axs[2]
    style(ax)
    g = DF["grid_log10_D"]
    for i, w in enumerate(WINDOWS):
        r = DF["runs"]["free"][w]
        ax.plot(g, np.minimum(r["two_delta_nll"], 30), color=COLORS[i], lw=2, marker=MARKERS[i], ms=4,
                markevery=4, label=LABEL[w])
    dcal = DF["runs"]["free"]["run1"]["D_calderoni_m2_s"]
    for x, t in ((np.log10(dcal), "Calderoni T"), (np.log10(DF["D_H2_ref_m2_s"]), "H2 7e-10")):
        ax.axvline(x, color=MUTED, lw=1)
        ax.text(x + 0.05, 27, t, color=MUTED, fontsize=7, rotation=90, va="top")
    ax.set_ylim(0, 30)
    ax.set_xlabel("log10 salt diffusivity D (m2/s); right edge = well mixed", color=INK, fontsize=9)
    ax.set_title("Profile of D, 1D column, TBR free", color=INK, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(HERE / "fig_profiles.png", dpi=150, facecolor="#fcfcfb")


def ci(block, name):
    p = block["profiles"][name]
    c = p["ci95"]
    lo = c["lo"] if c["lo"] is not None and c["status"] in ("identifiable", "lower bound only") else None
    hi = c["hi"] if c["hi"] is not None and c["status"] in ("identifiable", "upper bound only") else None
    return p["mle"], lo, hi


def fig_k():
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), facecolor="#fcfcfb")
    for ax, name, title, labkey in ((axs[0], "log10_ktop", "k_top (m/s)", "lab_k_top_m_s"),
                                    (axs[1], "log10_kwall", "k_wall (m/s)", "lab_k_wall_m_s")):
        ax.set_facecolor("#fcfcfb")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(True, axis="y", color=GRID, lw=0.6)
        for i, w in enumerate(WINDOWS):
            for dx, var, filled in ((-0.15, "fixedTBR", False), (0.15, "freeTBR", True)):
                m, lo, hi = ci(R["fits"][w][var], name)
                lo_v = lo if lo is not None else m - 1.5
                hi_v = hi if hi is not None else m + 1.5
                ax.plot([i + dx, i + dx], [10 ** lo_v, 10 ** hi_v], color=COLORS[i], lw=2,
                        ls="-" if (lo is not None and hi is not None) else ":")
                ax.plot(i + dx, 10 ** m, marker=MARKERS[i], ms=8, color=COLORS[i],
                        mfc=COLORS[i] if filled else "#fcfcfb", mec=COLORS[i], mew=2)
            lab = R["windows"][w][labkey]
            if lab and lab > 0:
                ax.plot(i, lab, marker="x", ms=8, color=INK, mew=1.5)
        ax.set_yscale("log")
        ax.set_xticks(range(len(WINDOWS)))
        ax.set_xticklabels([LABEL[w] for w in WINDOWS], fontsize=8, rotation=15, color=INK)
        ax.set_title(f"{title}: open = TBR fixed, filled = TBR free, x = lab value", color=INK,
                     fontsize=10, loc="left")
    axs[1].text(0.01, 0.02, "dotted = open side (bound only)", transform=axs[1].transAxes, color=MUTED,
                fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "fig_k_per_run.png", dpi=150, facecolor="#fcfcfb")


if __name__ == "__main__":
    fig_profiles()
    fig_k()
    print("wrote fig_profiles.png, fig_k_per_run.png")
