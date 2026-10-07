<!-- GENERATED FILE. Built by make_readme.py from the JSON files in this folder. Do not edit by hand. -->

# BABY-1L: the run-6 tSIE background, applied to runs 3 and 4

The per-vial data are rebuilt from the public raw LSC files with each run's own `tritium_model.py` and match `processed_data.json` exactly. Every number below is computed by the scripts in this folder. Fits, controls and caveats are in APPENDIX.md.

Up to run 4 each vial has its file's blank subtracted, counted at tSIE 280 to 292. In runs 3 and 4, 41 vials sit above tSIE 330, and 40 of them are the first bubbler vial of a sample, whose water evaporates between samples. The run-6 blank scan (`SLR_BLANK_SCAN.csv`) puts the background at those tSIE values well below the blank: 0.31 Bq at tSIE 278, 0.24 at 324, 0.17 at 489. From run 6 on, `tritium_model.py` subtracts `background_curve(tSIE)`. Applied to runs 3 and 4, with the curve scaled so it equals each file's own blank at the blank's tSIE:

| Vials corrected | Curve below tSIE 300 | Run 3 IV + OV (Bq) | Run 3 measured TBR | Run 3 OV soluble | Run 4 IV + OV (Bq) | Run 4 measured TBR | Run 4 OV soluble |
| --- | --- | --- | --- | --- | --- | --- | --- |
| none (as published) | | 29.65 ± 0.37 | 4.03e-03 | 17% | 54.90 ± 0.44 | 2.58e-03 | 1% |
| above tSIE 330 only | held flat below tSIE 300 | 30.45 ± 0.45 | 4.13e-03 (+2.7%) | 26% | 55.27 ± 0.49 | 2.59e-03 (+0.7%) | 3% |
| above tSIE 330 only | straight-line fit | 30.45 ± 0.41 | 4.13e-03 (+2.7%) | 27% | 55.10 ± 0.46 | 2.59e-03 (+0.4%) | 2% |
| above tSIE 330 only | quadratic fit | 30.92 ± 0.46 | 4.20e-03 (+4.3%) | 30% | 55.45 ± 0.47 | 2.60e-03 (+1.0%) | 4% |
| above tSIE 330 only | as interpolated (run 6) | 30.98 ± 0.47 | 4.21e-03 (+4.5%) | 30% | 55.41 ± 0.47 | 2.60e-03 (+0.9%) | 4% |
| all | held flat below tSIE 300 | 30.50 ± 0.48 | 4.14e-03 (+2.9%) | 26% | 55.36 ± 0.50 | 2.60e-03 (+0.8%) | 3% |
| all | straight-line fit | 31.01 ± 0.44 | 4.21e-03 (+4.6%) | 28% | 55.36 ± 0.47 | 2.60e-03 (+0.8%) | 3% |
| all | quadratic fit | 32.17 ± 0.59 | 4.37e-03 (+8.5%) | 32% | 56.07 ± 0.53 | 2.63e-03 (+2.1%) | 5% |
| all | as interpolated (run 6) | 33.10 ± 1.13 | 4.49e-03 (+11.6%) | 32% | 56.32 ± 0.61 | 2.64e-03 (+2.6%) | 6% |

Correcting only the vials above tSIE 330 raises the measured TBR by 2.7% to 4.5% in run 3 and 0.4% to 1.0% in run 4 across the four curve shapes. That is the part that holds up. Correcting every vial adds more, up to 11.6% and 2.6%, but that depends on the curve between tSIE 278 and 294, where the scan's first two points drop by 0.04 Bq (about 2 sigma) and the ordinary vials sit just above their blanks. Read the all-vials rows as the upper end of a range. ± is counting noise (vials, blanks, scan points, 2000 draws) and does not include the choice of curve shape. With the published background, low-activity vials above tSIE 330 average -0.046 Bq (run 3) and -0.030 Bq (run 4), below zero.

Runs 1 and 2 show the same sign (low-activity vials above tSIE 330 average -0.034 and -0.027 Bq). They were counted with the 3H-UG-LLT-NM quench set, not the 3H-UG-LLT-NM-C set of the scan, so the run-6 curve cannot be carried over. A blank scan with that quench set would settle them.

## A note on run 2

The run-2 script sets k_top to 0.7 times run 1's value and k_wall to 0 (`analysis/tritium/tritium_model.py`, lines 223 to 224 at v0.5). Fitted with one k_top at the lab's TBR, run 2's k_top is 0.48 of run 1's fitted value and 0.55 of run 1's script value. On the total release rate (top plus wall), the scripts imply 0.67 and the fits give 0.47. But one k_top does not describe run 2: its IV scatter floor is 1.15 Bq, against 0.14 to 0.48 in the other He runs. Letting k_top change once fits much better with the change at the day-13 sample (2ΔNLL 46, p = 0.005, the smallest possible with 200 simulations. The 5 days tried were picked after looking at the data, and the best is calibrated by simulation). The release time goes from 15.8 to 11.1 days, and the later k_top is 0.75 of run 1's script value (95%: 0.72 to 0.79), just above the 0.7 in the script. The paper suggests the salt freezing mid-run as a possible explanation for run 2's slower release (sec. 3.3). A partly frozen salt could look like that. A full stop in release fits worse than one constant k_top, and the run's files do not say when the salt froze, so this is not tested. Details: APPENDIX.md, Run 2.

## What this does not show

- The 0D model lumps salt mixing, the gas-side film and surface kinetics into k_top. The intervals are conditional on it, and synthetic tests put the nominal 95% intervals for k_top at 91% (APPENDIX.md, Controls).
- Refitting runs 3 and 4 with every vial corrected (the upper end) moves k_top by 2% and 0% (lab TBR). The background changes the collected totals and chemical form, not the release rates.
- These curves do not measure the salt diffusivity, and nothing here tests the paper's diffusion-limited reading (APPENDIX.md, Salt diffusivity).
- Runs 7 to 10 and the FLiBe-1L runs are not used. Paper: arXiv 2509.26174.
