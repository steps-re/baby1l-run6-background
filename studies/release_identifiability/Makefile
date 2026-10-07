# One-command reproduction of the release-identifiability study.
#   make venv             create .venv here from requirements.txt
#   make data             re-download the raw LSC files (sha256-pinned) and rebuild
#                         increments.json, quench.json and quench_variants.json with the lab's
#                         scripts; fails if either differs from the committed one
#   make reproduce-fast   every stage at low replicate counts.
#                         Overwrites the results, figures, README.md and APPENDIX.md
#                         with a build stamped DRAFT. Measured 2026-10-07: venv + data +
#                         reproduce-fast about 8 to 12 min on a 10-core laptop.
#   make reproduce-full   the counts the claims need (1000 / 500 / 100 / 200).
#                         About 2 h 10 min on 64 cores (n2d-highcpu-64), from the stage
#                         times measured 2026-10-07: controls 28 min, transfer bootstrap
#                         1 h 37 min, change.py 3 min on 10 cores. The diffusivity controls
#                         (not rerun then) took 12 min on 32 cores.
# git checkout (or re-extracting the bundle) restores the committed results.

PY ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)
WORKERS ?= $(shell python3 -c "import os; print(os.cpu_count())")
export OMP_NUM_THREADS = 1

.PHONY: venv data reproduce-fast reproduce-full

venv:
	python3 -m venv .venv
	.venv/bin/pip install --quiet -r requirements.txt

data:
	@python3 -c "import hashlib;print(' '.join(hashlib.sha256(open(f,'rb').read()).hexdigest() for f in ('increments.json','quench.json','quench_variants.json')))" > .data.before
	$(PY) fetch_lsc.py
	$(PY) build_increments.py
	$(PY) quench.py
	@python3 -c "import hashlib,sys;a=open('.data.before').read().split();fs=('increments.json','quench.json','quench_variants.json');b=[hashlib.sha256(open(f,'rb').read()).hexdigest() for f in fs];[print(f, 'rebuilt identical' if x==y else 'CHANGED') for f,x,y in zip(fs,a,b)];sys.exit(a!=b)"
	@rm -f .data.before

reproduce-fast:
	$(PY) analysis.py --part fits --workers $(WORKERS)
	$(PY) analysis.py --part controls --reps 10 --workers $(WORKERS)
	$(PY) analysis.py --part transfer --boot-reps 20 --workers $(WORKERS)
	$(PY) change.py --reps 20 --workers $(WORKERS)
	$(PY) festim_check.py
	$(PY) diffusivity.py --part windows --workers $(WORKERS)
	$(PY) diffusivity.py --part controls --control-reps 4 --workers $(WORKERS)
	$(PY) make_figures.py
	$(PY) make_readme.py --draft

reproduce-full:
	$(PY) analysis.py --part fits --workers $(WORKERS)
	$(PY) analysis.py --part controls --reps 1000 --workers $(WORKERS)
	$(PY) analysis.py --part transfer --boot-reps 500 --workers $(WORKERS)
	$(PY) change.py --reps 200 --workers $(WORKERS)
	$(PY) festim_check.py
	$(PY) diffusivity.py --part windows --workers $(WORKERS)
	$(PY) diffusivity.py --part controls --control-reps 100 --workers $(WORKERS)
	$(PY) make_figures.py
	$(PY) make_readme.py
