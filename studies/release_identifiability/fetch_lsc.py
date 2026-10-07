#!/usr/bin/env python3
"""Fetch the lab's raw LSC counter files for BABY-1L runs 1-4 and 6.

The processed release curves in upstream/ carry no counting uncertainty. The
per-vial counts do: every LSC export has the gross count rate (CPMA) and the
count time, so the Poisson error of each vial follows from the file itself.

This script downloads, at pinned refs, each run's
  data/general.json, data/processed_data.json,
  analysis/tritium/tritium_model.py, data/tritium_detection/*.csv
into .cache/release_identifiability/<repo>@<ref>/ with the lab's own folder
layout, so the lab's tritium_model.py can be executed unchanged against it.

Pins live in lsc_manifest.json next to this file. The first fetch writes it.
Every later fetch fails if a file's sha256 differs. For runs 1-4 the three
non-LSC files are also checked against upstream/MANIFEST.json, so this study
reads the same bytes as the rest of the repo.

Run 6 has no tag. It is pinned by commit sha (the pattern the repo already
uses for libra-toolbox).

Usage: python3 studies/release_identifiability/fetch_lsc.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CACHE = REPO / ".cache" / "release_identifiability"
MANIFEST = HERE / "lsc_manifest.json"
UPSTREAM_MANIFEST = REPO / "upstream" / "MANIFEST.json"

RUNS = {
    1: {"repo": "LIBRA-project/BABY-1L-run-1", "ref": "v0.6", "commit": "85b096f3021f30ed889af53a3645db4c547f6608"},
    2: {"repo": "LIBRA-project/BABY-1L-run-2", "ref": "v0.5", "commit": "67a2022589146c2260dc6892b764d0247962f262"},
    3: {"repo": "LIBRA-project/BABY-1L-run-3", "ref": "v0.2", "commit": "ba532615c673187e5ee1f41a0f6ac00d643c51a2"},
    4: {"repo": "LIBRA-project/BABY-1L-run-4", "ref": "v0.1", "commit": "cca792fb9b4bb3a0a1f49ef76ced9ebcb144f051"},
    6: {"repo": "LIBRA-project/BABY-1L-run-6", "ref": "82c5af9", "commit": "82c5af999e43d4855f1e1400c2c455c4fb10fe5e"},
}
FIXED_FILES = ["data/general.json", "data/processed_data.json", "analysis/tritium/tritium_model.py"]


def folder(run: int) -> Path:
    r = RUNS[run]
    return CACHE / f"{r['repo'].split('/')[1]}@{r['ref']}"


def list_lsc(run: int) -> list[str]:
    """LSC file names for a run. Taken from the committed pins when they exist,
    so a fresh clone needs no GitHub API call. Otherwise listed with the public
    GitHub contents API over plain HTTPS (no login, no gh CLI)."""
    if MANIFEST.exists():
        pinned = [f["path"].rsplit("/", 1)[1] for f in json.loads(MANIFEST.read_text())["files"]
                  if f["run"] == run and f["path"].startswith("data/tritium_detection/")]
        if pinned:
            return sorted(pinned)
    r = RUNS[run]
    url = f"https://api.github.com/repos/{r['repo']}/contents/data/tritium_detection?ref={r['commit']}"
    req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json",
                                               "User-Agent": "baby-benchmark-fetch"})
    items = json.loads(urllib.request.urlopen(req, timeout=60).read())
    return sorted(i["name"] for i in items if i["type"] == "file" and i["name"].lower().endswith(".csv"))


def main() -> int:
    old = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {"files": []}
    pinned = {(f["run"], f["path"]): f["sha256"] for f in old["files"]}
    up = json.loads(UPSTREAM_MANIFEST.read_text())["files"]
    up_sha = {(f["repo"], f["commit"], f["path"]): f["sha256"] for f in up}

    entries, bad = [], []
    for run, r in RUNS.items():
        paths = FIXED_FILES + [f"data/tritium_detection/{n}" for n in list_lsc(run)]
        for p in paths:
            url = f"https://raw.githubusercontent.com/{r['repo']}/{r['commit']}/{p}"
            data = urllib.request.urlopen(url, timeout=60).read()
            sha = hashlib.sha256(data).hexdigest()
            dest = folder(run) / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            if (run, p) in pinned and pinned[(run, p)] != sha:
                bad.append(f"run {run} {p}: sha changed")
            key = (r["repo"], r["commit"], p)
            if key in up_sha and up_sha[key] != sha:
                bad.append(f"run {run} {p}: differs from upstream/MANIFEST.json")
            entries.append({"run": run, "repo": r["repo"], "ref": r["ref"], "commit": r["commit"],
                            "path": p, "bytes": len(data), "sha256": sha,
                            "matches_upstream_manifest": (up_sha.get(key) == sha) if key in up_sha else None})
    if bad:
        print("\n".join(bad), file=sys.stderr)
        return 1
    MANIFEST.write_text(json.dumps({
        "_generated": "Written by fetch_lsc.py. Pins for the raw lab files this study reads. Do not edit by hand.",
        "files": entries}, indent=1) + "\n")
    n_up = sum(1 for e in entries if e["matches_upstream_manifest"])
    print(f"fetched {len(entries)} files into {CACHE}; {n_up} cross-checked against upstream/MANIFEST.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
