"""Single source of truth for where things live.

Nothing else in the repo hardcodes a path. Data lives outside git (§0.2);
override the root with the ARDS_DATA env var if it sits on another disk.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

DATA = Path(os.environ.get("ARDS_DATA", REPO / "data"))
CONFIG = REPO / "config"
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
LOGS = REPO / "logs"
SCRATCH = REPO / "scratch"

# --- data sources ------------------------------------------------------------
# Directory names are asserted against, not assumed: see DOWNLOADS.md and
# layer0_scoping/01_scope_all_sources.py, which report what is actually present.
CHALLENGE = DATA / "challenge2019"          # PhysioNet/CinC Challenge 2019
EICU_DEMO = DATA / "eicu_demo"              # eICU-CRD Demo v2.0.1
MIMIC_DEMO = DATA / "mimic_iv_demo"         # MIMIC-IV Clinical Database Demo
GEO = DATA / "geo"                          # GSE65682, GSE32707, GSE66890
OLINK = DATA / "olink_mgh"                  # MGH COVID Olink Explore 1536

DUCKDB = DATA / "ards.duckdb"

for _d in (RESULTS, FIGURES, LOGS, SCRATCH):
    _d.mkdir(parents=True, exist_ok=True)
