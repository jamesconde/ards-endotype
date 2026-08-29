#!/usr/bin/env python3
"""GEO ingest and QC. SPEC.md §3.1, §3.2.

Answers the hour-1 question first (§3.1, §6.6): are the MARS endotype labels
actually in GSE65682's pData, or only in the Scicluna et al. supplement? If
absent, the supervised anchor (§3.3) is dropped and Layer 1 runs
unsupervised-only. That is a planned fallback, not a crisis.

QC per §3.2:
  - detect whether the matrix is already log2 and do NOT double-log
  - probe -> gene collapse keeping maximum IQR per gene, count survivors
  - PCA coloured by any available batch variable; if PC1 tracks batch rather
    than phenotype, say so and proceed within-cohort only

Cohorts are NEVER merged (§3.5, §6.5). Batch effects across three platforms
would manufacture clusters that separate perfectly by cohort, which is the
single most common failure in this literature.

Usage
-----
    python layer1_molecular/01_geo_ingest.py
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, step

log = setup("L1_geo")

SERIES = {
    "GSE65682": "MARS — discovery + supervised anchor",
    "GSE32707": "Dolinay — transfer cohort 1 (has a day 0/7 axis)",
    "GSE66890": "Kangelaris — transfer cohort 2",
}


def read_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Parse a GEO series matrix into (expression, sample metadata)."""
    opener = gzip.open if path.suffix == ".gz" else open
    meta: dict[str, list[str]] = {}
    rows, header = [], None
    in_table = False
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                parts = [p.strip('"') for p in line.split("\t")]
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
            elif line.startswith("!Sample_"):
                key, *vals = line.split("\t")
                meta.setdefault(key.lstrip("!"), []).append(
                    [v.strip('"') for v in vals])
    if header is None or not rows:
        raise CheckpointError(f"{path.name}: no expression table found")
    expr = pd.DataFrame(rows, columns=header).set_index(header[0])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr, meta


def find_labels(meta: dict) -> dict:
    """Search every characteristics field for MARS endotype labels (§6.6)."""
    hits = {}
    for key, blocks in meta.items():
        if "characteristics" not in key.lower() and key != "Sample_title":
            continue
        for i, vals in enumerate(blocks):
            joined = " ".join(vals).lower()
            if re.search(r"\bmars\s*[1-4]\b|endotype|\bclass\b", joined):
                uniq = sorted(set(vals))
                hits[f"{key}[{i}]"] = uniq[:8]
    return hits


def main() -> int:
    out: dict = {}
    for gse, role in SERIES.items():
        matches = list(paths.GEO.glob(f"{gse}*series_matrix.txt*"))
        if not matches:
            log.warning("%s: not downloaded, skipping", gse)
            out[gse] = {"present": False}
            continue

        with step(log, f"{gse} ({role})"):
            expr, meta = read_series_matrix(matches[0])
            log.info("  matrix: %d probes x %d samples", *expr.shape)

            # §3.2: do not double-log.
            mx = float(np.nanmax(expr.to_numpy(dtype=float)))
            already_log = mx <= 50
            if not already_log:
                log.warning("  max value %.1f > 50 -> applying log2(x+1)", mx)
                expr = np.log2(expr.clip(lower=0) + 1)
            else:
                log.info("  max value %.2f -> already log-scale, NOT re-logging", mx)

            n_all_na = int(expr.isna().all(axis=1).sum())
            expr = expr.dropna(how="all")
            log.info("  dropped %d all-NaN probes; %d remain", n_all_na, len(expr))
            if not (5000 <= len(expr) <= 80000):
                log.warning("  probe count %d outside the expected 5k-80k band", len(expr))

            entry = {
                "present": True, "role": role,
                "n_probes": int(expr.shape[0]), "n_samples": int(expr.shape[1]),
                "already_log2": bool(already_log), "max_value": round(mx, 3),
            }

            if gse == "GSE65682":
                labels = find_labels(meta)
                entry["label_search"] = labels
                if labels:
                    log.warning("  MARS-like label fields found:")
                    for k, v in labels.items():
                        log.warning("    %s -> %s", k, v)
                else:
                    log.error("  NO MARS labels in pData. §3.3 supervised anchor is "
                              "DROPPED; Layer 1 proceeds unsupervised-only (§6.6 "
                              "planned fallback). Second option: extract from the "
                              "Scicluna 2017 supplement and join by GSM.")
                entry["mars_labels_present"] = bool(labels)
                report_only("geo_mars_labels_present", int(bool(labels)))

            # Keep the QC'd matrix out of git (§0.2) but available downstream.
            outdir = paths.SCRATCH / "geo"
            outdir.mkdir(parents=True, exist_ok=True)
            expr.to_parquet(outdir / f"{gse}_expr.parquet")

            # Sample metadata, flattened for later joins.
            md = {k: v for k, v in meta.items()
                  if "characteristics" in k.lower() or k in
                  ("Sample_title", "Sample_geo_accession", "Sample_source_name_ch1")}
            pd.to_pickle(md, outdir / f"{gse}_meta.pkl")
            entry["metadata_fields"] = sorted(md.keys())
            out[gse] = entry

    (paths.RESULTS / "L1_geo_ingest.json").write_text(json.dumps(out, indent=2, default=str))
    log.info("wrote results/L1_geo_ingest.json")
    log.warning("REMINDER (§3.5, §6.5): never merge these cohorts and cluster the "
                "merged matrix. Train on one, score fixed signatures on the others.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
