#!/usr/bin/env python3
"""Static and longitudinal features for the cohort. SPEC.md §4.4.

Static:       T0 .. T0+24h summaries.
Longitudinal: T0 .. T0+72h in 6-hour bins -> 12 timesteps.

**Every value channel gets a parallel binary missingness mask.** Whether a lab
was drawn encodes care process, not just absence of information, and
forward-filling without a mask fabricates trajectories. On this dataset the
labs are sparse enough that the mask may carry more signal than the values —
which is a reportable finding in its own right (§4.4), and is why the mask is
built here rather than left to the modelling step.

Nothing is imputed. Imputation is a modelling decision, and cohort-median
imputation in particular would put a population statistic into the inference
path (shared/no_population_stats.py).

Usage
-----
    python layer2_clinical/03_features.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, step

log = setup("L2_features")

VITALS = ["HR", "O2Sat", "Temp", "SBP", "MAP", "DBP", "Resp"]


def connect() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(str(paths.DUCKDB))
    con.execute("SET memory_limit='6GB'")
    con.execute(f"SET temp_directory='{paths.SCRATCH}'")
    return con


def build_static(con, labs: list[str], window: int) -> pd.DataFrame:
    """Per-patient summaries over T0 .. T0+window, plus a measured/not mask."""
    value_cols = labs + VITALS
    aggs = []
    for c in value_cols:
        aggs += [
            f"AVG({c}) AS {c}_mean",
            f"MIN({c}) AS {c}_min",
            f"MAX({c}) AS {c}_max",
            # The mask: was this ever measured in the window? (§4.4)
            f"CAST(COUNT({c}) > 0 AS INTEGER) AS {c}_measured",
            f"COUNT({c}) AS {c}_n_obs",
        ]
    con.execute(f"""
        CREATE OR REPLACE TABLE features_static AS
        SELECT c.patient_id, c.hospital_system, c.age, c.gender, c.hosp_adm_time,
               c.t0, c.min_sf, c.sepsis_any, c.n_hours,
               {", ".join(aggs)},
               AVG(h.O2Sat / NULLIF(h.FiO2, 0)) AS sf_ratio_mean,
               MIN(h.O2Sat / NULLIF(h.FiO2, 0)) AS sf_ratio_min
        FROM cohort c
        JOIN hourly h
          ON h.patient_id = c.patient_id
         AND h.ICULOS >= c.t0
         AND h.ICULOS <  c.t0 + {window}
        GROUP BY c.patient_id, c.hospital_system, c.age, c.gender,
                 c.hosp_adm_time, c.t0, c.min_sf, c.sepsis_any, c.n_hours
    """)
    df = con.execute("SELECT * FROM features_static").df()
    return df


def build_longitudinal(con, labs: list[str], window: int, bin_hours: int
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[int]]:
    """(n_patients, n_bins, n_features) value and mask tensors.

    Returns values with NaN where unobserved, and a separate binary mask, so a
    downstream model can consume either without the two ever being conflated.
    """
    value_cols = labs + VITALS
    n_bins = window // bin_hours
    aggs = ", ".join(
        [f"AVG({c}) AS {c}_mean" for c in value_cols]
        + [f"CAST(COUNT({c}) > 0 AS INTEGER) AS {c}_mask" for c in value_cols]
    )
    con.execute(f"""
        CREATE OR REPLACE TABLE features_long AS
        SELECT c.patient_id,
               CAST((h.ICULOS - c.t0) / {bin_hours} AS INTEGER) AS bin_idx,
               COUNT(*) AS n_rows_in_bin,
               {aggs}
        FROM cohort c
        JOIN hourly h
          ON h.patient_id = c.patient_id
         AND h.ICULOS >= c.t0
         AND h.ICULOS <  c.t0 + {window}
        GROUP BY c.patient_id, bin_idx
    """)
    long = con.execute("SELECT * FROM features_long ORDER BY patient_id, bin_idx").df()

    pids = sorted(long.patient_id.unique())
    pid_idx = {p: i for i, p in enumerate(pids)}
    vals = np.full((len(pids), n_bins, len(value_cols)), np.nan, dtype=np.float32)
    mask = np.zeros((len(pids), n_bins, len(value_cols)), dtype=np.int8)
    # in_record: was the patient still in the ICU during this bin at all?
    #
    # Without this, the mask conflates two completely different things: a lab
    # that was not ordered (informative missingness, the thing this project
    # studies) and a bin after the patient left the ICU (censoring, which
    # carries no information about ordering behaviour). Median stay is ~50h
    # against a 72h window, so a large share of late bins are censoring rather
    # than missingness, and a degradation analysis that mixed them would
    # attribute discharge to care process.
    in_record = np.zeros((len(pids), n_bins), dtype=np.int8)

    rows = long.bin_idx.between(0, n_bins - 1)
    sub = long[rows]
    ri = sub.patient_id.map(pid_idx).to_numpy()
    bi = sub.bin_idx.to_numpy().astype(int)
    in_record[ri, bi] = 1
    for j, c in enumerate(value_cols):
        vals[ri, bi, j] = sub[f"{c}_mean"].to_numpy(dtype=np.float32)
        mask[ri, bi, j] = sub[f"{c}_mask"].to_numpy(dtype=np.int8)

    return vals, mask, in_record, value_cols, pids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()

    p = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    idx = p["index_time"]
    # Extract the union of the §4.4 static list and the D12 transferable list.
    # This is not a change to the analysis plan: both lists were frozen before
    # any model was fitted, and the transferable set needs Potassium, Hgb and
    # Hct, which the §4.4 list happens not to contain. Extracting a superset
    # keeps one feature table serving both experiments.
    labs = list(dict.fromkeys(
        p["static_features"]["labs"] + p["transferable_features"]["labs"]))
    log.info("extracting %d labs (union of §4.4 static and D12 transferable)", len(labs))
    con = connect()

    n_cohort = con.execute("SELECT COUNT(*) FROM cohort").fetchone()[0]
    log.info("cohort: %d patients", n_cohort)

    with step(log, "static features") as r:
        static = build_static(con, labs, idx["static_window_hours"])
        r["n_out"] = len(static)
        if len(static) != n_cohort:
            raise CheckpointError(
                f"static features cover {len(static)} of {n_cohort} cohort patients — "
                "the join dropped patients, which means some have no rows in "
                "T0..T0+24h. Investigate before proceeding.")

        meas = [c for c in static.columns if c.endswith("_measured")]
        cov = static[meas].mean().sort_values()
        log.info("static feature availability (fraction of cohort with >=1 value):")
        for c, v in cov.items():
            log.info("    %-24s %.1%%" if False else "    %-24s %5.1f%%", c[:-9], 100 * v)

        # A feature nobody has is not a feature. Flag rather than silently keep.
        dead = cov[cov < 0.05]
        if len(dead):
            log.warning("features observed in <5%% of the cohort: %s", list(dead.index))

    with step(log, "longitudinal tensor") as r:
        vals, mask, in_record, cols, pids = build_longitudinal(
            con, labs, idx["longitudinal_window_hours"], idx["bin_hours"])
        r["n_out"] = len(pids)
        log.info("tensor: %s (patients x bins x features)", vals.shape)

        if vals.shape[1] != 12:
            raise CheckpointError(f"expected 12 timesteps, got {vals.shape[1]}")
        all_nan = [cols[j] for j in range(vals.shape[2])
                   if np.isnan(vals[:, :, j]).all()]
        if all_nan:
            raise CheckpointError(f"all-NaN feature columns in tensor: {all_nan}")
        if mask.min() == 1:
            raise CheckpointError(
                "mask channel is all ones — missingness was lost, and the mask is "
                "the point of this step (§4.4)")

        # Density must be computed WITHIN the record, or discharge is counted
        # as a missing measurement.
        density = mask.mean(axis=(0, 1))
        live = in_record.astype(bool)
        density_in = np.array([mask[:, :, j][live].mean() for j in range(mask.shape[2])])
        log.info("in-record bin coverage: %.1f%% of patient-bins "
                 "(the rest are post-discharge, not missing)", 100 * in_record.mean())
        log.info("observation density, sparsest first (raw vs within-record):")
        for j in np.argsort(density_in)[:8]:
            log.info("    %-20s raw %5.1f%%   within-record %5.1f%%",
                     cols[j], 100 * density[j], 100 * density_in[j])
        log.info("  overall: raw %.1f%%, within-record %.1f%%",
                 100 * mask.mean(), 100 * mask[live].mean())
        report_only("tensor_mask_density_raw", round(float(mask.mean()), 4))
        report_only("tensor_mask_density_in_record", round(float(mask[live].mean()), 4),
                    "the honest figure — excludes post-discharge bins")

    out = paths.SCRATCH / "features"
    out.mkdir(parents=True, exist_ok=True)
    static.to_parquet(out / "static.parquet", index=False)
    np.savez_compressed(out / "longitudinal.npz",
                        values=vals, mask=mask, in_record=in_record,
                        feature_names=np.array(cols), patient_ids=np.array(pids))
    log.info("wrote %s and %s", out / "static.parquet", out / "longitudinal.npz")

    # Committed summary only — aggregate statistics, never row-level (§8.4).
    summary = pd.DataFrame({
        "feature": cols,
        "bin_density_raw": density.round(4),
        "bin_density_in_record": density_in.round(4),
        "static_measured_frac": [
            round(float(static.get(f"{c}_measured", pd.Series([np.nan])).mean()), 4)
            for c in cols],
    }).sort_values("bin_density_in_record")
    summary.to_csv(paths.RESULTS / "L2_feature_availability.csv", index=False)
    log.info("\n%s", summary.to_string(index=False))

    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
