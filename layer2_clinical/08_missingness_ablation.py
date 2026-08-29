#!/usr/bin/env python3
"""Deployment stress-test: how gracefully does it degrade? SPEC.md §4.5c, C5.

Three ablations, because they answer different questions:

  A. PROGRESSIVE REMOVAL — drop variables one at a time, least-available first,
     and measure agreement. This is the degradation curve.
  B. NATIVE-SPARSE vs NATIVE-COMPLETE — split patients by how much was actually
     measured on them and compare performance. No simulation involved.
  C. MASK-ONLY — train on the missingness pattern alone, discarding every
     value. If this performs well, missingness is carrying signal that a
     values-only model would attribute to physiology.

C is the one worth designing for. §4.4 predicted the mask "may carry more
signal than the values", and the only way to know is to fit it alone.

Removal order is by AVAILABILITY, not by importance. Dropping the least-useful
variable first flatters the curve; a real deployment loses whatever the site
does not measure, which is what availability encodes.

Usage
-----
    python layer2_clinical/08_missingness_ablation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import cohen_kappa_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step

log = setup("L2_ablation")

FEATURES = [
    "HCO3_mean", "Creatinine_mean", "WBC_mean", "Platelets_mean",
    "BUN_mean", "Glucose_mean", "pH_mean",
    "MAP_mean", "HR_mean", "Resp_mean", "Temp_mean", "sf_ratio_min",
]


def fit_score(X, y, seed) -> tuple[float, float]:
    """Out-of-fold kappa and AUC. NaN is passed through to the model.

    HistGradientBoosting handles NaN natively by learning a default branch
    direction per split, so no imputation is needed — which matters, because
    cohort-median imputation would put a population statistic into the inference
    path and quietly break portability (shared/no_population_stats).

    Note: plain GradientBoostingClassifier does NOT accept NaN. It appeared to
    work in an earlier version of this script only because that version ran on
    complete cases, where there were no NaNs to reject.
    """
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    m = HistGradientBoostingClassifier(max_iter=150, max_depth=3, random_state=seed)
    proba = cross_val_predict(m, X, y, cv=cv, method="predict_proba")[:, 1]
    pred = (proba >= 0.5).astype(int)
    return float(cohen_kappa_score(y, pred)), float(roc_auc_score(y, proba))


def main() -> int:
    cfg = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    seed = set_seeds(cfg["clustering"]["seed"])

    static = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")
    con = duckdb.connect(str(paths.DUCKDB))
    ref = con.execute(
        "SELECT patient_id, is_hyper_like_k2 FROM reference_assignment").df()
    con.close()
    df = static.merge(ref, on="patient_id", how="inner")
    log.info("labelled patients: %d", len(df))

    y = df.is_hyper_like_k2.to_numpy(int)

    # Availability MUST come from the full cohort, not the labelled subset.
    #
    # The labelled patients are exactly the complete cases the reference mixture
    # was fitted on, so within them every feature is present by construction and
    # every availability is 1.0. Computing the removal order there produced an
    # arbitrary order, and — worse — made the mask-only model degenerate: an
    # all-ones mask carries no information, so it scored AUC 0.499 and would have
    # been reported as "missingness is uninformative" when in fact the analysis
    # had conditioned the missingness away. Availability now comes from `static`.
    avail = {}
    for f in FEATURES:
        mask_col = f[:-5] + "_measured" if f.endswith("_mean") else None
        avail[f] = (float(static[mask_col].mean()) if mask_col in static.columns
                    else float(static[f].notna().mean()))
    order = sorted(FEATURES, key=lambda f: avail[f])   # least available first

    within_label = {f: float(df[f].notna().mean()) for f in FEATURES}
    if min(within_label.values()) > 0.999:
        log.warning(
            "Every feature is 100%% present among labelled patients, as expected: "
            "they ARE the reference model's complete cases. Ablations B and C "
            "therefore run on the FULL cohort against sepsis onset, which exists "
            "for every patient regardless of what was measured.")
    log.info("removal order (least available first): %s",
             [(f, round(avail[f], 3)) for f in order])

    # --- A. progressive removal ---------------------------------------------
    rows = []
    with step(log, "A: progressive variable removal"):
        remaining = list(FEATURES)
        while len(remaining) >= 2:
            X = df[remaining].to_numpy(float)
            k, auc = fit_score(X, y, seed)
            rows.append({
                "n_features": len(remaining),
                "kappa": round(k, 4), "auc": round(auc, 4),
                "removed_so_far": len(FEATURES) - len(remaining),
                "min_availability": round(min(avail[f] for f in remaining), 4),
                "features": ",".join(remaining),
            })
            log.info("  %2d features  kappa=%.3f  AUC=%.3f", len(remaining), k, auc)
            drop = next(f for f in order if f in remaining)
            remaining.remove(drop)
    curve = pd.DataFrame(rows)
    curve.to_csv(paths.RESULTS / "L2_ablation_curve.csv", index=False)

    full_k = curve.iloc[0].kappa
    half = curve[curve.kappa <= full_k * 0.5]
    n_half = int(half.n_features.max()) if len(half) else None
    log.warning("full-feature kappa %.3f; falls below half that at %s features",
                full_k, n_half)
    report_only("ablation_full_kappa", round(float(full_k), 4))

    # --- B & C on the FULL cohort, against an endpoint that always exists ----
    #
    # Endotype labels exist only for complete cases, so they cannot support any
    # analysis ABOUT missingness. Sepsis onset is recorded for all 5,687 cohort
    # patients whether or not a single lab was drawn, which makes it the only
    # target that can answer "does the ordering pattern itself carry signal?"
    meas_cols = [f[:-5] + "_measured" for f in FEATURES
                 if f[:-5] + "_measured" in static.columns]
    if not meas_cols:
        raise CheckpointError("no _measured columns; rerun 03_features.py")

    full = static.copy()
    full["n_measured"] = full[meas_cols].sum(axis=1)
    y_full = full.sepsis_any.to_numpy(int)
    log.info("full cohort: %d patients, measured-variable count ranges %d-%d "
             "(median %.0f)", len(full), int(full.n_measured.min()),
             int(full.n_measured.max()), full.n_measured.median())

    native = []
    with step(log, "B: native-sparse vs native-complete (full cohort)"):
        med = full.n_measured.median()
        for name, g in [("native_sparse", full[full.n_measured < med]),
                        ("native_complete", full[full.n_measured >= med])]:
            if len(g) < 100 or g.sepsis_any.nunique() < 2:
                log.warning("  %s: n=%d, skipping", name, len(g))
                continue
            k, auc = fit_score(g[FEATURES].to_numpy(float),
                               g.sepsis_any.to_numpy(int), seed)
            native.append({"group": name, "n": int(len(g)),
                           "mean_measured": round(float(g.n_measured.mean()), 2),
                           "sepsis_rate": round(float(g.sepsis_any.mean()), 4),
                           "kappa_sepsis": round(k, 4), "auc_sepsis": round(auc, 4)})
            log.info("  %-16s n=%4d  mean measured %.1f  sepsis %.1f%%  AUC=%.3f",
                     name, len(g), g.n_measured.mean(),
                     100 * g.sepsis_any.mean(), auc)
        pd.DataFrame(native).to_csv(
            paths.RESULTS / "L2_ablation_native_sparse.csv", index=False)

    with step(log, "C: mask-only vs values-only, predicting sepsis onset"):
        k_mask, auc_mask = fit_score(full[meas_cols].to_numpy(float), y_full, seed)
        k_val, auc_val = fit_score(full[FEATURES].to_numpy(float), y_full, seed)
        k_both, auc_both_ = fit_score(
            full[FEATURES + meas_cols].to_numpy(float), y_full, seed)
        log.warning("  values only : kappa=%.3f  AUC=%.3f", k_val, auc_val)
        log.warning("  MASK ONLY   : kappa=%.3f  AUC=%.3f", k_mask, auc_mask)
        log.warning("  values+mask : kappa=%.3f  AUC=%.3f", k_both, auc_both_)
        ratio = auc_mask / auc_val if auc_val else float("nan")
        log.warning("  the mask alone recovers %.0f%% of the values-only AUC", 100 * ratio)
        if auc_mask > auc_val:
            log.warning("  THE MASK BEATS THE VALUES. Which labs were ordered predicts "
                        "the endpoint better than what those labs said. §4.4 "
                        "anticipated this; it is the clearest possible statement that "
                        "missingness in ICU data encodes care process.")
        report_only("mask_only_auc_sepsis", round(auc_mask, 4),
                    "missingness pattern alone, full cohort, no values")

    out = {
        "removal_order_rationale": "least-available-first, not least-important-first",
        "availability": {k: round(v, 4) for k, v in avail.items()},
        "full_feature_kappa": round(float(full_k), 4),
        "n_features_at_half_kappa": n_half,
        "availability_within_labelled_subset": {k: round(v, 4) for k, v in within_label.items()},
        "native_sparse_vs_complete": native,
        "mask_vs_values_target": "sepsis onset (full cohort)",
        "mask_only": {"kappa": round(k_mask, 4), "auc": round(auc_mask, 4)},
        "values_only": {"kappa": round(k_val, 4), "auc": round(auc_val, 4)},
        "values_plus_mask": {"kappa": round(k_both, 4), "auc": round(auc_both_, 4)},
        "mask_share_of_values_auc": round(float(ratio), 4),
        "note": ("Ablation A runs on labelled patients, who are the reference "
                 "model's complete cases, so it measures loss from REMOVING "
                 "variables rather than from natural sparsity. B and C run on the "
                 "full cohort against sepsis onset, because endotype labels exist "
                 "only where data is complete and therefore cannot support any "
                 "analysis about missingness."),
    }
    (paths.RESULTS / "L2_ablation_summary.json").write_text(json.dumps(out, indent=2))
    log.info("wrote results/L2_ablation_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
