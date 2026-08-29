#!/usr/bin/env python3
"""The portable classifier — §4.5a2, the headline experiment.

Reproduce the reference assignment (§4.5a) using ONLY per-patient values at
inference: raw measurements, absolute thresholds, and ratios between a
patient's own variables. No z-scores against a cohort, no population means,
medians or quantiles anywhere in the inference path.

The claim is made executable rather than asserted (shared/no_population_stats):
  - static introspection rejects fitted scalers, imputers and PCA
  - a functional batch-invariance probe rejects anything that normalises across
    rows inside a custom predict function

Three candidate forms, in increasing order of deployability, per §4.5a2:
  1. logistic regression on raw, untransformed values
  2. a shallow decision tree with absolute cut-points — readable as a bedside rule
  3. gradient boosting on raw values — the performance ceiling

Reported: kappa against the reference, confusion matrix, calibration, and the
parsimony curve (§4.5a3) cross-referenced against variable availability, because
a variable present in 12% of patients is not deployable however much it helps.

Usage
-----
    python layer2_clinical/05_portable_classifier.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml
from sklearn.calibration import calibration_curve
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (brier_score_loss, cohen_kappa_score, confusion_matrix,
                             roc_auc_score)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.tree import DecisionTreeClassifier, export_text

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step
from shared.no_population_stats import (assert_batch_invariant,
                                        assert_no_population_stats,
                                        describe_inference_path)

log = setup("L2_portable")

# Raw, per-patient features. Deliberately the same variables the reference model
# used, so the comparison isolates PORTABILITY rather than information content.
FEATURES = [
    "HCO3_mean", "Creatinine_mean", "WBC_mean", "Platelets_mean",
    "BUN_mean", "Glucose_mean", "pH_mean",
    "MAP_mean", "HR_mean", "Resp_mean", "Temp_mean", "sf_ratio_min",
]


def load() -> pd.DataFrame:
    static = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")
    con = duckdb.connect(str(paths.DUCKDB))
    ref = con.execute("SELECT * FROM reference_assignment").df()
    con.close()
    df = static.merge(ref, on="patient_id", how="inner", suffixes=("", "_ref"))
    if df.empty:
        raise CheckpointError("no overlap between static features and reference_assignment")
    return df


def evaluate(name, model, X, y, feats, seed) -> dict:
    """Out-of-fold evaluation, then a portability audit of the refit model."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
    pred = (proba >= 0.5).astype(int)

    model.fit(X, y)
    # The portability claim, checked two ways (§4.5a2).
    assert_no_population_stats(model, name)
    assert_batch_invariant(lambda Z: model.predict_proba(Z)[:, 1], X, label=name)

    cm = confusion_matrix(y, pred)
    out = {
        "model": name,
        "n": int(len(y)),
        "kappa": round(float(cohen_kappa_score(y, pred)), 4),
        "auc": round(float(roc_auc_score(y, proba)), 4),
        "brier": round(float(brier_score_loss(y, proba)), 4),
        "accuracy": round(float((pred == y).mean()), 4),
        "confusion_matrix": cm.tolist(),
        "portable": True,
    }
    log.info("%-12s kappa=%.3f  AUC=%.3f  Brier=%.3f  acc=%.3f",
             name, out["kappa"], out["auc"], out["brier"], out["accuracy"])
    return out, proba


def parsimony_curve(X, y, feats, seed, availability, avail_ab) -> pd.DataFrame:
    """§4.5a3 — agreement as a function of variable count, greedily selected.

    Reported against BOTH variable count and variable availability: a variable
    that lifts kappa but is recorded for a small minority of patients is not a
    deployable gain.
    """
    remaining, chosen, rows = list(range(len(feats))), [], []
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    while remaining:
        best = None
        for j in remaining:
            cols = chosen + [j]
            m = DecisionTreeClassifier(max_depth=4, random_state=seed,
                                       min_samples_leaf=25)
            pr = cross_val_predict(m, X[:, cols], y, cv=cv, method="predict")
            k = cohen_kappa_score(y, pr)
            if best is None or k > best[0]:
                best = (k, j)
        kappa, j = best
        chosen.append(j)
        remaining.remove(j)
        ab = avail_ab.get(feats[j], {})
        rows.append({
            "n_features": len(chosen),
            "added": feats[j],
            "kappa": round(float(kappa), 4),
            "added_availability": round(float(availability.get(feats[j], np.nan)), 4),
            "avail_setA": round(float(ab.get("A", np.nan)), 4),
            "avail_setB": round(float(ab.get("B", np.nan)), 4),
            "transferable": bool(abs(ab.get("A", 0) - ab.get("B", 0)) < 0.10),
            "min_availability_in_set": round(
                float(min(availability.get(feats[c], 1.0) for c in chosen)), 4),
            "min_avail_setB_in_set": round(
                float(min(avail_ab.get(feats[c], {}).get("B", 1.0) for c in chosen)), 4),
        })
        log.info("  +%-18s n=%2d  kappa=%.3f   A=%3.0f%% B=%3.0f%%  %s",
                 feats[j], len(chosen), kappa, 100 * ab.get("A", float("nan")),
                 100 * ab.get("B", float("nan")),
                 "" if rows[-1]["transferable"] else "<-- NOT TRANSFERABLE")
    return pd.DataFrame(rows)


def main() -> int:
    p = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    seed = set_seeds(p["clustering"]["seed"])

    df = load()
    log.info("patients with both features and a reference label: %d", len(df))

    # Target: the k=2 hyperinflammatory-like class, which is the one comparable
    # to the published dichotomy (§4.5a).
    target = "is_hyper_like_k2"
    if target not in df.columns:
        raise CheckpointError(f"{target} missing — rerun 04_reference_assignment.py")

    sub = df.dropna(subset=FEATURES + [target])
    log.info("complete cases for the portable model: %d of %d", len(sub), len(df))
    X = sub[FEATURES].to_numpy(dtype=float)
    y = sub[target].to_numpy(dtype=int)
    log.info("class balance: %.1f%% hyperinflammatory-like", 100 * y.mean())

    # Availability, for the parsimony cross-reference (§4.5a3).
    static_all = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")

    def availability_of(f: str) -> float:
        """Fraction of the COHORT for which this feature has a value.

        The `_measured` mask is authoritative where it exists. Derived features
        such as sf_ratio_min have no mask, so availability is non-nullness.
        Guard against the earlier bug: a name that does not end in `_mean` used
        to fall through to `static_all[f].mean()`, which returned the mean S/F
        RATIO (≈150) as though it were a proportion.
        """
        mask_col = f[:-5] + "_measured" if f.endswith("_mean") else None
        if mask_col and mask_col in static_all.columns:
            return float(static_all[mask_col].mean())
        if f in static_all.columns:
            return float(static_all[f].notna().mean())
        return float("nan")

    availability = {f: availability_of(f) for f in FEATURES}

    # Per-hospital availability. Cohort-wide availability hides the thing that
    # actually decides deployability: a variable can be near-universal overall
    # and absent at one site (D12). This is the number a trial designer needs.
    def availability_by_system(f: str) -> dict[str, float]:
        mask_col = f[:-5] + "_measured" if f.endswith("_mean") else None
        col = mask_col if (mask_col and mask_col in static_all.columns) else None
        out = {}
        for sysname, g in static_all.groupby("hospital_system"):
            out[sysname] = float(g[col].mean()) if col else float(g[f].notna().mean())
        return out

    avail_ab = {f: availability_by_system(f) for f in FEATURES}
    bad = {k: v for k, v in availability.items() if not (0.0 <= v <= 1.0)}
    if bad:
        raise CheckpointError(
            f"availability must be a proportion in [0,1]; got {bad}. "
            "A value above 1 means a raw measurement was read as a fraction.")

    results, probas = [], {}
    with step(log, "portable candidates"):
        for name, model in [
            ("logistic", LogisticRegression(max_iter=5000, C=1.0)),
            ("tree", DecisionTreeClassifier(max_depth=4, min_samples_leaf=25,
                                            random_state=seed)),
            ("gbm", GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                               random_state=seed)),
        ]:
            res, proba = evaluate(name, model, X, y, FEATURES, seed)
            results.append(res)
            probas[name] = proba
            if name == "tree":
                rule = export_text(model, feature_names=FEATURES, decimals=1)
                (paths.RESULTS / "L2_portable_tree_rule.txt").write_text(rule)
                log.info("bedside rule written to results/L2_portable_tree_rule.txt\n%s",
                         rule[:900])
                log.info("inference path ships only: %s",
                         describe_inference_path(model, "tree"))

    # --- calibration, not just discrimination (§4.5d, C4) --------------------
    calib = []
    for name, proba in probas.items():
        frac_pos, mean_pred = calibration_curve(y, proba, n_bins=10, strategy="quantile")
        calib += [{"model": name, "bin": i,
                   "mean_predicted": round(float(mp), 4),
                   "observed_frac": round(float(fp), 4)}
                  for i, (mp, fp) in enumerate(zip(mean_pred, frac_pos))]
    pd.DataFrame(calib).to_csv(paths.RESULTS / "L2_portable_calibration.csv", index=False)

    with step(log, "parsimony curve"):
        curve = parsimony_curve(X, y, FEATURES, seed, availability, avail_ab)
        curve.to_csv(paths.RESULTS / "L2_parsimony_curve.csv", index=False)
        best_k = curve.kappa.max()
        within = curve[curve.kappa >= best_k - 0.02]
        smallest = int(within.n_features.min())
        log.info("smallest set within 0.02 kappa of the best: %d variables (kappa %.3f)",
                 smallest, float(within[within.n_features == smallest].kappa.iloc[0]))
        report_only("parsimony_smallest_set", smallest, "within 0.02 kappa of best")

        # The deployment question: of the variables the parsimonious model needs,
        # how many are actually recorded at the OTHER hospital?
        best_set = curve.head(smallest)
        broken = best_set[~best_set.transferable]
        log.warning(
            "PARSIMONY vs DEPLOYMENT: the %d-variable set is %s. %d of %d of its "
            "variables are not recorded comparably at set B: %s",
            smallest,
            "fully transferable" if broken.empty else "NOT deployable across sites",
            len(broken), smallest,
            {r.added: f"A={r.avail_setA:.0%} B={r.avail_setB:.0%}"
             for r in broken.itertuples()} or "none")
        report_only("parsimony_untransferable_vars", len(broken),
                    f"of the {smallest}-variable optimal set")

    pd.DataFrame(results).to_csv(paths.RESULTS / "L2_portable_performance.csv", index=False)
    (paths.RESULTS / "L2_portable_summary.json").write_text(json.dumps({
        "target": target,
        "n": int(len(sub)),
        "class_balance": round(float(y.mean()), 4),
        "features": FEATURES,
        "feature_availability": {k: round(v, 4) for k, v in availability.items()},
        "results": results,
        "parsimony_smallest_within_0.02": smallest,
        "availability_by_hospital": {k: {s2: round(v2, 4) for s2, v2 in v.items()}
                                     for k, v in avail_ab.items()},
        "portability": ("every model passed both static introspection and the "
                        "batch-invariance probe: no cohort statistic appears in "
                        "any inference path"),
    }, indent=2))
    log.info("wrote results/L2_portable_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
