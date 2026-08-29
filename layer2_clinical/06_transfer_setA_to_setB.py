#!/usr/bin/env python3
"""Cross-hospital transfer, A -> B. SPEC.md §4.5d, restructured by DECISIONS.md D12.

WHAT THIS CAN AND CANNOT MEASURE — read before interpreting anything here.

The spec asked for "agreement of the portable classifier with a set-B-native
reference assignment". **That is impossible on this data.** The reference model
needs bicarbonate and chloride, which site B does not chart: only 16 of set B's
1,007 cohort patients have the variables it requires (D15). There is no
set-B-native reference to agree with, and no amount of modelling creates one.

So this script measures the three things that ARE measurable, and says so:

1. **Scoreability** — what fraction of set B can the frozen model score at all?
   This is the deployment question, and it is answered by counting, not fitting.
2. **Prevalence shift** — does the predicted class prevalence hold across sites?
3. **Outcome gradient** — does the predicted class still stratify a REAL,
   independently recorded endpoint (sepsis onset) at site B? This validates the
   transfer without needing endotype labels at B, because the endpoint exists
   whether or not the classifier can be scored.

Two feature sets, and the contrast between them is the case-mix versus
measurement-practice decomposition:

  RESTRICTED — the D12 transferable set, coverage matched across sites by
               construction. A drop here is genuine case mix / distribution shift.
  NAIVE      — the full set including the blood-gas block. Expected to fail on
               scoreability. That failure is the deployment reality, not a bug.

Models are FROZEN after fitting on A. Nothing is refit on B.

Usage
-----
    python layer2_clinical/06_transfer_setA_to_setB.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step
from shared.no_population_stats import assert_batch_invariant, assert_no_population_stats

log = setup("L2_transfer")

FULL_FEATURES = [
    "HCO3_mean", "Creatinine_mean", "WBC_mean", "Platelets_mean",
    "BUN_mean", "Glucose_mean", "pH_mean",
    "MAP_mean", "HR_mean", "Resp_mean", "Temp_mean", "sf_ratio_min",
]


def restricted_features(cfg: dict, columns) -> list[str]:
    tf = cfg["transferable_features"]
    wanted = [f"{c}_mean" for c in tf["vitals"] + tf["labs"]]
    have = [c for c in wanted if c in columns]
    missing = set(wanted) - set(have)
    if missing:
        log.warning("transferable features absent from the static table: %s", missing)
    return have


def wilson(k: int, n: int) -> tuple[float, float]:
    """Wilson interval — honest at small n, which set B is."""
    if n == 0:
        return (float("nan"), float("nan"))
    lo, hi = stats.binomtest(k, n).proportion_ci(confidence_level=0.95)
    return float(lo), float(hi)


def outcome_gradient(df: pd.DataFrame, pred_col: str, label: str) -> dict:
    """Does the predicted class stratify sepsis onset? Works without endotype labels."""
    out = {"site": label, "n": int(len(df))}
    for cls in (0, 1):
        g = df[df[pred_col] == cls]
        k = int(g.sepsis_any.sum())
        lo, hi = wilson(k, len(g))
        out[f"class{cls}_n"] = int(len(g))
        out[f"class{cls}_sepsis"] = round(float(g.sepsis_any.mean()), 4) if len(g) else None
        out[f"class{cls}_ci"] = [round(lo, 4), round(hi, 4)]
    if out["class0_n"] and out["class1_n"]:
        a, b = df[df[pred_col] == 1].sepsis_any, df[df[pred_col] == 0].sepsis_any
        out["risk_ratio"] = round(float(a.mean() / b.mean()), 3) if b.mean() else None
        out["auc_sepsis"] = round(float(roc_auc_score(df.sepsis_any, df[pred_col])), 4) \
            if df.sepsis_any.nunique() > 1 else None
        # A preserved DIRECTION is not a preserved EFFECT. Set B is small, so the
        # test is reported alongside the ratio to stop the two being conflated.
        k1, n1 = int(a.sum()), int(len(a))
        k0, n0 = int(b.sum()), int(len(b))
        p = float(stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])[1])
        out["fisher_p"] = round(p, 4)
        out["significant_at_0.05"] = bool(p < 0.05)
    return out


def main() -> int:
    cfg = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    seed = set_seeds(cfg["clustering"]["seed"])

    static = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")
    con = duckdb.connect(str(paths.DUCKDB))
    ref = con.execute(
        "SELECT patient_id, is_hyper_like_k2 FROM reference_assignment").df()
    con.close()

    df = static.merge(ref, on="patient_id", how="left")
    A = df[df.hospital_system == "A"]
    B = df[df.hospital_system == "B"]
    log.info("set A: %d cohort patients (%d with reference labels)",
             len(A), int(A.is_hyper_like_k2.notna().sum()))
    log.info("set B: %d cohort patients (%d with reference labels)",
             len(B), int(B.is_hyper_like_k2.notna().sum()))

    n_b_ref = int(B.is_hyper_like_k2.notna().sum())
    report_only("setB_with_reference_labels", n_b_ref, f"of {len(B)} — see D15")
    if n_b_ref > 100:
        log.warning("set B now has %d reference labels; the D15 premise has changed "
                    "and a native-reference comparison may be possible", n_b_ref)

    feature_sets = {
        "restricted": restricted_features(cfg, static.columns),
        "naive": [f for f in FULL_FEATURES if f in static.columns],
    }
    for k, v in feature_sets.items():
        log.info("%s feature set: %d variables — %s", k, len(v), v)

    results, rows_gradient = [], []
    for setname, feats in feature_sets.items():
        with step(log, f"transfer [{setname}]"):
            train = A.dropna(subset=feats + ["is_hyper_like_k2"])
            if len(train) < 200:
                raise CheckpointError(
                    f"{setname}: only {len(train)} labelled complete cases on set A")
            Xtr = train[feats].to_numpy(float)
            ytr = train.is_hyper_like_k2.to_numpy(int)

            model = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                               random_state=seed).fit(Xtr, ytr)
            assert_no_population_stats(model, f"{setname}_gbm")
            assert_batch_invariant(lambda Z: model.predict_proba(Z)[:, 1], Xtr,
                                   label=f"{setname}_gbm")
            log.info("  fitted on %d set-A patients, then FROZEN", len(train))

            # --- 1. scoreability, the deployment question -------------------
            scoreable_A = A.dropna(subset=feats)
            scoreable_B = B.dropna(subset=feats)
            frac_A = len(scoreable_A) / len(A)
            frac_B = len(scoreable_B) / len(B)
            log.info("  scoreable: A %d/%d (%.1f%%)   B %d/%d (%.1f%%)",
                     len(scoreable_A), len(A), 100 * frac_A,
                     len(scoreable_B), len(B), 100 * frac_B)

            entry = {
                "feature_set": setname,
                "n_features": len(feats),
                "n_train_A": int(len(train)),
                "scoreable_A": int(len(scoreable_A)),
                "scoreable_B": int(len(scoreable_B)),
                "frac_scoreable_A": round(frac_A, 4),
                "frac_scoreable_B": round(frac_B, 4),
            }

            if len(scoreable_B) < 30:
                log.error("  set B is effectively UNSCOREABLE with this feature set "
                          "(%d patients). This is the result, not a failure to report "
                          "around.", len(scoreable_B))
                entry.update({"prevalence_A": None, "prevalence_B": None,
                              "unscoreable": True})
                results.append(entry)
                continue

            # --- 2. prevalence shift ----------------------------------------
            pA = model.predict(scoreable_A[feats].to_numpy(float))
            pB = model.predict(scoreable_B[feats].to_numpy(float))
            entry["prevalence_A"] = round(float(pA.mean()), 4)
            entry["prevalence_B"] = round(float(pB.mean()), 4)
            entry["unscoreable"] = False
            log.info("  predicted hyper-like prevalence: A %.1f%%  B %.1f%%",
                     100 * pA.mean(), 100 * pB.mean())

            # --- 3. outcome gradient on a real, independent endpoint --------
            a_df = scoreable_A.assign(pred=pA)
            b_df = scoreable_B.assign(pred=pB)
            gA = outcome_gradient(a_df, "pred", f"A [{setname}]")
            gB = outcome_gradient(b_df, "pred", f"B [{setname}]")
            rows_gradient += [gA, gB]
            for g in (gA, gB):
                log.info("  %-18s sepsis: class0 %.1f%% (n=%d)  class1 %.1f%% (n=%d)  "
                         "RR=%s  AUC=%s  p=%s %s", g["site"],
                         100 * (g["class0_sepsis"] or 0), g["class0_n"],
                         100 * (g["class1_sepsis"] or 0), g["class1_n"],
                         g.get("risk_ratio"), g.get("auc_sepsis"), g.get("fisher_p"),
                         "" if g.get("significant_at_0.05") else "<-- NOT SIGNIFICANT")
            entry["auc_sepsis_A"] = gA.get("auc_sepsis")
            entry["auc_sepsis_B"] = gB.get("auc_sepsis")
            entry["rr_A"] = gA.get("risk_ratio")
            entry["rr_B"] = gB.get("risk_ratio")
            entry["p_A"] = gA.get("fisher_p")
            entry["p_B"] = gB.get("fisher_p")
            results.append(entry)

    res = pd.DataFrame(results)
    res.to_csv(paths.RESULTS / "L2_transfer_results.csv", index=False)
    pd.DataFrame(rows_gradient).to_csv(
        paths.RESULTS / "L2_transfer_outcome_gradient.csv", index=False)
    log.info("\n%s", res.to_string(index=False))

    # --- the decomposition ---------------------------------------------------
    r = res.set_index("feature_set")
    log.warning("=" * 74)
    log.warning("DECOMPOSITION — what the A->B drop is actually made of")
    if "naive" in r.index and "restricted" in r.index:
        log.warning("  scoreable at B:  restricted %.1f%%   naive %.1f%%",
                    100 * r.loc["restricted", "frac_scoreable_B"],
                    100 * r.loc["naive", "frac_scoreable_B"])
        log.warning("  The gap between those two numbers is MEASUREMENT PRACTICE.")
        log.warning("  Any residual change within the restricted set is CASE MIX, "
                    "because its coverage is matched across sites by construction.")
    log.warning("=" * 74)

    (paths.RESULTS / "L2_transfer_summary.json").write_text(json.dumps({
        "design": "fit on set A, frozen, applied to set B; no refitting on B",
        "what_cannot_be_measured": (
            "agreement against a set-B-native reference assignment — set B has "
            f"{n_b_ref} of {len(B)} patients with the variables the reference model "
            "needs (D15), so no such reference exists"),
        "feature_sets": {k: v for k, v in feature_sets.items()},
        "results": results,
        "outcome_gradient": rows_gradient,
    }, indent=2))
    log.info("wrote results/L2_transfer_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
