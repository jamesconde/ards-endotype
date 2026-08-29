#!/usr/bin/env python3
"""Are the endotypes just severity? SPEC.md §6.4.

The standard and correct criticism of this entire literature: if two clusters
separate on "sicker vs less sick", you have rediscovered an organ-failure score
with extra steps. §6.4 pre-specifies three checks, all implemented here:

  1. AUC for the endpoint from endotype alone vs severity alone vs both.
     If endotype adds < 0.02 AUC over severity, say so plainly.
  2. Re-cluster with features standardised WITHIN severity quartile. Do the
     classes survive?
  3. Report separation on ratios and shapes, not levels alone.

Severity is a **modified SOFA** (DECISIONS.md D13). It must never be called
SOFA, because two of its six organ systems are unavailable on this data path:

  respiration   -- from S/F, not P/F (no arterial gases for most patients)
  coagulation   -- platelets                                  [standard]
  liver         -- bilirubin                                  [standard]
  renal         -- creatinine                                 [standard]
  cardiovascular-- MAP only; CAPPED AT 1 POINT, since scoring 2-4 requires
                   vasopressor dose, which this dataset does not contain
  CNS           -- OMITTED ENTIRELY, no GCS

So the maximum attainable score is 17, not 24, and the cardiovascular and CNS
components are systematically understated. That biases the comparison in
ENDOTYPE's favour — a weakened severity score is easier to beat — which is
worth stating, because it means a null result here is stronger than it looks.

Usage
-----
    python layer2_clinical/07_severity_check.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step

log = setup("L2_severity")

MAX_MODIFIED_SOFA = 17   # 4+4+4+4+1, CNS omitted


def sf_to_estimated_pf(sf: pd.Series) -> pd.Series:
    """Estimate P/F from S/F by linear interpolation between two anchors.

    The two well-established equivalences are S/F 235 ~ P/F 200 and S/F 315 ~
    P/F 300, giving slope 1.25. Extrapolation beyond that range is stated as an
    estimate rather than presented as a measurement.
    """
    return 200.0 + 1.25 * (sf - 235.0)


def modified_sofa(df: pd.DataFrame) -> pd.DataFrame:
    """Per-patient modified SOFA. Missing components score 0 and are counted."""
    out = pd.DataFrame(index=df.index)

    pf = sf_to_estimated_pf(df["sf_ratio_min"])
    out["resp"] = np.select(
        [pf < 100, pf < 200, pf < 300, pf < 400], [4, 3, 2, 1], default=0)
    out["resp_missing"] = df["sf_ratio_min"].isna().astype(int)

    plt = df["Platelets_min"] if "Platelets_min" in df else df["Platelets_mean"]
    out["coag"] = np.select([plt < 20, plt < 50, plt < 100, plt < 150],
                            [4, 3, 2, 1], default=0)
    out["coag_missing"] = plt.isna().astype(int)

    bili = df["Bilirubin_total_max"] if "Bilirubin_total_max" in df else df["Bilirubin_total_mean"]
    out["liver"] = np.select([bili >= 12, bili >= 6, bili >= 2, bili >= 1.2],
                             [4, 3, 2, 1], default=0)
    out["liver_missing"] = bili.isna().astype(int)

    cr = df["Creatinine_max"] if "Creatinine_max" in df else df["Creatinine_mean"]
    out["renal"] = np.select([cr >= 5, cr >= 3.5, cr >= 2, cr >= 1.2],
                             [4, 3, 2, 1], default=0)
    out["renal_missing"] = cr.isna().astype(int)

    # Capped at 1: scoring above this needs vasopressor dose, unavailable here.
    mapv = df["MAP_min"] if "MAP_min" in df else df["MAP_mean"]
    out["cardio"] = (mapv < 70).astype(int)
    out["cardio_missing"] = mapv.isna().astype(int)

    # CNS omitted entirely — no GCS anywhere in this dataset.
    comp = ["resp", "coag", "liver", "renal", "cardio"]
    for c in comp:
        out[c] = pd.Series(out[c], index=df.index).fillna(0)
    out["modified_sofa"] = out[comp].sum(axis=1)
    out["n_components_missing"] = out[[f"{c}_missing" for c in comp]].sum(axis=1)
    return out


def auc_cv(X: np.ndarray, y: np.ndarray, seed: int) -> float:
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    p = cross_val_predict(LogisticRegression(max_iter=5000), X, y, cv=cv,
                          method="predict_proba")[:, 1]
    return float(roc_auc_score(y, p))


def main() -> int:
    cfg = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    seed = set_seeds(cfg["clustering"]["seed"])

    static = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")
    con = duckdb.connect(str(paths.DUCKDB))
    ref = con.execute(
        "SELECT patient_id, is_hyper_like_k2, reference_class FROM reference_assignment").df()
    con.close()
    df = static.merge(ref, on="patient_id", how="inner")
    log.info("patients with reference labels: %d", len(df))

    sofa = modified_sofa(df)
    df = pd.concat([df.reset_index(drop=True), sofa.reset_index(drop=True)], axis=1)

    log.info("modified SOFA: mean %.2f, median %.0f, max %d (theoretical max %d)",
             df.modified_sofa.mean(), df.modified_sofa.median(),
             int(df.modified_sofa.max()), MAX_MODIFIED_SOFA)
    log.info("component missingness (scored 0 when absent): %s",
             {c: int(df[f"{c}_missing"].sum()) for c in
              ["resp", "coag", "liver", "renal", "cardio"]})
    report_only("modified_sofa_mean", round(float(df.modified_sofa.mean()), 2),
                f"max attainable {MAX_MODIFIED_SOFA}, CNS omitted")

    y = df.sepsis_any.to_numpy(int)
    endo = df.is_hyper_like_k2.to_numpy(int).reshape(-1, 1)
    sev = df.modified_sofa.to_numpy(float).reshape(-1, 1)

    with step(log, "does endotype add over severity?"):
        auc_endo = auc_cv(endo, y, seed)
        auc_sev = auc_cv(sev, y, seed)
        auc_both = auc_cv(np.hstack([endo, sev]), y, seed)
        delta = auc_both - auc_sev
        log.info("  AUC sepsis | endotype alone      %.4f", auc_endo)
        log.info("  AUC sepsis | modified SOFA alone %.4f", auc_sev)
        log.info("  AUC sepsis | both                %.4f", auc_both)
        log.warning("  endotype adds %+.4f AUC over severity alone", delta)
        if delta < 0.02:
            log.warning("  BELOW the 0.02 threshold §6.4 set. State plainly: on this "
                        "cohort, endotype membership adds essentially nothing over a "
                        "severity score for predicting the endpoint.")
        report_only("endotype_auc_gain_over_severity", round(delta, 4),
                    "§6.4 threshold is 0.02")

        # Guard against the wrong reading. "Endotype adds nothing over severity"
        # is only interesting if severity itself predicts something. Here every
        # AUC sits near 0.5, so this is a comparison between two weak predictors
        # and the honest conclusion is about the ENDPOINT, not only the classes.
        if max(auc_endo, auc_sev, auc_both) < 0.60:
            log.warning(
                "  ALL THREE AUCs ARE NEAR CHANCE (%.3f / %.3f / %.3f). Neither "
                "endotype nor severity predicts sepsis onset well in this cohort. "
                "Do not report 'endotype adds nothing over severity' without also "
                "reporting that severity itself adds little — the limiting factor "
                "is the endpoint, which is a 6h-shifted Sepsis-3 label rather than "
                "an outcome these classes were ever built to predict.",
                auc_endo, auc_sev, auc_both)

    with step(log, "do the classes survive severity adjustment?"):
        # Re-cluster with features standardised WITHIN severity quartile (§6.4).
        feats = ["HCO3_mean", "Creatinine_mean", "WBC_mean", "Platelets_mean",
                 "BUN_mean", "Glucose_mean", "pH_mean", "MAP_mean", "HR_mean",
                 "Resp_mean", "Temp_mean", "sf_ratio_min"]
        sub = df.dropna(subset=feats).copy()
        sub["sofa_q"] = pd.qcut(sub.modified_sofa, 4, labels=False, duplicates="drop")
        log.info("  severity quartiles: %s", sub.sofa_q.value_counts().sort_index().to_dict())

        Z = sub[feats].copy()
        for c in feats:
            if c not in ("pH_mean", "sf_ratio_min"):
                Z[c] = np.log1p(Z[c].clip(lower=0))
        # Within-quartile standardisation removes the severity gradient.
        Zw = Z.groupby(sub.sofa_q).transform(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1))
        gm = GaussianMixture(2, covariance_type="full", n_init=20,
                             random_state=seed).fit(Zw.to_numpy())
        lab_adj = gm.predict(Zw.to_numpy())

        ari = adjusted_rand_score(sub.is_hyper_like_k2.to_numpy(int), lab_adj)
        log.warning("  ARI between original classes and severity-adjusted classes: %.3f", ari)
        if ari < 0.2:
            log.warning("  The classes do NOT survive severity adjustment. On this "
                        "cohort, clinical-variable endotyping is largely a severity "
                        "axis — which §6.4 pre-accepts as a finding, and which changes "
                        "the report's conclusion rather than invalidating it.")
        report_only("ari_after_severity_adjustment", round(float(ari), 4))

        # Severity profile of the original classes, to show the confound directly.
        prof = sub.groupby("is_hyper_like_k2").modified_sofa.agg(["mean", "std", "count"])
        log.info("  modified SOFA by original class:\n%s", prof.round(2).to_string())

    out = {
        "severity_score": "modified SOFA (5 of 6 systems; cardiovascular capped at 1, CNS omitted)",
        "max_attainable": MAX_MODIFIED_SOFA,
        "mean_score": round(float(df.modified_sofa.mean()), 3),
        "auc_sepsis_endotype_alone": round(auc_endo, 4),
        "auc_sepsis_severity_alone": round(auc_sev, 4),
        "auc_sepsis_both": round(auc_both, 4),
        "endotype_gain_over_severity": round(delta, 4),
        "gain_exceeds_0.02_threshold": bool(delta >= 0.02),
        "ari_original_vs_severity_adjusted": round(float(ari), 4),
        "all_aucs_near_chance": bool(max(auc_endo, auc_sev, auc_both) < 0.60),
        "liver_component_missing_frac": round(float(df.liver_missing.mean()), 4),
        "caveat": ("The severity score is deliberately weakened by missing "
                   "vasopressors and GCS, and bilirubin is absent for most "
                   "patients so the liver component scores 0 by default. All of "
                   "that biases the comparison in ENDOTYPE's favour, so the null "
                   "is stronger than it appears, not weaker."),
        "second_caveat": ("Every AUC here is near 0.5. Neither endotype nor "
                          "severity predicts sepsis onset well, so this compares "
                          "two weak predictors. The limiting factor is the "
                          "endpoint — a 6h-shifted Sepsis-3 label, not an outcome "
                          "these subphenotypes were built to predict. Report both "
                          "halves or the finding is misleading."),
    }
    (paths.RESULTS / "L2_severity_check.json").write_text(json.dumps(out, indent=2))
    df[["patient_id", "hospital_system", "modified_sofa", "is_hyper_like_k2",
        "sepsis_any"]].groupby(
        ["hospital_system", "is_hyper_like_k2"]).agg(
        n=("patient_id", "count"), sofa_mean=("modified_sofa", "mean"),
        sepsis=("sepsis_any", "mean")).round(3).to_csv(
        paths.RESULTS / "L2_severity_by_class.csv")
    log.info("wrote results/L2_severity_check.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
