#!/usr/bin/env python3
"""Reference endotype assignment — the deliberately UNUSABLE baseline. §4.5a.

Gaussian mixture on log-transformed, **cohort-standardised** static features.
It needs the cohort mean and SD of every variable at inference time, so it
cannot be deployed on one patient at a time. That is the point: it exists to
generate the reference labels that the portable classifier (§4.5a2) has to
reproduce, and to make the deployment problem concrete.

This is NOT presented as a result. It is the baseline whose deployment failure
the project is about.

Expected structure, per the published literature (Calfee 2014 and after): a
hyperinflammatory-like class with lower bicarbonate, higher lactate, higher
creatinine and bilirubin, and a higher event rate. Recovering that is a
REPRODUCTION of established work, and is labelled as such (I1).

k is selected by BIC over the pre-registered range in cohort_params.yaml.

Usage
-----
    python layer2_clinical/04_reference_assignment.py --fit-on A
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step

log = setup("L2_reference")

# Variables that are (a) available in most of the cohort and (b) the ones the
# published class definitions are built on. Chosen from observed availability
# in results/L2_feature_availability.csv, before any clustering was run.
CANDIDATE_FEATURES = [
    "HCO3_mean", "Creatinine_mean", "WBC_mean", "Platelets_mean",
    "BUN_mean", "Glucose_mean", "pH_mean",
    "MAP_mean", "HR_mean", "Resp_mean", "Temp_mean", "sf_ratio_min",
]
# Strictly positive variables get a log transform; pH and S/F do not.
LOG_FEATURES = {"Creatinine_mean", "WBC_mean", "Platelets_mean",
                "BUN_mean", "Glucose_mean", "HR_mean", "Resp_mean"}


def load_static() -> pd.DataFrame:
    f = paths.SCRATCH / "features" / "static.parquet"
    if not f.exists():
        raise CheckpointError(f"{f} missing — run 03_features.py first")
    return pd.read_parquet(f)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fit-on", default="A", choices=["A", "B", "both"],
                    help="hospital system to FIT on (§4.5d freezes this model)")
    args = ap.parse_args()

    p = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    cl = p["clustering"]
    seed = set_seeds(cl["seed"])

    static = load_static()
    log.info("static features: %s", static.shape)

    feats = [c for c in CANDIDATE_FEATURES if c in static.columns]
    missing = set(CANDIDATE_FEATURES) - set(feats)
    if missing:
        raise CheckpointError(f"expected features absent from static table: {missing}")

    # Complete cases only. Imputing here would inject a cohort statistic into a
    # model that is already non-portable, and would blur the very thing §4.5a2
    # is trying to measure. The cost is reported, not hidden.
    sub = static.dropna(subset=feats)
    log.info("complete cases on %d features: %d of %d (%.1f%%)",
             len(feats), len(sub), len(static), 100 * len(sub) / len(static))
    report_only("reference_complete_cases", len(sub),
                f"of {len(static)} cohort patients")
    if len(sub) < 500:
        raise CheckpointError(
            f"only {len(sub)} complete cases — too few to fit a stable mixture. "
            "Reduce the feature set rather than imputing.")

    # Who gets dropped by requiring complete data? This is a selection on top of
    # a selection (§6.2) and belongs in the writeup.
    dropped = static.loc[~static.index.isin(sub.index)]
    if len(dropped):
        sel = pd.DataFrame({
            "in_reference": [sub.age.mean(), sub.sepsis_any.mean(),
                             sub.n_hours.mean(), (sub.hospital_system == "A").mean()],
            "dropped": [dropped.age.mean(), dropped.sepsis_any.mean(),
                        dropped.n_hours.mean(), (dropped.hospital_system == "A").mean()],
        }, index=["age", "sepsis_rate", "iculos", "frac_setA"]).round(3)
        log.info("selection induced by requiring complete features:\n%s", sel.to_string())
        sel.to_csv(paths.RESULTS / "L2_reference_completecase_selection.csv")

    fit_mask = (sub.hospital_system == args.fit_on) if args.fit_on != "both" \
        else pd.Series(True, index=sub.index)
    fit = sub.loc[fit_mask]
    log.info("fitting on set %s: %d patients", args.fit_on, len(fit))

    X_raw = fit[feats].copy()
    for c in LOG_FEATURES & set(feats):
        X_raw[c] = np.log1p(X_raw[c].clip(lower=0))

    # THE NON-PORTABLE STEP, made explicit. These two vectors are cohort
    # statistics; a deployed model would not have them.
    mu, sd = X_raw.mean(), X_raw.std(ddof=0).replace(0, 1.0)
    X = ((X_raw - mu) / sd).to_numpy()
    log.warning("standardised against COHORT mean/SD — this model is not deployable "
                "by construction (§4.5a). It exists to label, not to ship.")

    with step(log, "BIC model selection"):
        rows = []
        best = None
        for k in range(min(cl["k_range"]), max(cl["k_range"]) + 1):
            gm = GaussianMixture(n_components=k, covariance_type="full",
                                 n_init=cl["n_restarts"], random_state=seed,
                                 max_iter=500)
            gm.fit(X)
            bic = gm.bic(X)
            rows.append({"k": k, "bic": round(bic, 1),
                         "converged": bool(gm.converged_)})
            log.info("  k=%d  BIC=%.1f  converged=%s", k, bic, gm.converged_)
            if best is None or bic < best[1]:
                best, best_gm = (k, bic), gm
        sel_df = pd.DataFrame(rows)
        sel_df.to_csv(paths.RESULTS / "L2_reference_bic.csv", index=False)
        k_best = best[0]
        log.info("BIC selects k=%d", k_best)
        report_only("reference_k_selected", k_best, "by BIC, not by eye")

        # A boundary solution is not a selection. If BIC is still falling at the
        # top of the pre-registered range, the data may prefer more classes and
        # the range — not the data — chose k. Say so rather than presenting k as
        # though the criterion had found an interior optimum.
        if k_best == max(cl["k_range"]):
            log.warning(
                "k=%d is the MAXIMUM of the pre-registered range %s and BIC was still "
                "falling. This is a boundary solution: the range constrained k, not the "
                "criterion. Report it as such. The range is NOT widened after seeing "
                "this — that would be selecting k by result (I4).", k_best, cl["k_range"])
            report_only("reference_k_at_range_boundary", 1,
                        "BIC still falling at k_max — range-constrained, not data-selected")

    labels = best_gm.predict(X)
    post = best_gm.predict_proba(X)
    sizes = pd.Series(labels).value_counts().sort_index()
    log.info("class sizes: %s", sizes.to_dict())
    if sizes.min() < 30:
        log.warning("a class has n=%d; §3.4 says a cluster with n<30 is not a finding",
                    sizes.min())

    # --- which class is the hyperinflammatory-like one? ----------------------
    # Decided by the published direction (lower HCO3), not by which looks nicer.
    prof = fit[feats].copy()
    prof["class"] = labels
    means = prof.groupby("class").mean()
    hyper = int(means["HCO3_mean"].idxmin()) if "HCO3_mean" in means else 0
    log.info("class %d has the lowest bicarbonate -> labelled hyperinflammatory-like "
             "on the published direction", hyper)

    means["n"] = sizes
    means["sepsis_rate"] = prof.assign(s=fit.sepsis_any.values).groupby("class").s.mean()
    means = means.round(3)
    means.to_csv(paths.RESULTS / "L2_reference_class_profiles.csv")
    log.info("\n%s", means.to_string())

    # Check the reproduction against the published direction, and say so plainly.
    d = means.loc[hyper] - means.drop(index=hyper).mean()
    directions = {
        "HCO3 lower": bool(d.get("HCO3_mean", 0) < 0),
        "Creatinine higher": bool(d.get("Creatinine_mean", 0) > 0),
        "sepsis rate higher": bool(d.get("sepsis_rate", 0) > 0),
        "WBC higher": bool(d.get("WBC_mean", 0) > 0),
    }
    log.info("agreement with the published hyper/hypo direction: %s", directions)
    n_agree = sum(directions.values())
    log.info("%d of %d published directions reproduced", n_agree, len(directions))

    # --- the published two-class structure, for comparison only --------------
    # NOT a second selection. BIC chose k above; this is fitted unconditionally
    # so that the reproduction can be compared against the published hyper/hypo
    # dichotomy on its own terms, and it is reported as a comparison.
    gm2 = GaussianMixture(n_components=2, covariance_type="full",
                          n_init=cl["n_restarts"], random_state=seed, max_iter=500).fit(X)
    lab2 = gm2.predict(X)
    prof2 = fit[feats].copy()
    prof2["class"] = lab2
    means2 = prof2.groupby("class").mean()
    hyper2 = int(means2["HCO3_mean"].idxmin())
    means2["n"] = pd.Series(lab2).value_counts().sort_index()
    means2["sepsis_rate"] = prof2.assign(s=fit.sepsis_any.values).groupby("class").s.mean()
    means2.round(3).to_csv(paths.RESULTS / "L2_reference_class_profiles_k2.csv")
    d2 = means2.loc[hyper2] - means2.drop(index=hyper2).mean()
    dir2 = {
        "HCO3 lower": bool(d2.get("HCO3_mean", 0) < 0),
        "Creatinine higher": bool(d2.get("Creatinine_mean", 0) > 0),
        "sepsis rate higher": bool(d2.get("sepsis_rate", 0) > 0),
        "WBC higher": bool(d2.get("WBC_mean", 0) > 0),
    }
    log.info("k=2 comparison (published structure): hyper-like class %d, n=%d/%d, "
             "directions %s", hyper2, int(means2.loc[hyper2, "n"]), len(lab2), dir2)
    hyper2_prev = float((lab2 == hyper2).mean())
    log.info("k=2 hyperinflammatory-like prevalence: %.1f%% "
             "(published cohorts report roughly 20-35%%)", 100 * hyper2_prev)

    out = pd.DataFrame({
        "patient_id": fit.patient_id.values,
        "reference_class_k2": lab2,
        "is_hyper_like_k2": (lab2 == hyper2).astype(int),
        "hospital_system": fit.hospital_system.values,
        "reference_class": labels,
        "is_hyper_like": (labels == hyper).astype(int),
        "max_posterior": post.max(axis=1),
    })
    con = duckdb.connect(str(paths.DUCKDB))
    con.execute("DROP TABLE IF EXISTS reference_assignment")
    con.execute("CREATE TABLE reference_assignment AS SELECT * FROM out")
    con.close()

    # Assignment confidence matters: a mixture that assigns everyone at p=0.51
    # has not found classes.
    log.info("posterior confidence: median %.3f, frac >0.9 = %.1f%%",
             np.median(post.max(axis=1)), 100 * (post.max(axis=1) > 0.9).mean())

    (paths.RESULTS / "L2_reference_summary.json").write_text(json.dumps({
        "fit_on": args.fit_on,
        "n_fit": int(len(fit)),
        "n_complete_cases": int(len(sub)),
        "n_cohort": int(len(static)),
        "features": feats,
        "k_selected": int(k_best),
        "class_sizes": {int(k): int(v) for k, v in sizes.items()},
        "hyper_like_class": hyper,
        "published_directions_reproduced": directions,
        "median_posterior": float(np.median(post.max(axis=1))),
        "k_at_range_boundary": bool(k_best == max(cl["k_range"])),
        "k2_comparison": {
            "hyper_like_class": hyper2,
            "hyper_like_prevalence": round(hyper2_prev, 4),
            "published_directions_reproduced": dir2,
        },
        "note": ("Reproduction of established subphenotype structure, not a "
                 "discovery. Cohort-standardised and therefore not deployable "
                 "by construction — see §4.5a and DECISIONS.md."),
    }, indent=2))
    log.info("wrote results/L2_reference_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
