#!/usr/bin/env python3
"""Iterative refinement loop. SPEC.md §4.7, C6.

The JD line nothing else in this project covers: pipelines where predictions,
observed outcomes and new data feed back into subsequent model refinement.

THE STRUCTURE, which is the deliverable more than the performance change is:

  Round 1   train on set A, predict on a DIAGNOSIS slice of set B
  Confront  compare predictions with observed outcomes; build an error taxonomy
            stratified by missingness, case mix and predicted class
  Derive    a refinement traceable to a SPECIFIC round-1 finding, chosen by a
            PRE-SPECIFIED rule rather than by inspecting the holdout
  Round 2   retrain with the refinement, evaluate on a FINAL HOLDOUT that was
            never touched during diagnosis
  Report    the delta honestly, including if it is null or negative

GUARD AGAINST THE OBVIOUS FAILURE. This is a loop, not a licence to tune on the
test set. Set B is split once, at the top, with a fixed seed. The holdout is
scored exactly once, at the end. The refinement rules below are written down
BEFORE any diagnosis runs, so the choice of refinement cannot be steered by
what would help the holdout.

THE TASK BEING REFINED is enrichment quality: does the predicted-positive group
have a materially higher event rate than the predicted-negative group? That is
what a trial screening deployment actually needs (§4.6), and it is measurable
on set B without endotype labels, which set B does not have (D15).

Usage
-----
    python layer2_clinical/11_refinement_loop.py
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step

log = setup("L2_refinement")

HOLDOUT_FRAC = 0.20


@dataclass
class RefinementRound:
    """One cycle of the loop, kept as data so the architecture is legible."""
    round_id: int
    description: str
    features: list[str]
    threshold: float
    trained_on: str
    evaluated_on: str
    n_train: int = 0
    n_eval: int = 0
    event_rate_positive: float = float("nan")
    event_rate_negative: float = float("nan")
    risk_ratio: float = float("nan")
    auc: float = float("nan")
    fisher_p: float = float("nan")
    triggering_observation: str = ""
    refinement_applied: str = "none (baseline round)"

    def evaluate(self, y_pred: np.ndarray, y_prob: np.ndarray,
                 y_true: np.ndarray) -> "RefinementRound":
        self.n_eval = int(len(y_true))
        pos, neg = y_true[y_pred == 1], y_true[y_pred == 0]
        self.event_rate_positive = round(float(pos.mean()), 4) if len(pos) else float("nan")
        self.event_rate_negative = round(float(neg.mean()), 4) if len(neg) else float("nan")
        if len(pos) and len(neg) and neg.mean() > 0:
            self.risk_ratio = round(float(pos.mean() / neg.mean()), 3)
            self.fisher_p = round(float(stats.fisher_exact(
                [[int(pos.sum()), len(pos) - int(pos.sum())],
                 [int(neg.sum()), len(neg) - int(neg.sum())]])[1]), 4)
        if len(np.unique(y_true)) > 1:
            self.auc = round(float(roc_auc_score(y_true, y_prob)), 4)
        return self


# ---------------------------------------------------------------------------
# PRE-SPECIFIED REFINEMENT RULES — written before any diagnosis was run.
# Each names its trigger, its action, and the round-1 observation that would
# fire it. Nothing here may be edited after the diagnosis is seen.
# ---------------------------------------------------------------------------
RULES = [
    {"id": "R1",
     "trigger": "predicted-positive prevalence differs from set A by more than 5 points",
     "action": "recalibrate the decision threshold so predicted prevalence on B "
               "matches the prevalence the model produced on A"},
    {"id": "R2",
     "trigger": "enrichment (risk ratio) differs by more than 0.3 between "
                "natively-sparse and natively-complete patients",
     "action": "add the count of measured variables as an explicit feature, so "
               "the model can condition on care intensity rather than confounding "
               "it with physiology"},
    {"id": "R3",
     "trigger": "any input feature is available for less than 50% of set B",
     "action": "drop that feature and retrain, since it cannot be relied on at "
               "the deployment site"},
]


def fit_predict(train: pd.DataFrame, evalset: pd.DataFrame, feats: list[str],
                threshold: float, seed: int):
    tr = train.dropna(subset=feats + ["is_hyper_like_k2"])
    m = HistGradientBoostingClassifier(max_iter=200, max_depth=3,
                                       random_state=seed).fit(
        tr[feats].to_numpy(float), tr.is_hyper_like_k2.to_numpy(int))
    prob = m.predict_proba(evalset[feats].to_numpy(float))[:, 1]
    return m, (prob >= threshold).astype(int), prob, len(tr)


def main() -> int:
    cfg = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    seed = set_seeds(cfg["clustering"]["seed"])

    static = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")
    con = duckdb.connect(str(paths.DUCKDB))
    ref = con.execute(
        "SELECT patient_id, is_hyper_like_k2 FROM reference_assignment").df()
    con.close()
    df = static.merge(ref, on="patient_id", how="left")

    tf = cfg["transferable_features"]
    feats0 = [f"{c}_mean" for c in tf["vitals"] + tf["labs"]
              if f"{c}_mean" in static.columns]
    # Measured-count must span ALL labs, not just the model's inputs. Set B is
    # filtered to complete cases on feats0, so a count over feats0 alone is
    # constant by construction and the sparse/complete stratification would be
    # empty — which silently disables rule R2 rather than testing it.
    meas_cols = [c for c in static.columns if c.endswith("_measured")]
    model_meas = [f.replace("_mean", "_measured") for f in feats0
                  if f.replace("_mean", "_measured") in static.columns]
    log.info("measured-count spans %d lab/vital masks (model uses %d of them)",
             len(meas_cols), len(model_meas))

    A = df[df.hospital_system == "A"]
    B = df[df.hospital_system == "B"].dropna(subset=feats0).copy()
    B["n_measured"] = B[meas_cols].sum(axis=1)

    # THE SPLIT. Once, here, before anything is looked at.
    rng = np.random.default_rng(seed)
    is_holdout = rng.random(len(B)) < HOLDOUT_FRAC
    B_diag, B_hold = B[~is_holdout], B[is_holdout]
    log.info("set B scoreable: %d -> diagnosis %d, FINAL HOLDOUT %d (touched once)",
             len(B), len(B_diag), len(B_hold))
    if len(B_hold) < 80:
        raise CheckpointError(f"holdout is only {len(B_hold)} patients — too small")

    rounds: list[RefinementRound] = []
    decisions: list[dict] = []

    # ---------------- ROUND 1 ------------------------------------------------
    with step(log, "round 1: train on A, predict on B-diagnosis"):
        r1 = RefinementRound(1, "baseline: transferable features, threshold 0.5",
                             feats0, 0.5, "set A", "set B (diagnosis slice)")
        _, pred, prob, n_tr = fit_predict(A, B_diag, feats0, 0.5, seed)
        r1.n_train = n_tr
        r1.evaluate(pred, prob, B_diag.sepsis_any.to_numpy(int))
        rounds.append(r1)
        log.info("  positive rate %.1f%%, event rate pos %.1f%% vs neg %.1f%%, "
                 "RR=%.2f, AUC=%.3f, p=%s", 100 * pred.mean(),
                 100 * r1.event_rate_positive, 100 * r1.event_rate_negative,
                 r1.risk_ratio, r1.auc, r1.fisher_p)

        # prevalence the same model produces on its own training site
        _, predA, _, _ = fit_predict(A, A.dropna(subset=feats0), feats0, 0.5, seed)
        prevalence_A = float(predA.mean())
        prevalence_B = float(pred.mean())

    # ---------------- CONFRONT: error taxonomy -------------------------------
    with step(log, "confront predictions with outcomes"):
        diag = B_diag.assign(pred=pred, prob=prob)
        med = diag.n_measured.median()
        strata = {
            "natively_sparse": diag[diag.n_measured < med],
            "natively_complete": diag[diag.n_measured >= med],
            "younger_than_65": diag[diag.age < 65],
            "65_and_older": diag[diag.age >= 65],
        }
        tax = []
        for name, g in strata.items():
            pos, neg = g[g.pred == 1].sepsis_any, g[g.pred == 0].sepsis_any
            rr = float(pos.mean() / neg.mean()) if len(neg) and neg.mean() > 0 else np.nan
            tax.append({"stratum": name, "n": int(len(g)),
                        "positive_rate": round(float(g.pred.mean()), 4),
                        "event_pos": round(float(pos.mean()), 4) if len(pos) else None,
                        "event_neg": round(float(neg.mean()), 4) if len(neg) else None,
                        "risk_ratio": round(rr, 3) if np.isfinite(rr) else None})
            log.info("  %-20s n=%4d  RR=%s", name, len(g), tax[-1]["risk_ratio"])
        pd.DataFrame(tax).to_csv(paths.RESULTS / "L2_refinement_error_taxonomy.csv",
                                 index=False)

        avail_B = {f: float(B_diag[f].notna().mean()) for f in feats0}
        rr_sparse = next(t["risk_ratio"] for t in tax if t["stratum"] == "natively_sparse")
        rr_complete = next(t["risk_ratio"] for t in tax if t["stratum"] == "natively_complete")

    # ---------------- DERIVE: apply pre-specified rules ----------------------
    with step(log, "derive refinement from round-1 findings"):
        applied, new_feats, new_threshold = [], list(feats0), 0.5

        prev_gap = abs(prevalence_B - prevalence_A)
        if prev_gap > 0.05:
            new_threshold = float(np.quantile(prob, 1 - prevalence_A))
            applied.append({
                "rule": "R1", "fired": True,
                "observation": f"predicted-positive prevalence {prevalence_B:.1%} on B "
                               f"vs {prevalence_A:.1%} on A, gap {prev_gap:.1%}",
                "action": f"threshold moved 0.5 -> {new_threshold:.3f} to match A's prevalence"})
        else:
            applied.append({"rule": "R1", "fired": False,
                            "observation": f"prevalence gap {prev_gap:.1%} within 5 points"})

        rr_gap = (abs(rr_sparse - rr_complete)
                  if (rr_sparse and rr_complete) else 0.0)
        if rr_gap > 0.3:
            new_feats = feats0 + ["n_measured_feat"]
            applied.append({
                "rule": "R2", "fired": True,
                "observation": f"risk ratio {rr_sparse} in sparse vs {rr_complete} in "
                               f"complete patients, gap {rr_gap:.2f}",
                "action": "added measured-variable count as an explicit feature"})
        else:
            applied.append({"rule": "R2", "fired": False,
                            "observation": f"risk-ratio gap {rr_gap:.2f} within 0.3"})

        low = [f for f, a in avail_B.items() if a < 0.50]
        if low:
            new_feats = [f for f in new_feats if f not in low]
            applied.append({"rule": "R3", "fired": True,
                            "observation": f"features below 50% availability on B: {low}",
                            "action": f"dropped {low}"})
        else:
            applied.append({"rule": "R3", "fired": False,
                            "observation": "all features available in >=50% of set B"})

        for a in applied:
            log.info("  %s fired=%s: %s", a["rule"], a["fired"], a["observation"])
        decisions = applied
        if not any(a["fired"] for a in applied):
            log.warning("  NO RULE FIRED. Round 2 will be identical to round 1, and the "
                        "honest report is that diagnosis found nothing worth changing.")

    # ---------------- ROUND 2 on the FINAL HOLDOUT ---------------------------
    with step(log, "round 2: evaluate on the final holdout (touched once)"):
        A2, Bh = A.copy(), B_hold.copy()
        if "n_measured_feat" in new_feats:
            A2["n_measured_feat"] = A2[meas_cols].sum(axis=1)
            Bh["n_measured_feat"] = Bh[meas_cols].sum(axis=1)

        # Round 1 baseline, scored on the SAME holdout, for a fair delta.
        r1h = RefinementRound(1, "baseline scored on final holdout", feats0, 0.5,
                              "set A", "set B (FINAL HOLDOUT)")
        _, p1, pr1, n1 = fit_predict(A, Bh, feats0, 0.5, seed)
        r1h.n_train = n1
        r1h.evaluate(p1, pr1, Bh.sepsis_any.to_numpy(int))

        r2 = RefinementRound(2, "refined per pre-specified rules", new_feats,
                             new_threshold, "set A", "set B (FINAL HOLDOUT)")
        r2.triggering_observation = "; ".join(
            a["observation"] for a in decisions if a["fired"]) or "no rule fired"
        r2.refinement_applied = "; ".join(
            a.get("action", "") for a in decisions if a["fired"]) or "none"
        _, p2, pr2, n2 = fit_predict(A2, Bh, new_feats, new_threshold, seed)
        r2.n_train = n2
        r2.evaluate(p2, pr2, Bh.sepsis_any.to_numpy(int))
        rounds += [r1h, r2]

        log.warning("  holdout n=%d", len(Bh))
        log.warning("  round 1: RR=%.2f  AUC=%.3f  event pos %.1f%% / neg %.1f%%  p=%s",
                    r1h.risk_ratio, r1h.auc, 100 * r1h.event_rate_positive,
                    100 * r1h.event_rate_negative, r1h.fisher_p)
        log.warning("  round 2: RR=%.2f  AUC=%.3f  event pos %.1f%% / neg %.1f%%  p=%s",
                    r2.risk_ratio, r2.auc, 100 * r2.event_rate_positive,
                    100 * r2.event_rate_negative, r2.fisher_p)
        d_rr = r2.risk_ratio - r1h.risk_ratio
        d_auc = r2.auc - r1h.auc
        log.warning("  DELTA: risk ratio %+.3f, AUC %+.4f", d_rr, d_auc)
        if d_rr <= 0 and d_auc <= 0:
            log.warning("  THE REFINEMENT DID NOT HELP. Reported as a null, which §4.7 "
                        "explicitly requires. The loop's value is the traceable "
                        "decision record, not a guaranteed improvement.")
        report_only("refinement_delta_risk_ratio", round(d_rr, 4))
        report_only("refinement_delta_auc", round(d_auc, 4))

    # ---------------- the deliverable: the decision log ----------------------
    lines = ["# Refinement log — C6\n",
             "*Generated by `layer2_clinical/11_refinement_loop.py`. This log, not the "
             "performance change, is the deliverable for C6 (SPEC.md §4.7).*\n",
             "## Design\n",
             f"- Set B scoreable: {len(B)} patients, split once with seed {seed} into "
             f"{len(B_diag)} diagnosis / {len(B_hold)} final holdout.",
             "- The holdout was scored **exactly once**, after the refinement was fixed.",
             "- Refinement rules were written down before any diagnosis ran, so the "
             "choice of refinement could not be steered by the holdout.\n",
             "## Round 1 — predictions confronted with outcomes\n",
             f"Trained on set A (n={rounds[0].n_train}), predicted on the set B diagnosis "
             f"slice (n={rounds[0].n_eval}).\n",
             f"- Enrichment: {rounds[0].event_rate_positive:.1%} event rate in "
             f"predicted-positive vs {rounds[0].event_rate_negative:.1%} in "
             f"predicted-negative (RR {rounds[0].risk_ratio}, p={rounds[0].fisher_p})\n",
             "### Error taxonomy\n",
             pd.DataFrame(tax).to_markdown(index=False), "",
             "\n## Refinement decisions, each traced to its trigger\n",
             "| Rule | Fired | Observation that triggered it | Action |",
             "|---|---|---|---|"]
    for a in decisions:
        lines.append(f"| {a['rule']} | {'YES' if a['fired'] else 'no'} | "
                     f"{a['observation']} | {a.get('action', '—')} |")
    lines += ["", "\n## Round 2 — evaluated on the final holdout\n",
              "| Round | n | RR | AUC | event rate (pos) | event rate (neg) | p |",
              "|---|---|---|---|---|---|---|",
              f"| 1 (baseline) | {r1h.n_eval} | {r1h.risk_ratio} | {r1h.auc} | "
              f"{r1h.event_rate_positive:.1%} | {r1h.event_rate_negative:.1%} | {r1h.fisher_p} |",
              f"| 2 (refined) | {r2.n_eval} | {r2.risk_ratio} | {r2.auc} | "
              f"{r2.event_rate_positive:.1%} | {r2.event_rate_negative:.1%} | {r2.fisher_p} |",
              "",
              f"\n**Delta: risk ratio {d_rr:+.3f}, AUC {d_auc:+.4f}.**\n",
              (("The refinement did not improve the holdout on either metric. Reported "
                "as a null because §4.7 requires it, and because a loop that only ever "
                "reports improvements is not a loop, it is a selection procedure.")
               if (d_rr <= 0 and d_auc <= 0) else
               ("The refinement improved BOTH enrichment and discrimination."
                if (d_rr > 0 and d_auc > 0) else
                f"**Mixed, and the distinction matters.** Enrichment improved "
                f"(risk ratio {r1h.risk_ratio} -> {r2.risk_ratio}) while discrimination "
                f"did not (AUC {r1h.auc} -> {r2.auc}). Both refinements act on the "
                f"OPERATING POINT — a threshold move cannot change ranking at all, and "
                f"AUC is threshold-invariant — so this is the expected shape rather than "
                f"a contradiction. The deployable quantity got better; the model did not. "
                f"Reporting only the risk ratio would overstate what happened.")),
              "\n## What this loop cannot do\n",
              "- Set B has no endotype labels (D15), so the loop optimises **enrichment "
              "quality** against sepsis onset rather than agreement with a reference "
              "assignment. That is the deployable objective, but it is not the same "
              "objective the round-1 model was trained for.",
              "- One cycle only. A second cycle would need a second holdout."]
    (paths.RESULTS / "refinement_log.md").write_text("\n".join(lines))

    (paths.RESULTS / "L2_refinement_rounds.json").write_text(json.dumps(
        {"rules_prespecified": RULES, "decisions": decisions,
         "rounds": [asdict(r) for r in rounds],
         "delta_risk_ratio": round(d_rr, 4), "delta_auc": round(d_auc, 4)}, indent=2))
    log.info("wrote results/refinement_log.md — the C6 deliverable")
    return 0


if __name__ == "__main__":
    sys.exit(main())
