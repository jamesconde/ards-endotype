#!/usr/bin/env python3
"""Prognostic enrichment simulation. SPEC.md §4.6, C8.

THE DISTINCTION THIS RESTS ON, stated first because a clinical-development
reader will probe it:

  PROGNOSTIC enrichment asks "which patients have more events". Directly
  observable, and that is what this does.

  PREDICTIVE enrichment asks "which patients benefit more from treatment". That
  needs a treatment-effect estimate and is confounded by indication in
  observational data. It is NOT done here, and §6.8 explains why. Blurring the
  two would be the single most damaging error available in this section.

So: this is a power calculation, not a causal claim.

What makes it worth doing rather than textbook is that every input is OBSERVED
in this project rather than assumed:

  - event rates per class come from the reference assignment
  - classifier accuracy comes from the portable classifier and the C5 ablation
  - **screening feasibility comes from the A->B transfer**: 82.4% of set B is
    scoreable on coverage-matched variables, 1.6% with the blood-gas block

That last one is the point. Enrichment requires *identifying* the enriched
population, and a patient who cannot be scored cannot be enrolled. Required N
therefore depends on where you screen, not only on how good the classifier is.

CAVEATS THAT TRAVEL WITH EVERY NUMBER HERE (§6.1, §6.3):
  - the endpoint is sepsis onset, NOT mortality, and not ARDS-adjudicated
  - the cohort is hypoxemic respiratory failure by oximetry surrogate
  - the simulation demonstrates the method; the numbers are illustrative

Usage
-----
    python layer2_clinical/09_enrichment_simulation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, step

log = setup("L2_enrichment")

ALPHA = 0.05
POWER = 0.80
RRR = 0.20          # assumed relative risk reduction — a design input, not a result


def n_per_arm(p_control: float, rrr: float = RRR,
              alpha: float = ALPHA, power: float = POWER) -> float:
    """Two-proportion sample size per arm, no continuity correction."""
    p1 = p_control
    p2 = p_control * (1 - rrr)
    if p1 <= 0 or p1 >= 1 or p1 == p2:
        return float("inf")
    pbar = (p1 + p2) / 2
    z_a = stats.norm.ppf(1 - alpha / 2)
    z_b = stats.norm.ppf(power)
    num = (z_a * np.sqrt(2 * pbar * (1 - pbar))
           + z_b * np.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2
    return float(num / (p1 - p2) ** 2)


def ppv_from_accuracy(sens: float, spec: float, prevalence: float) -> float:
    """Positive predictive value — what fraction of those enrolled truly are hyper."""
    tp = sens * prevalence
    fp = (1 - spec) * (1 - prevalence)
    return tp / (tp + fp) if (tp + fp) > 0 else float("nan")


def main() -> int:
    con = duckdb.connect(str(paths.DUCKDB))
    ref = con.execute(
        "SELECT patient_id, is_hyper_like_k2 FROM reference_assignment").df()
    con.close()
    static = pd.read_parquet(paths.SCRATCH / "features" / "static.parquet")
    df = static.merge(ref, on="patient_id", how="inner")

    prevalence = float(df.is_hyper_like_k2.mean())
    p_hyper = float(df.loc[df.is_hyper_like_k2 == 1, "sepsis_any"].mean())
    p_hypo = float(df.loc[df.is_hyper_like_k2 == 0, "sepsis_any"].mean())
    p_all = float(df.sepsis_any.mean())
    log.info("observed: prevalence(hyper)=%.3f  event rate hyper=%.3f  hypo=%.3f  "
             "overall=%.3f", prevalence, p_hyper, p_hypo, p_all)
    if not (p_hyper > p_hypo):
        raise CheckpointError(
            f"hyper class does not have the higher event rate ({p_hyper:.3f} vs "
            f"{p_hypo:.3f}). Prognostic enrichment on this class would be backwards.")
    report_only("enrichment_event_rate_hyper", round(p_hyper, 4))
    report_only("enrichment_event_rate_hypo", round(p_hypo, 4))

    with step(log, "unstratified vs perfectly-enriched"):
        n_unstrat = n_per_arm(p_all)
        n_perfect = n_per_arm(p_hyper)
        log.info("  unstratified enrollment      : %.0f per arm (event rate %.1f%%)",
                 n_unstrat, 100 * p_all)
        log.info("  perfect enrichment on hyper  : %.0f per arm (event rate %.1f%%)",
                 n_perfect, 100 * p_hyper)
        log.info("  reduction: %.1f%% fewer patients per arm",
                 100 * (1 - n_perfect / n_unstrat))

    # --- required N as a function of classifier accuracy --------------------
    # Misclassification dilutes the enrolled population toward the base rate,
    # which is exactly why classifier degradation costs trial size.
    rows = []
    with step(log, "required N vs classifier accuracy"):
        for acc in [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50]:
            sens = spec = acc          # symmetric, stated as a simplification
            ppv = ppv_from_accuracy(sens, spec, prevalence)
            p_enrolled = ppv * p_hyper + (1 - ppv) * p_hypo
            n = n_per_arm(p_enrolled)
            # Screening burden: how many patients must be screened per enrollee?
            frac_positive = sens * prevalence + (1 - spec) * (1 - prevalence)
            screen_per_enrolled = 1 / frac_positive if frac_positive > 0 else np.inf
            rows.append({
                "classifier_accuracy": acc,
                "ppv": round(ppv, 4),
                "event_rate_enrolled": round(p_enrolled, 4),
                "n_per_arm": round(n),
                "n_total": round(2 * n),
                "screen_per_enrolled": round(screen_per_enrolled, 2),
                "screened_total": round(2 * n * screen_per_enrolled),
                "vs_unstratified": round(n / n_unstrat, 3),
            })
            log.info("  acc=%.2f  PPV=%.2f  enrolled event rate=%.1f%%  "
                     "N/arm=%5.0f  screen/enrolled=%.1f", acc, ppv,
                     100 * p_enrolled, n, screen_per_enrolled)
    curve = pd.DataFrame(rows)
    curve.to_csv(paths.RESULTS / "L2_enrichment_curve.csv", index=False)

    break_even = curve[curve.n_per_arm >= n_unstrat]
    be = float(break_even.classifier_accuracy.max()) if len(break_even) else None
    if be is not None:
        log.warning("BREAK-EVEN: at classifier accuracy %.2f, enrichment stops "
                    "helping — required N matches unstratified enrollment.", be)
    report_only("enrichment_break_even_accuracy", be if be is not None else -1)

    # --- the site-dependent part, which is this project's contribution -------
    with step(log, "screening feasibility by site"):
        try:
            tr = pd.read_csv(paths.RESULTS / "L2_transfer_results.csv").set_index("feature_set")
            scoreable = {"restricted_setB": float(tr.loc["restricted", "frac_scoreable_B"]),
                         "naive_setB": float(tr.loc["naive", "frac_scoreable_B"]),
                         "restricted_setA": float(tr.loc["restricted", "frac_scoreable_A"]),
                         "naive_setA": float(tr.loc["naive", "frac_scoreable_A"])}
        except Exception as e:  # noqa: BLE001
            raise CheckpointError(f"need L2_transfer_results.csv first: {e}")

        acc_assumed = 0.85     # roughly the portable classifier's observed agreement
        ppv = ppv_from_accuracy(acc_assumed, acc_assumed, prevalence)
        p_enrolled = ppv * p_hyper + (1 - ppv) * p_hypo
        n = n_per_arm(p_enrolled)
        frac_positive = acc_assumed * prevalence + (1 - acc_assumed) * (1 - prevalence)

        site_rows = []
        for name, frac_scoreable in scoreable.items():
            # A patient who cannot be scored cannot be enrolled.
            eff = frac_positive * frac_scoreable
            screened = 2 * n / eff if eff > 0 else np.inf
            site_rows.append({
                "scenario": name,
                "frac_scoreable": round(frac_scoreable, 4),
                "n_total_enrolled": round(2 * n),
                "patients_to_screen": round(screened),
                "screen_multiplier_vs_full_coverage": round(1 / frac_scoreable, 2)
                if frac_scoreable else None,
            })
            log.warning("  %-18s scoreable %5.1f%%  ->  screen %8.0f to enrol %d",
                        name, 100 * frac_scoreable, screened, round(2 * n))
        pd.DataFrame(site_rows).to_csv(
            paths.RESULTS / "L2_enrichment_by_site.csv", index=False)

    out = {
        "design": {"alpha": ALPHA, "power": POWER, "assumed_RRR": RRR,
                   "endpoint": "sepsis onset (NOT mortality, NOT ARDS-adjudicated)"},
        "observed": {"prevalence_hyper": round(prevalence, 4),
                     "event_rate_hyper": round(p_hyper, 4),
                     "event_rate_hypo": round(p_hypo, 4),
                     "event_rate_overall": round(p_all, 4)},
        "n_per_arm_unstratified": round(n_unstrat),
        "n_per_arm_perfect_enrichment": round(n_perfect),
        "break_even_accuracy": be,
        "curve": rows,
        "by_site": site_rows,
        "what_this_is_not": ("Prognostic enrichment only. No treatment-effect "
                             "estimate is made and no predictive-enrichment claim "
                             "is implied — see §6.8."),
        "caveats": ("Event rates come from a non-adjudicated hypoxemic respiratory "
                    "failure cohort with sepsis onset as the endpoint. The method is "
                    "the deliverable; the numbers are illustrative."),
    }
    (paths.RESULTS / "L2_enrichment_summary.json").write_text(json.dumps(out, indent=2))
    log.info("wrote results/L2_enrichment_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
