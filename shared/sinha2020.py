"""The Sinha 2020 parsimonious 3-variable classifier — frozen external reference.

Source: Sinha, Delucchi, McAuley et al., Lancet Respir Med 2020;8:247-257,
"Development and validation of parsimonious algorithms to classify ARDS
phenotypes." Derived on ARMA/ALVEOLI/FACTT (n=2,022), validated on SAILS (n=715).

PRIOR_ART.md §4: implement exactly as printed. Do not refit, do not "improve"
the coefficients, do not restandardise the inputs. The design depends on this
model being external. This module therefore has no `fit` and never will.

---------------------------------------------------------------------------
THREE THINGS TO VERIFY AGAINST THE PAPER ON DAY 1 (Tier 1 reading, PMC7543720).
Each is flagged VERIFY in the code and repeated in results/L0_scoping.md.
---------------------------------------------------------------------------

VERIFY-1  RESOLVED (2026-08-27). The source line was printed as

              P(hyperinflammatory) = -11.9593 + 1.0138*ln(IL6+1)
                                     - 0.2436*HCO3 + 1.1913*ln(sTNFR1+1)

          but that right-hand side is a logistic model's *linear predictor*,
          not a probability — it is unbounded and negative for plausible
          inputs. The transcription dropped the transform. Secondary
          descriptions confirm the model is applied by computing predicted
          probabilities via logistic regression and assigning the
          hyperinflammatory class at a 0.5 probability cut-off. This module
          applies P = 1/(1+exp(-x)) accordingly.

VERIFY-2  Units. Not stated in the source that supplied the equation.
          This module requires IL-6 and sTNFR1 in **pg/mL** and HCO3 in
          **mmol/L**, on the following evidence rather than on the paper:

            archetypal hypo  (IL6 100 pg/mL, sTNFR1 5,000 pg/mL, HCO3 22)
                -> linear predictor -2.49, P = 0.08
            archetypal hyper (IL6 1500 pg/mL, sTNFR1 15,000 pg/mL, HCO3 18)
                -> linear predictor +2.53, P = 0.93

          Those straddle 0.5 sensibly. Putting sTNFR1 in ng/mL instead drives
          both archetypes to P < 0.01, which would classify essentially every
          patient hypoinflammatory and is clearly wrong. `self_test()` below
          locks this in, so a unit error fails loudly rather than producing a
          plausible-looking prevalence. Confirm the units in the paper anyway —
          a consistency argument is not a citation.

          Supporting, still secondary: the critical-care literature reports
          IL-6 in pg/mL throughout (day-1 ICU median ~200 pg/mL, IQR 43-1576;
          hyperinflammatory defined at >=100 pg/mL), which is consistent with
          the pg/mL reading and inconsistent with ng/mL.

VERIFY-3  The classification cut-point. 0.5 is used here and secondary
          descriptions of the model's application report a 0.5 probability
          cut-off. Confirm against the paper.

VERIFY-4  **WHICH VARIANT IS THIS?** Several three-variable parsimonious
          classifiers for these subphenotypes are in circulation, and they do
          not use the same variables:

            - IL-6 / sTNFR1 / bicarbonate   <- the equation implemented here
            - IL-8 / bicarbonate / protein C
            - IL-8 / bicarbonate / sTNFR1   (used in PHARLAP)

          The Sinha 2020 Lancet Respir Med citation may attach to a different
          variant than the coefficients above. Citing the wrong paper for the
          model you actually ran is precisely the sloppiness this project
          cannot afford, since its whole claim is careful reproduction.
          `MODEL_VARIANT` below names the implemented variant explicitly so
          the code cannot be ambiguous even while the citation is pending.

          Until resolved, `CITATION` is marked PROVISIONAL and the README must
          not name a specific paper for this equation.

---------------------------------------------------------------------------
WHY THIS MODULE MATTERS TWICE
---------------------------------------------------------------------------
It is the shared spine of SPEC.md §5: the same published classifier scores
Layer 1 in protein space (IL6, TNFRSF1A on Olink Explore 1536) and Layer 2 in
clinical-lab space. HCO3 is the variable both cohorts have.

But see `NPX_WARNING` below before applying it to Olink data. Olink NPX is not
pg/mL and the substitution is not free.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# --- the frozen equation -----------------------------------------------------
INTERCEPT = -11.9593
BETA_LN_IL6 = 1.0138        # ln(IL-6 pg/mL + 1)
BETA_HCO3 = -0.2436         # HCO3 mmol/L, untransformed
BETA_LN_STNFR1 = 1.1913     # ln(sTNFR1 pg/mL + 1)
THRESHOLD = 0.5             # VERIFY-3

# Name the variant in code, so that an unresolved citation cannot become an
# ambiguous implementation. See VERIFY-4.
MODEL_VARIANT = "IL6_sTNFR1_HCO3"
VARIANTS_IN_CIRCULATION = {
    "IL6_sTNFR1_HCO3": "implemented here",
    "IL8_HCO3_proteinC": "not implemented — reportedly the Lancet paper's 3-variable model",
    "IL8_HCO3_sTNFR1": "not implemented — reportedly used in PHARLAP",
}

CITATION_IS_PROVISIONAL = True      # VERIFY-4
CITATION = ("PROVISIONAL — the primary source for THESE coefficients is not yet confirmed. "
            "Commonly attributed to: Sinha P, Delucchi KL, McAuley DF, et al. Development "
            "and validation of parsimonious algorithms to classify ARDS phenotypes. "
            "Lancet Respir Med. 2020;8(3):247-257. doi:10.1016/S2213-2600(19)30369-8 — but "
            "that paper's three-variable model may use IL-8/bicarbonate/protein C rather "
            "than the IL-6/sTNFR1/bicarbonate form implemented here. Do not print this "
            "citation in the README until VERIFY-4 closes.")

# Plausible physiological ranges, used to catch unit errors at call time.
RANGES = {
    "il6_pg_ml": (0.0, 5.0e5),
    "stnfr1_pg_ml": (0.0, 2.0e5),
    "hco3_mmol_l": (2.0, 60.0),
}

# An upper bound alone does not catch the failure that actually matters. NPX
# values (~0-15) sit comfortably *inside* the pg/mL range and would sail
# through, producing numbers that look like probabilities and are not — which
# is precisely the silent failure DECISIONS.md D3 exists to prevent.
#
# sTNFR1 is the robust discriminator: circulating sTNFR1 is in the thousands of
# pg/mL in every human, healthy or critically ill. A cohort whose sTNFR1 maximum
# is below 100 is not measured in pg/mL. IL-6 cannot serve here — genuinely
# healthy controls do sit at a few pg/mL — so it warns rather than raises.
STNFR1_FLOOR_PG_ML = 100.0
IL6_SUSPICION_CEILING_PG_ML = 20.0

NPX_WARNING = """\
Olink NPX is a relative, log2, arbitrary-unit scale, normalised per-panel and
per-plate. It is NOT pg/mL. Substituting NPX into this equation silently
produces numbers that look like probabilities and are not — the intercept and
the ln() transform both assume an absolute concentration scale.

For the Layer 1 / Olink comparison (SPEC.md §3.7), use `rank_surrogate()`:
it applies the published coefficients' *signs and relative weights* to
within-cohort ranks, giving an ordering that can be cross-tabulated against
the unsupervised clusters. That is a weaker claim than "we applied the Sinha
classifier" and the README must say which one was done.

If a published NPX->pg/mL conversion exists for these specific Olink assays,
use it and cite it. Do not invent one.
"""


@dataclass(frozen=True)
class SinhaResult:
    linear_predictor: np.ndarray
    probability: np.ndarray
    label: np.ndarray            # "hyperinflammatory" / "hypoinflammatory"
    n_complete: int
    n_incomplete: int


def _check_units(name: str, x: np.ndarray) -> None:
    lo, hi = RANGES[name]
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return
    if finite.min() < lo or finite.max() > hi:
        raise ValueError(
            f"{name}: observed range [{finite.min():.4g}, {finite.max():.4g}] falls outside "
            f"the plausible range [{lo}, {hi}]. This is almost always a unit error. "
            f"See VERIFY-2 in shared/sinha2020.py."
        )


def _check_not_log_scale(il6: np.ndarray, stnfr1: np.ndarray) -> None:
    """Refuse NPX and other relative log-scale inputs. DECISIONS.md D3."""
    finite = stnfr1[np.isfinite(stnfr1)]
    if finite.size and finite.max() < STNFR1_FLOOR_PG_ML:
        raise ValueError(
            f"sTNFR1 maximum is {finite.max():.4g}, below the {STNFR1_FLOOR_PG_ML:g} pg/mL "
            f"floor that any human cohort clears. This is a unit error — the input looks "
            f"like a relative log-scale (Olink NPX), not pg/mL.\n\n{NPX_WARNING}"
        )
    il6_finite = il6[np.isfinite(il6)]
    if il6_finite.size > 2 and il6_finite.max() < IL6_SUSPICION_CEILING_PG_ML:
        log.warning(
            "IL-6 maximum is %.4g pg/mL across %d values. That is plausible only for a "
            "cohort of healthy controls; in acute illness it suggests a log or relative "
            "scale. Confirm the units before trusting these probabilities.",
            il6_finite.max(), il6_finite.size,
        )


def linear_predictor(
    il6_pg_ml: np.ndarray | pd.Series,
    hco3_mmol_l: np.ndarray | pd.Series,
    stnfr1_pg_ml: np.ndarray | pd.Series,
    check_units: bool = True,
) -> np.ndarray:
    """The published equation. Per-patient values only; nothing cohort-derived."""
    il6 = np.asarray(il6_pg_ml, dtype=float)
    hco3 = np.asarray(hco3_mmol_l, dtype=float)
    stnfr1 = np.asarray(stnfr1_pg_ml, dtype=float)

    if check_units:
        _check_units("il6_pg_ml", il6)
        _check_units("hco3_mmol_l", hco3)
        _check_units("stnfr1_pg_ml", stnfr1)
        _check_not_log_scale(il6, stnfr1)

    return (INTERCEPT
            + BETA_LN_IL6 * np.log(il6 + 1.0)
            + BETA_HCO3 * hco3
            + BETA_LN_STNFR1 * np.log(stnfr1 + 1.0))


def classify(
    il6_pg_ml, hco3_mmol_l, stnfr1_pg_ml,
    threshold: float = THRESHOLD,
    check_units: bool = True,
) -> SinhaResult:
    """Score patients. Rows with any missing input return NaN — never imputed.

    Imputing a missing input from a cohort median would put a population
    statistic back in the inference path (see shared/no_population_stats.py),
    and how often these three variables are *jointly* observed is itself one of
    the project's headline numbers (PRIOR_ART.md §5, on §1.5 scoping).
    """
    lp = linear_predictor(il6_pg_ml, hco3_mmol_l, stnfr1_pg_ml, check_units)
    p = 1.0 / (1.0 + np.exp(-lp))          # VERIFY-1
    complete = np.isfinite(lp)

    label = np.full(lp.shape, "unclassifiable", dtype=object)
    label[complete & (p >= threshold)] = "hyperinflammatory"
    label[complete & (p < threshold)] = "hypoinflammatory"

    return SinhaResult(
        linear_predictor=lp,
        probability=p,
        label=label,
        n_complete=int(complete.sum()),
        n_incomplete=int((~complete).sum()),
    )


def rank_surrogate(
    df: pd.DataFrame,
    il6_col: str,
    hco3_col: str,
    stnfr1_col: str,
) -> pd.Series:
    """Sign-and-weight-preserving ordering for scales that are not pg/mL (Olink NPX).

    NOT the Sinha classifier. Report it as "an ordering derived from the
    published coefficients", never as "the Sinha 2020 classification". It uses
    within-cohort ranks, so it is deliberately NOT portable and must never be
    fed to assert_no_population_stats as if it were an inference model.
    """
    log.warning(NPX_WARNING)
    r = df[[il6_col, hco3_col, stnfr1_col]].rank(pct=True)
    return (BETA_LN_IL6 * r[il6_col]
            + BETA_HCO3 * r[hco3_col]
            + BETA_LN_STNFR1 * r[stnfr1_col])


# Reported hyperinflammatory prevalence clusters around 20-35% across published
# cohorts. 0.30 is a stated assumption, not a measurement, and the concordance
# result must be reported as sensitive to it.
PUBLISHED_HYPER_PREVALENCE = 0.30


def rank_surrogate_classes(
    df: pd.DataFrame,
    il6_col: str,
    hco3_col: str,
    stnfr1_col: str,
    prevalence: float = PUBLISHED_HYPER_PREVALENCE,
) -> pd.Series:
    """Assign classes on a relative scale by splitting the rank score at `prevalence`.

    The honest way to use a fixed-coefficient absolute-concentration model on a
    relative scale (Olink NPX): keep the coefficients' directions and relative
    weights, discard the intercept — which is the part that encodes the units —
    and take the top `prevalence` fraction as hyperinflammatory.

    What this supports: class concordance and rank agreement against the
    unsupervised clusters (SPEC.md §3.7, §5).
    What it does NOT support: any statement of the form "P(hyperinflammatory) =
    x" for a patient, or a claim that the published classifier was applied.
    Label it "an ordering derived from the published coefficients, split at the
    published prevalence" — see DECISIONS.md D3.
    """
    if not 0.0 < prevalence < 1.0:
        raise ValueError(f"prevalence must be in (0,1), got {prevalence}")
    score = rank_surrogate(df, il6_col, hco3_col, stnfr1_col)
    cut = score.quantile(1.0 - prevalence)
    labels = pd.Series(np.where(score >= cut, "hyperinflammatory", "hypoinflammatory"),
                       index=df.index, dtype=object)
    labels[score.isna()] = "unclassifiable"
    log.info("rank surrogate: %d hyper / %d hypo at assumed prevalence %.0f%% "
             "(cut=%.4f). This is NOT the published classifier.",
             int((labels == "hyperinflammatory").sum()),
             int((labels == "hypoinflammatory").sum()), prevalence * 100, cut)
    return labels


def self_test() -> None:
    """Lock in VERIFY-2. Runs in CI and at import time in the test suite."""
    hypo = classify(il6_pg_ml=[100.0], hco3_mmol_l=[22.0], stnfr1_pg_ml=[5000.0])
    hyper = classify(il6_pg_ml=[1500.0], hco3_mmol_l=[18.0], stnfr1_pg_ml=[15000.0])

    assert hypo.label[0] == "hypoinflammatory", hypo
    assert hyper.label[0] == "hyperinflammatory", hyper
    assert hypo.probability[0] < 0.2 < 0.8 < hyper.probability[0], (
        f"archetypes no longer straddle the threshold: "
        f"hypo={hypo.probability[0]:.4f} hyper={hyper.probability[0]:.4f}"
    )

    # Monotonicity in each variable, in the published direction.
    up_il6 = classify([100.0, 5000.0], [22.0, 22.0], [5000.0, 5000.0]).probability
    assert up_il6[1] > up_il6[0], "P should rise with IL-6"
    up_hco3 = classify([100.0, 100.0], [18.0, 30.0], [5000.0, 5000.0]).probability
    assert up_hco3[1] < up_hco3[0], "P should fall as HCO3 rises"

    # Missing input propagates rather than imputing.
    miss = classify([np.nan], [22.0], [5000.0])
    assert miss.n_incomplete == 1 and miss.label[0] == "unclassifiable"

    log.info("✓ sinha2020 self-test passed (units and directions consistent)")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    self_test()
    print(f"\n{CITATION}\n")
    for name, (lo, hi) in RANGES.items():
        print(f"  {name:16s} plausible range [{lo}, {hi}]")
    print("\nVERIFY-1/2/3 in the module docstring are open until the paper is read.")
