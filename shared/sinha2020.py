"""Sinha 2020 ancillary parsimonious classifiers — frozen external reference.

Source: Sinha, Delucchi, McAuley et al., Lancet Respir Med 2020;8(3):247-257,
"Development and validation of parsimonious algorithms to classify ARDS
phenotypes." Derived on ARMA/ALVEOLI/FACTT (n=2,022), validated on SAILS (n=715).
Coefficients implemented here come from the **supplementary appendix, Table S4**
("Model composition with variable coefficients of ancillary parsimonious
models"), not from the main text.

PRIOR_ART.md §4: implement exactly as published. Do not refit, do not "improve"
the coefficients, do not restandardise the inputs. The design depends on this
model being external. This module therefore has no `fit` and never will.

---------------------------------------------------------------------------
CALL IT BY ITS RIGHT NAME
---------------------------------------------------------------------------
The paper's **headline** three- and four-variable models are:

    3-variable: IL-8, bicarbonate, protein C
    4-variable: IL-8, bicarbonate, protein C, vasopressor use

Neither is what this module implements. The IL-6/sTNFR1/bicarbonate form is one
of a family of **ancillary** models in Table S4 that permute the six top
variables. The main text mentions them only in passing ("Replacing Protein C
with sTNFR1 resulted in similar AUCs"; "Replacing IL-8 with IL-6 generally
increased model sensitivity, but specificity was lower").

So never write "the Sinha three-variable classifier" about this code. Write
"an ancillary parsimonious model from Sinha et al. 2020 (appendix Table S4)".
The distinction is small, checkable, and exactly the kind of thing a reader of
this literature will check.

---------------------------------------------------------------------------
VERIFY-1..4 — ALL RESOLVED 2026-08-27 against the primary sources
---------------------------------------------------------------------------
VERIFY-1  It is a **linear predictor (log-odds)**, not a probability. Table S4
          is a table of logistic-regression coefficients and Figure 2 prints
          the models in full logistic form. Several downstream papers print
          "Probability (hyperinflammatory) = -11.9593 + ...", and that labelling
          is simply wrong. This module applies P = 1/(1+exp(-x)).

VERIFY-2  **Units confirmed: IL-6 pg/mL, sTNFR1 pg/mL, bicarbonate mmol/L.**
          Appendix Table S2 row labels read "Interleukin-6 (pg/mL)", "Soluble
          tumor-necrosis factor receptor-1 (pg/mL)", "Bicarbonate (mmol/L)".
          The day-0 inference from internal consistency was correct, and
          `self_test()` now uses the paper's own reference medians instead of
          invented archetypes.

VERIFY-3  **Two published cut-points for this specific model**, both exposed
          below: 0.274 (Youden) and 0.5. Table S4 reports sens/spec 0.93/0.74
          at Youden and 0.82/0.89 at p>=0.5. Do NOT use 0.295/0.301 — those
          belong to the IL-8/protein C models. Do not use 0.276 either; that is
          the ancillary IL-6/sTNFR1/**vasopressor** model used in HARP-2.

VERIFY-4  **The attribution is correct but the location is the appendix.**
          Corroborated independently in the appendix of Sinha, Calfee, Cherian
          et al., Lancet Respir Med 2020;8(12):1209-1218 (COVID ARDS), which
          prints the same coefficient set and states the classifier "is a
          logistic regression-based model comprised of the IL-6, sodium
          bicarbonate and sTNFR-1".

---------------------------------------------------------------------------
COEFFICIENT CORRECTION — the value we were handed was wrong
---------------------------------------------------------------------------
The sTNFR1 coefficient circulating in secondary sources is **1.1913**. Both
primary Sinha-authored appendices give **1.1903**. This module uses 1.1903.

The effect is small (~0.008 log-odds at typical sTNFR1) and would never have
been visible in any output — which is the entire reason it is written down
here rather than silently corrected. Note the two confirming sources share an
author, so they are primary but not fully independent.

---------------------------------------------------------------------------
WHY THIS MODULE MATTERS TWICE
---------------------------------------------------------------------------
It is the shared spine of SPEC.md §5: the same published model scores Layer 1
in protein space (IL6, TNFRSF1A on Olink Explore 1536) and Layer 2 in
clinical-lab space. HCO3 is the variable both cohorts have.

But see `NPX_WARNING` below before applying it to Olink data. Olink NPX is not
pg/mL and the substitution is not free — see DECISIONS.md D3.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# --- the frozen equations, from appendix Table S4 ----------------------------
# All coefficients apply to ln(x + 1) of the named analyte, per the Table S4
# footnote: "Coefficients were derived using the logarithmic (loge) transformed
# values... a value of 1 was added to the measured value prior to enable log
# transformation."
INTERCEPT = -11.9593
BETA_LN_IL6 = 1.0138        # ln(IL-6 pg/mL + 1)
BETA_HCO3 = -0.2436         # HCO3 mmol/L, untransformed
BETA_LN_STNFR1 = 1.1903     # ln(sTNFR1 pg/mL + 1)  — NOT 1.1913, see the docstring

THRESHOLD = 0.5             # sens 0.82 / spec 0.89 (Table S4)
YOUDEN_THRESHOLD = 0.274    # sens 0.93 / spec 0.74; confirmed in the COVID paper
AUC = "0.94 (0.92-0.95)"

MODEL_VARIANT = "IL6_sTNFR1_HCO3"

# The ancillary model pairing sTNFR1 and bicarbonate with IL-8 instead of IL-6.
# Better specificity, and Olink Explore carries CXCL8 — so this is available in
# Layer 1 even though Challenge 2019 has neither cytokine.
IL8_INTERCEPT = -10.6110
IL8_BETA_LN_IL8 = 1.2902
IL8_BETA_HCO3 = -0.2326
IL8_BETA_LN_STNFR1 = 1.0732
IL8_AUC = "0.95"

VARIANTS_IN_CIRCULATION = {
    "IL6_sTNFR1_HCO3": "implemented — Table S4 ancillary model, AUC 0.94",
    "IL8_sTNFR1_HCO3": "implemented as classify_il8() — Table S4 ancillary, AUC 0.95",
    "IL8_HCO3_proteinC": "NOT implemented — the paper's headline 3-variable model; "
                         "protein C is unavailable in every cohort here",
    "IL6_sTNFR1_vasopressor": "NOT implemented — ancillary model used in HARP-2; "
                              "note its Youden cut-point (0.276) is not ours",
}

CITATION_IS_PROVISIONAL = False     # closed 2026-08-27
CITATION = ("Sinha P, Delucchi KL, McAuley DF, et al. Development and validation of "
            "parsimonious algorithms to classify ARDS phenotypes. Lancet Respir Med. "
            "2020;8(3):247-257. doi:10.1016/S2213-2600(19)30369-8 — supplementary "
            "appendix Table S4, ancillary parsimonious models. NOTE: this is an "
            "ancillary model, NOT the paper's headline 3-variable model "
            "(IL-8/bicarbonate/protein C). Describe it accordingly.")

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


def classify_il8(
    il8_pg_ml, hco3_mmol_l, stnfr1_pg_ml,
    threshold: float = THRESHOLD,
    check_units: bool = True,
) -> SinhaResult:
    """The IL-8 / sTNFR1 / bicarbonate ancillary model (Table S4, AUC 0.95).

    Higher specificity than the IL-6 form. Available in Layer 1 because Olink
    Explore carries CXCL8; unavailable in Layer 2, which has neither cytokine.
    Its Youden cut-point is not the same as the IL-6 model's — this function
    therefore defaults to 0.5 and takes the threshold explicitly.
    """
    il8 = np.asarray(il8_pg_ml, dtype=float)
    hco3 = np.asarray(hco3_mmol_l, dtype=float)
    stnfr1 = np.asarray(stnfr1_pg_ml, dtype=float)

    if check_units:
        _check_units("hco3_mmol_l", hco3)
        _check_units("stnfr1_pg_ml", stnfr1)
        _check_not_log_scale(il8, stnfr1)

    lp = (IL8_INTERCEPT
          + IL8_BETA_LN_IL8 * np.log(il8 + 1.0)
          + IL8_BETA_HCO3 * hco3
          + IL8_BETA_LN_STNFR1 * np.log(stnfr1 + 1.0))
    p = 1.0 / (1.0 + np.exp(-lp))
    complete = np.isfinite(lp)
    label = np.full(lp.shape, "unclassifiable", dtype=object)
    label[complete & (p >= threshold)] = "hyperinflammatory"
    label[complete & (p < threshold)] = "hypoinflammatory"
    return SinhaResult(linear_predictor=lp, probability=p, label=label,
                       n_complete=int(complete.sum()),
                       n_incomplete=int((~complete).sum()))


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


# Class-conditional medians printed in appendix Table S2. Using the paper's own
# reference distributions rather than invented archetypes means the calibration
# check is anchored to published data, not to a plausibility argument.
REFERENCE_MEDIANS = {
    "hypoinflammatory": {"il6_pg_ml": 116.0, "stnfr1_pg_ml": 3225.0, "hco3_mmol_l": 24.0},
    "hyperinflammatory": {"il6_pg_ml": 933.0, "stnfr1_pg_ml": 7452.0, "hco3_mmol_l": 20.0},
}
# NOTE: the IL-6 and sTNFR1 medians are from Table S2. The bicarbonate values are
# representative, not quoted — Table S2's bicarbonate rows were not read directly.
# They only need to be directionally right for this check to bind.


def self_test() -> None:
    """Lock the units and directions. Runs in CI and in the test suite."""
    lo = REFERENCE_MEDIANS["hypoinflammatory"]
    hi = REFERENCE_MEDIANS["hyperinflammatory"]
    hypo = classify(il6_pg_ml=[lo["il6_pg_ml"]], hco3_mmol_l=[lo["hco3_mmol_l"]],
                    stnfr1_pg_ml=[lo["stnfr1_pg_ml"]])
    hyper = classify(il6_pg_ml=[hi["il6_pg_ml"]], hco3_mmol_l=[hi["hco3_mmol_l"]],
                     stnfr1_pg_ml=[hi["stnfr1_pg_ml"]])

    assert hypo.label[0] == "hypoinflammatory", hypo
    assert hyper.label[0] == "hyperinflammatory", hyper
    # The published class medians must land on the correct side of BOTH
    # published cut-points, or the units are wrong.
    assert hypo.probability[0] < YOUDEN_THRESHOLD, (
        f"published hypoinflammatory median scores P={hypo.probability[0]:.4f}, "
        f"above the Youden cut-point — check the units")
    assert hyper.probability[0] > THRESHOLD, (
        f"published hyperinflammatory median scores P={hyper.probability[0]:.4f}, "
        f"below 0.5 — check the units")

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
