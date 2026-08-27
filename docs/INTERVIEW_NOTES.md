# Interview material

Companion to `SPEC.md` §11. Drafted as the work happens, so the answers are made of real
events rather than reconstructed at the end. Revise with real numbers once they exist.

---

## "How do you connect ML results to biology?" (C7)

**The strongest answer available in this project is not a result. It is a bug that a domain
fact caught.**

The plan was to apply a published three-variable classifier — IL-6, sTNFR1, bicarbonate —
to Olink Explore proteomics, since two of its three inputs are on the panel. That is a
category error. Olink reports NPX: a relative, log2, arbitrary-unit scale normalised within
the dataset. It tells you patient A has more IL-6 than patient B; it does not tell you how
many pg/mL, and the mapping to concentration differs per assay. The published equation's
intercept and `ln()` transform both assume absolute concentration.

The part worth telling is **why it nearly got through**. NPX values run about 0–15, and
`ln(concentration in pg/mL)` for real values runs about 5–10. They overlap. A plausibility
range check on the inputs passes NPX cleanly, and the model then emits numbers between 0
and 1 that look exactly like probabilities. Nothing errors. Prevalence looks reasonable.

The fix is domain knowledge, not statistics: circulating sTNFR1 is in the thousands of
pg/mL in every human, healthy or critically ill, so a cohort whose sTNFR1 maximum is below
100 is not measured in pg/mL. That single floor check is a reliable detector of a relative
scale, and no purely statistical guard would have found it.

It was caught by a test, not by inspection — the test asserted that NPX input must be
*refused*, and the first implementation failed it.

**What it demonstrates:** knowing what an assay actually measures, rather than treating
columns with matching names as interchangeable. Also the general habit — when a model
consumes fixed published coefficients, the units are part of the model, and a unit error is
invisible in the output.

**The size of the error, if asked:** one patient at IL-6 500 pg/mL and sTNFR1 10,000 pg/mL
differs by roughly 14 log-odds units between the pg/mL and ng/mL readings. Every
classification flips. The prevalence still looks plausible.

---

## "Tell me about a time you were wrong" / "how do you handle uncertainty in a source?"

**Four separate defects in the reference model, all found before any analysis ran, and all
invisible in output.** The equation came to me second-hand, as it does in most real work.

1. **It was labelled a probability but was a linear predictor.** Caught by noticing the
   right-hand side is unbounded and negative for every plausible patient. The logistic
   transform had been dropped in transcription. Several *published* downstream papers print
   it the same wrong way.
2. **The units were not stated.** Rather than pick one, I locked a calibration test so a
   unit error fails loudly. Later confirmed pg/mL from the appendix — and the test now
   anchors on the paper's own class-conditional medians (IL-6 116 vs 933 pg/mL; sTNFR1
   3,225 vs 7,452), which must land on the correct side of both published cut-points.
   They do: P = 0.034 and 0.672.
3. **It is not the model everyone thinks it is.** The paper's headline three-variable model
   uses IL-8, bicarbonate and protein C. The IL-6/sTNFR1/bicarbonate form is an *ancillary*
   model from supplementary Table S4. Calling it "the Sinha three-variable classifier" —
   which is how it is usually described — is a checkable error in front of readers who know
   this literature.
4. **One coefficient was wrong.** The circulating sTNFR1 coefficient is 1.1913; both
   primary appendices give 1.1903. The difference is ~0.009 log-odds, far too small to flip
   any classification, and therefore permanently invisible. It has propagated into at least
   one published paper. It is now pinned by a test.

**The point:** in a reproduction study the reference standard is the one thing that cannot
be allowed to be approximately right, because every downstream number is measured against
it. None of these four would have produced a wrong-looking result. Three cut-points also
circulate for nearby models — 0.295, 0.301, 0.276 — none of which belong to this one.

---

## "How do you validate a model?" (C4)

*To be written after the transfer experiment runs.* The distinction to lead with: whether
set A → set B measures genuine case-mix transfer or is confounded by differences in
measurement practice. Scoping quantifies the per-variable coverage gap between the two
hospital systems before the transfer result is interpreted, because a performance drop on B
means something completely different depending on the answer.

## "How do you handle missing data?" (C5)

*To be written after the ablation runs.* Mask channels rather than silent imputation,
because missingness in ICU data encodes ordering behaviour. Note also that cohort-median
imputation puts a population statistic back into the inference path, which is why
`shared/no_population_stats.py` rejects `SimpleImputer` alongside the scalers.

## "How would you improve a model in production?" (C6)

*To be written after the refinement loop runs.* See `results/refinement_log.md`.
