# Decision log

Every deviation from `SPEC.md`, and every choice frozen before results were seen (I4).
Append-only. Each entry states the trigger, the decision, and what would reverse it.

---

## D1 — Layer 1 runs in Python, not R (day 0)

**Spec says:** `layer1_molecular/02_qc_normalize.R`, using GEOquery, limma and
ConsensusClusterPlus, on the grounds that "R is less painful here."

**Trigger:** there is no R on this machine (`which R` → nothing). Standing up
r-base plus a Bioconductor stack under WSL2 is hours of install, and it adds a second
language to `make all` and to `environment.yml` for the sake of three functions.

**Decision:** Python throughout.

| Spec (R) | Used instead | Note |
|---|---|---|
| `GEOquery::getGEO` | `GEOparse` | Parses the same series-matrix files, exposes the same `pData` characteristics fields |
| `ConsensusClusterPlus` | own implementation in `shared/consensus.py` | ~80 lines: resample 80%, cluster, accumulate co-assignment, compute PAC. The PAC criterion (§3.4) is explicit in the spec, so it must be computed rather than eyeballed either way |
| `GSVA::gsva(method="ssgsea")` | `gseapy.ssgsea` | Same Barbie et al. algorithm |
| `limma` | `scipy.stats` + Benjamini-Hochberg via `statsmodels` | Only needed for per-gene differential expression between discovery clusters (§3.5); no complex design matrices in this project |

**Cost, stated honestly:** the consensus-clustering implementation is mine, not the
canonical package, so it needs its own test against a synthetic dataset with known
structure (`tests/test_consensus.py`). That test is the price of this decision.

**What would reverse it:** if the PAC curve looks pathological and the implementation is
suspect, install R and cross-check one cohort against `ConsensusClusterPlus`. Half a day.

---

## D2 — `shared/sinha2020.py` infers its units rather than citing them (day 0)

**Trigger:** `PRIOR_ART.md` §1 prints the 3-variable equation but not the units of IL-6
and sTNFR1, and the equation is unusable without them — sTNFR1 in ng/mL instead of pg/mL
drives every plausible patient to P < 0.01.

**Decision:** implement with IL-6 and sTNFR1 in **pg/mL**, HCO3 in **mmol/L**, justified
by an internal consistency argument (archetypal hyper/hypo patients straddle P = 0.5
sensibly; see VERIFY-2 in the module). `self_test()` locks the calibration so a unit error
fails loudly instead of producing a plausible-looking prevalence.

**Status: CLOSED (2026-08-27)**, verified against the paper's supplementary appendix.

- **VERIFY-1 (logit vs probability): confirmed a transcription error in the source
  document, not a misreading.** Table S4 is a table of logistic-regression coefficients;
  the right-hand side is log-odds. Several *published* downstream papers print it as
  "Probability (hyperinflammatory) = …", and that labelling is wrong. `1/(1+exp(-x))` stands.
- **VERIFY-2 (units): CONFIRMED pg/mL for both cytokines, mmol/L for bicarbonate**, from
  appendix Table S2 row labels. The day-0 inference from internal consistency was correct.
  `self_test()` now anchors on the paper's own class-conditional medians (IL-6 116 vs 933
  pg/mL; sTNFR1 3,225 vs 7,452 pg/mL) rather than invented archetypes, which is a
  materially stronger check: the published medians must land on the correct side of both
  published cut-points, and they do (P = 0.034 and 0.672).
- **VERIFY-3 (threshold): two published cut-points**, both now exposed — 0.274 (Youden,
  sens 0.93 / spec 0.74) and 0.5 (sens 0.82 / spec 0.89). Two nearby values belong to
  *different* models and must not be used: 0.295/0.301 (the IL-8/protein C models) and
  0.276 (the IL-6/sTNFR1/vasopressor model used in HARP-2).
- **VERIFY-4: resolved.** See D7.
- **NEW — the coefficient we were given was wrong.** See D8.

**Why this is logged rather than quietly handled:** a wrong-unit reference classifier
would silently miscalibrate every downstream comparison in both layers, and it would be
invisible — the outputs would still look like probabilities. Concretely, one patient at
IL-6 500 pg/mL and sTNFR1 10,000 pg/mL differs by roughly 14 log-odds units between the
pg/mL and ng/mL readings. Every classification flips; the prevalence still looks plausible.

---

## D3 — Olink NPX cannot be fed to the Sinha equation (day 0)

**Trigger:** SPEC.md §3.7 proposes implementing the 3-variable classifier on Olink
Explore data because IL6 and TNFRSF1A are on the panel. But NPX is a relative, log2,
arbitrary-unit scale — not pg/mL. The equation's intercept and `ln()` transform both
assume absolute concentration.

**Decision:** `classify()` range-checks its inputs and refuses NPX. The Layer 1 comparison
uses `rank_surrogate()` instead — the published coefficients' signs and relative weights
applied to within-cohort ranks. The README must call this *"an ordering derived from the
published coefficients"*, never *"the Sinha 2020 classification"*.

**Confirmed on review (2026-08-27).** The spec's §3.7 plan is a category error, and the
reason the first version of the guard missed it is the part worth keeping: NPX values sit
around 0–15, and `ln(pg/mL)` for real concentrations sits around 5–10, so a range check
cannot separate them. The fix is domain-grounded rather than statistical — circulating
sTNFR1 is in the thousands of pg/mL in every human, so `max(sTNFR1) < 100` is a reliable
detector of a relative scale. A test caught this, not inspection.

**Resolution adopted:** `rank_surrogate_classes()` scores by rank, discards the intercept
(the part that encodes the units), and splits at the published hyperinflammatory prevalence
(assumed 0.30) to assign classes. Concordance against the unsupervised clusters is then
reportable. The assumed prevalence is a stated assumption, not a measurement, and the
concordance result must be reported as sensitive to it.

**Rejected alternative:** fitting our own logistic model on NPX with the same three
variables. Legitimate, but it is then *our* model, not a published external reference, and
the entire design depends on the reference being external. If it is ever done, it must be
labelled as ours.

**Consequence for §5:** the cross-layer bridge is weaker than the spec assumes. It rests
on concordance of *structure* and *ordering*, not on the same classifier producing
comparable probabilities in both layers. State it as weaker.

**What would reverse it:** a published NPX → pg/mL conversion for these specific assays.
If one exists, use it and cite it. Do not invent one.

---

## D7 — Which three-variable classifier is this? (day 0, RESOLVED same day)

**Resolution:** the attribution is correct, but **the location is the supplementary
appendix, Table S4 — not the main text**, and the model is an **ancillary** one rather
than the paper's headline three-variable model.

- The paper's headline models are **IL-8 / bicarbonate / protein C** (3-variable) and the
  same plus vasopressor use (4-variable). Figure 2's legend states this verbatim.
- Table S4 ("Model composition with variable coefficients of ancillary parsimonious
  models") contains a family of models permuting the six top variables. Our
  IL-6/sTNFR1/bicarbonate coefficients are one row of it, AUC 0.94 (0.92–0.95).
- Independently corroborated in the appendix of Sinha, Calfee, Cherian et al., *Lancet
  Respir Med* 2020;8(12):1209–1218 (COVID ARDS), which prints the same coefficients and
  describes the classifier as "comprised of the IL-6, sodium bicarbonate and sTNFR-1".

**Consequence for wording, and it is not cosmetic:** never write "the Sinha three-variable
classifier" about this code. Write "an ancillary parsimonious model from Sinha et al. 2020
(appendix Table S4)". A reader of this literature will know the headline 3-variable model
uses IL-8 and protein C, and calling ours by that name is a checkable error.

**Also implemented as a result:** the Table S4 ancillary model pairing sTNFR1 and
bicarbonate with **IL-8** (intercept −10.6110, IL-8 1.2902, HCO3 −0.2326, sTNFR1 1.0732,
AUC 0.95) is now available as `classify_il8()`. It has better specificity, and Olink
Explore carries CXCL8, so Layer 1 can use it. Layer 2 has neither cytokine.

**Carry into the writeup (verify before citing):** ATS 2026 abstract C24-04, Tripathi et
al., *AJRCCM* 2026;212(Suppl 1) — an improved parsimonious model derived on EARLI+FACTT
(n=1,604) and validated on ROSE (n=1,006) and CLOVERS (n=1,563). The augmented model is
IL-8 / angiopoietin-2 / sTNFR-1 / bicarbonate, AUROC 0.97 and 0.95. Its stated critique of
the model we implement, verbatim: *"This PCM may not fully capture the biology defined by
LCA, given the importance of IL-8 and endothelial markers in capturing the
Hyperinflammatory pathophysiology."* That the IL-6 model is in live prospective trial use
strengthens this project's relevance; the critique is a ready-made limitations paragraph.

---

## D11 — Proteomics cohort landscape, and three traps (day 1, partly OPEN)

Prompted by the `olink.com/mgh-covid-study` link in DOWNLOADS.md returning no accessible
page. A verification pass produced the following. Items marked **verified** were confirmed
by downloading and parsing the file; the rest are second-hand and marked as such.

### Traps that would have cost real time

1. **In Olink Explore 1536, TNFRSF1A sits on the NEUROLOGY panel (OID21155)** — not
   Inflammation and not Cardiometabolic. A reasonable-looking decision to subset to the
   Inflammation panel silently drops one of the two proteins the whole §3.7 comparator
   depends on. **Never subset Explore panels before checking which panel each target is on.**
2. **Assay names are not gene symbols.** The same protein appears as `TNFRSF1A`, `TNF-R1`
   (Olink Target 96, Gisby) and `TNF sR-I` (SomaScan). Match on **UniProt (P19438)** or
   gene ID, never on the assay label, or the join returns zero rows and looks like absence.
3. **Olink Target 96 CVD-II does not contain TNFRSF1A** — it is on CVD-III. Relevant only
   if a Target 96 cohort is used, but it is the kind of assumption that goes unchecked.

### Substitutes, if MGH access does not resolve

**CounterCOVID (Duijvelaar et al., Nat Commun 2024;15:744)** — **verified**, no
registration, direct download:
`https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-44986-w/MediaObjects/41467_2024_44986_MOESM6_ESM.xlsx`
(~29.7 MB). 332 patients, 535 samples, 203 with two timepoints, SomaScan ~7.3k. IL-6 and
TNFRSF1A both present. Crucially it carries **hard outcomes**: death, invasive ventilation,
day-28 and day-90 status. Caveat: it is the CounterCOVID imatinib RCT, so the randomised
arm is a covariate that must be adjusted for, not ignored.

**Gisby et al., eLife 2021** — Olink Target 96, downloadable from GitHub without a Dryad
bot wall. Only 55 cases, which is thin for clustering.

### The unit problem does NOT go away by switching cohorts

SomaScan reports RFU — relative fluorescence units — which is no more a concentration than
NPX is. D3 applies unchanged to CounterCOVID. Any cohort on a relative platform supports
rank-based concordance only, never evaluation of the published fixed-coefficient equation.
The only thing that would change this is an absolute-concentration immunoassay cohort.

### Still open

Where the MGH Explore 1536 data is actually downloadable today, and whether bicarbonate
exists anywhere in its clinical file (D9 says it is absent from the data dictionary in
hand). Both are being verified. Until they resolve, **Layer 1's proteomic arm has no
confirmed data source**, and Layer 2 — which needs neither — remains the critical path.

---

## D10 — FiO2 needs a physiological validity window (day 1)

**Trigger:** scoping set A found `FiO2` ranging from 0.000 to 10.000 on a fraction scale
where room air is 0.209 and pure oxygen is 1.0. 321 of 112,155 recorded values (0.29%) are
impossible: 81 exact zeros, ~238 below room air, and 2 above 1.0.

**Decision:** pre-registered validity window `[0.21, 1.0]` in `config/cohort_params.yaml`.
Hours with an out-of-range FiO2 are rejected as S/F candidates; **the patient is not
excluded**, only the impossible hour.

**Why it matters more than 0.29% suggests — the errors are not symmetric.** A zero sends
S/F to infinity. A value below room air inflates S/F and would falsely *exclude*. A value
above 1.0 *deflates* S/F and would falsely *admit* a patient as severely hypoxemic, adding
false positives to the cohort. Only two such values exist in set A, so the practical effect
here is negligible — the cohort count is 4,682 either way — but the asymmetry is the reason
to filter explicitly rather than rely on the count being small.

**Note against SPEC.md §6.2**, which argued this dataset avoids the MIMIC problem of "P/F
ratios of 0.4 and 3,000 in the same table." The same class of error exists here, at roughly
0.3% rather than MIMIC's rate. Much smaller, but not zero, and the spec's claim should be
softened rather than repeated.

---

## D9 — The MGH clinical file has no bicarbonate, which breaks the §5 bridge (day 0)

**Trigger:** the MGH cohort's clinical data dictionary (Life Science Alliance 2023
supplementary material 2, 40 variables) was obtained and read in full.

**What is there, and it is genuinely useful:**

- **WHO ordinal severity at days 0, 3, 7 and 28, plus WHO max.** A real longitudinal
  outcome for the trajectory work (C3), scoreable with a proportional-odds model.
- **A mortality outcome.** WHO score 1 is defined as "death within 28 days". Layer 1
  therefore has the outcome Layer 2 lacks (§6.1) — worth stating in the README, since it
  partially offsets the open-data cost.
- **A COVID 0/1 flag**, which makes the aetiology-versus-host-response comparison in §6.7
  buildable: do the endotype clusters mix COVID-positive and COVID-negative patients, or
  separate them? Mixing is evidence the endotype reflects host response, which is the
  premise a host-directed p38 modulator depends on.
- Comorbidity flags and presenting symptoms for case mix.

**What is NOT there, and it matters: bicarbonate.** SPEC.md §3.7 assumed "the third
[Calfee/Sinha variable] is in the clinical file if present; otherwise substitute an
available acidosis marker." It is not present, and **there is no acidosis marker to
substitute** — no bicarbonate, no lactate, no pH, no base excess. The labs are
neutrophils, lymphocytes, monocytes, creatinine, CRP, D-dimer, LDH, and a troponin flag.

**Consequence, and it is not small.** SPEC.md §5 states that the shared bicarbonate
variable "is what makes this more than an analogy." That variable does not exist on the
Layer 1 side. Combined with D3 (Olink NPX is not pg/mL, so the published equation cannot be
evaluated at all), the Layer 1 arm of the cross-layer comparator is blocked twice over:

1. the two proteins are on a relative scale the fixed coefficients cannot consume, and
2. the third variable is absent entirely.

**Decision:** the §5 bridge reduces to concordance of *structure* — both layers
independently recovering a two-class split with a consistent inflammatory direction. That
is a real but weaker claim, it must be written as weaker, and the paragraph in §5 asserting
the bicarbonate link has to go.

**Rejected:** dropping the bicarbonate term and running the remaining two coefficients.
Removing a term from a fixed-coefficient logistic model does not yield "the model minus one
variable" — the intercept was fitted in its presence, so the result is a different model
with no published provenance. If a two-protein ordering is used, it is ours and is labelled
ours.

**Second-order limitation:** every lab in the file is binned into five ordinal categories
rather than reported continuously, which caps the resolution of anything built on them.

**Data-quality flag:** the dictionary's own row for `ddimer_7_cat` describes it as "D-dimer
day 3 category" — an apparent copy-paste error in the published supplement. Verify against
the data before using the day-7 D-dimer variable.

---

## D8 — The sTNFR1 coefficient we were handed was wrong (day 0)

**Trigger:** the equation supplied in `PRIOR_ART.md` gives the sTNFR1 coefficient as
**1.1913**. Both primary Sinha-authored appendices give **1.1903**.

**Decision:** use 1.1903, pinned by a test.

**Why this is logged rather than silently corrected.** The difference is about 0.009
log-odds at a typical sTNFR1 value — far too small to change any classification, and
therefore permanently invisible in any output. It is exactly the class of error that
survives review indefinitely because nothing ever looks wrong. The typo has propagated into
at least one published paper (BMC Pulm Med 2025).

**Caveat on the evidence:** the two confirming sources share an author, so they are primary
but not fully independent. No third-party re-derivation was found.

---

## D7-ORIGINAL — the open question as first recorded (superseded above)

**Trigger:** at least three parsimonious three-variable classifiers for these subphenotypes
circulate, and they do not share variables:

| Variant | Status here |
|---|---|
| IL-6 / sTNFR1 / bicarbonate | the coefficients we were given, and what is implemented |
| IL-8 / bicarbonate / protein C | reportedly the Sinha 2020 *Lancet Respir Med* three-variable model |
| IL-8 / bicarbonate / sTNFR1 | reportedly used in PHARLAP |

**The risk:** citing Sinha 2020 *Lancet Respir Med* for coefficients that came from a
different paper. For a project whose entire claim is careful reproduction of published
work, attaching the wrong citation to the reference model is worse than any analysis error
— it is the exact failure the reader is checking for.

**Decision:** name the variant in code (`MODEL_VARIANT = "IL6_sTNFR1_HCO3"`), mark
`CITATION_IS_PROVISIONAL = True`, and **do not print any specific citation for this
equation in the README until the primary source is confirmed.** A test asserts the
provisional marker is present. The implementation is unambiguous even while the
attribution is not.

**Closes when:** the primary source for these specific coefficients is read.

**Also worth carrying into the writeup:** the IL-6/sTNFR1/bicarbonate model is reportedly
in live prospective trial use, and a 2026 ATS abstract reports an improved parsimonious
model (derived on EARLI+FACTT, validated on ROSE and CLOVERS) noting it may not fully
capture the biology LCA defines. Live prospective use strengthens this project's
relevance; the stated critique is a ready-made limitations paragraph. Verify both before
citing.

---

## D4 — MIMIC-IV stays out (day 0, standing)

SPEC.md §6.12 predicts drift toward MIMIC-shaped code, itemids and mimic-code concepts,
none of which exist on the open path. The MIMIC-IV demo is downloaded for schema
familiarity and smoke tests only. If MIMIC-IV appears anywhere outside §10, that is drift.
`grep -rn "itemid\|mimic_code\|d_items" layer2_clinical/` should stay empty.

---

## D5 — `config/genesets.yaml` deliberately not written on day 0 (CLOSED day 1)

**Closed 2026-08-27.** Written and committed in its own commit after scoping completed and
before any expression data was loaded or any clustering run. 17 genes in the p38/MK2 set,
3 mandatory comparators including an interferon negative control, 10 Olink assays for the
protein arm, and four pre-committed interpretation rules so the conclusion cannot be chosen
after seeing the numbers. The commit ordering is what makes the freeze auditable.

### Original entry

SPEC.md §3.6 requires the pathway gene sets frozen before §3.4 runs, and the day-1 line in
§7 says they are written *after* scoping and *before* any modelling. Writing them on day 0
would be harmless; writing them after seeing a clustering result would not be, and the
ordering is what makes the freeze meaningful. The file is created on day 1 and committed in
its own commit so the timestamp is auditable.

---

## D6 — pre-commit guards enforce three spec rules mechanically (day 0)

`tools/repo_guard.py` blocks commits that (a) contain row-level data files, (b) apply the
syndrome acronym to the unadjudicated primary cohort (§6.3), or (c) use discovery language
about established subphenotypes (§6.13, I1).

Guard (b) matches on co-occurrence within a line, so this file tripped it on first run —
while describing the guard. The rule says reword rather than widen the exemption list, so
the sentence above avoids the literal pairing. That is the intended behaviour of a guard
that errs toward false positives.

These are the three integrity rules most likely to erode under deadline pressure, and all
three are pattern-matchable. Guard (b) and (c) will produce occasional false positives;
the correct response is to reword, not to widen the exemption list.
