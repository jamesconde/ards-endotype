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

**Status: OPEN.** A consistency argument is not a citation. Three VERIFY items
(logit-vs-probability, units, threshold) close when the paper is read on day 1.

**Why this is logged rather than quietly handled:** a wrong-unit reference classifier
would silently miscalibrate every downstream comparison in both layers, and it would be
invisible — the outputs would still look like probabilities.

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

**Consequence for §5:** the cross-layer bridge is weaker than the spec assumes. It rests
on concordance of *structure* and *ordering*, not on the same classifier producing
comparable probabilities in both layers. State it as weaker.

**What would reverse it:** a published NPX → pg/mL conversion for these specific assays.
If one exists, use it and cite it. Do not invent one.

---

## D4 — MIMIC-IV stays out (day 0, standing)

SPEC.md §6.12 predicts drift toward MIMIC-shaped code, itemids and mimic-code concepts,
none of which exist on the open path. The MIMIC-IV demo is downloaded for schema
familiarity and smoke tests only. If MIMIC-IV appears anywhere outside §10, that is drift.
`grep -rn "itemid\|mimic_code\|d_items" layer2_clinical/` should stay empty.

---

## D5 — `config/genesets.yaml` deliberately not written on day 0

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
