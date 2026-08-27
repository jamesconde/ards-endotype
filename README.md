# Endotyping for precision immunomodulation in acute hypoxemic respiratory failure

**Status: day 0 of 24. Scaffolding and guardrails only — no results yet.**
This README is a placeholder for the deliverable described in `SPEC.md` §9. Nothing below
is a finding.

---

## What this is, stated before anything else

The hyperinflammatory / hypoinflammatory subphenotypes of the acute respiratory distress
syndrome are **well established and independently replicated**, and this project
**reproduces them rather than discovering them**. They were identified by latent class
analysis of two randomised trials (Calfee et al., *Lancet Respir Med* 2014) and have since
been recovered across five or more cohorts, shown to carry differential mortality and
differential treatment response, and reduced to parsimonious classifiers usable at the
bedside (Sinha et al., *Lancet Respir Med* 2020; *AJRCCM* 2020). Full citations in
[`CITATIONS.md`](CITATIONS.md).

What is being tested here is not whether those subphenotypes exist. It is:

> **Do published subphenotype classifiers survive the conditions of routine, sparsely
> measured, unadjudicated ICU data — and how gracefully do they degrade as the variables
> they need go missing?**

Every published validation of these classifiers ran in curated research cohorts — ARDSNet
trials, EARLI, VALID, LUNG SAFE — where variables were measured on protocol and missingness
was low by design. Real screening does not look like that. That gap is what this project
occupies, and `PRIOR_ART.md` documents how the question was narrowed twice as the
literature closed off earlier framings.

## What the data is, and is not

The primary cohort is the PhysioNet/CinC Challenge 2019 corpus: 40,336 ICU patients, one
row per hour, drawn from two hospital systems. Two constraints are stated here rather than
buried in a limitations section:

- **It carries no mortality outcome.** The only label is sepsis onset. Mortality is carried
  separately by the eICU demo (~2,500 stays), with wide confidence intervals reported as
  wide. No composite "severity outcome" is constructed and presented as though it were
  mortality.
- **It is not an adjudicated ARDS cohort, and is never called one.** There is no imaging,
  no PEEP, and no ventilation flag on this data path, so Berlin criteria cannot be
  approximated even loosely. The cohort is defined as hypoxemic respiratory failure by an
  SpO₂/FiO₂ surrogate, and that is the term used in every figure, caption and variable
  name. A pre-commit guard enforces it.

## Repository

```
config/          pre-registered parameters, frozen before results were seen
shared/          assertions, portability checks, the external reference model
layer0_scoping/  day-1 descriptive scoping — no models
layer1_molecular/  transcriptomics + plasma proteomics
layer2_clinical/   ICU trajectories, classification, transfer, ablation, enrichment
tools/           pre-commit integrity guards
tests/           synthetic fixtures; the guards are themselves tested
```

Three pieces carry more weight than their size suggests:

- **`shared/no_population_stats.py`** — makes "this classifier is portable" an executable
  claim rather than an assertion. A deployed screening model sees one patient at a time and
  has no cohort to standardise against, so the check is functional: every patient's
  prediction must be identical whether computed alone or inside a batch. A fitted scaler,
  a cohort quantile, or median imputation all fail it.
- **`shared/sinha2020.py`** — the published 3-variable classifier, implemented from its
  printed equation and never refitted. It is a fixed external reference; the design depends
  on it being external. Three open verification items are recorded in the module rather
  than resolved by assumption.
- **`tools/repo_guard.py`** — blocks commits containing row-level data, applying the word
  "ARDS" to the unadjudicated cohort, or using discovery language about established
  findings.

## Reproducing

```bash
conda env create -f environment.yml && conda activate ards
# then follow DOWNLOADS.md — data acquisition is manual and deliberate
make check    # which sources are present
make test     # unit tests + integrity guards
make scope    # day 1 descriptive scoping
```

No data is committed. `data/` is gitignored before the first commit. The licences would
permit redistribution; a repo bloated with re-hosted public data is worse engineering
regardless, and the habit is the one worth having if restricted data is ever involved.

These are real patients, de-identified but not anonymous. There are no individual patient
walkthroughs and no per-patient trajectory plots anywhere in this repository.

## Documents

| File | What it holds |
|---|---|
| [`SPEC.md`](SPEC.md) | the full project specification (v1.4) |
| [`PRIOR_ART.md`](PRIOR_ART.md) | what is already known, and the two theses it closed off |
| [`DECISIONS.md`](DECISIONS.md) | every deviation from the spec, with its trigger |
| [`DOWNLOADS.md`](DOWNLOADS.md) | manual data acquisition |
| [`CITATIONS.md`](CITATIONS.md) | attribution — a licence condition for CC BY 4.0 sources |
