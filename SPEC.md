# ARDS Endotyping for Precision Immunomodulation — Project Spec (Layers 1 & 2)

**Target:** GEn1E Lifesciences — ML Engineer, Computational Biology (Palo Alto). Apply by 2026-09-19.
**Author:** James
**Spec version:** 1.4 — written to be handed to Claude Code as the driving document.

> **Changes in 1.4 — the purpose changed, so the design did.** v1.0–1.3 optimized for a novel scientific claim, and narrowed twice as prior art closed off each thesis. **That was the wrong objective.** This project's job is to be the evidence base for a résumé line and an interview, mapped to a specific job description — so it optimizes for **demonstrated competency coverage**, with scientific integrity as a constraint rather than novelty as a goal. Breadth is restored (v1.3's narrowing shrank the skill surface). New **§0.4 competency matrix** is the organizing device. §1 success criteria rewritten as competency demonstrations. Prognostic **enrichment simulation** added (§4.6) — the drug-development bullet, and safe because it is a power calculation, not a causal claim. **Iterative refinement loop** added (§4.7) to hit a JD bullet nothing else addresses. New **§11 interview material**. Missingness robustness retained as one component, not the whole thesis.
>
> **Changes in 1.3 — the question changed.** Portability spine; day-1 scoping; p38 demoted. *Superseded in part: see `prior-art-and-reading-list.md` §1 — the portability thesis was largely solved by Sinha et al. 2020. Retained as a validation step, not a claim.*
>
> **Changes in 1.2 — open data only.** Layer 2 moved off credentialed MIMIC-IV/eICU onto **PhysioNet/CinC Challenge 2019** (Open Access, CC BY 4.0, 40,336 ICU patients, 42 MB) plus the open eICU and MIMIC demos. §4 rewritten end to end. Credentialing risk eliminated (old §6.1). New §6.1 covers the real cost: **Challenge 2019 has no mortality outcome.** Cohort weakened from a Berlin proxy to hypoxemic respiratory failure by S/F ratio (§4.3, §6.3). External validation is now nearly free via the dataset's built-in two-hospital split (§4.5d) and is no longer on the cut list. Storage now trivial (§0.3, §6.10). Timeline gains ~4 days (§7).
>
> **Changes in 1.1:** storage budget corrected downward. Table skip-list and filtered-ingest pattern. BigQuery as a zero-download alternative. Downgraded the SQL-dialect porting risk given `mimic-code`'s official DuckDB build. *Most of 1.1 applied to MIMIC-IV and is superseded by 1.2; it is retained only for the optional stretch in §10.*

---

## 0. Read this first

### 0.1 What this project is

GEn1E's lead asset (GEn-1124 IV) is a p38α:MK2 dual signal modulator in Phase 2 for ARDS. Their GRID platform page names four applications this project deliberately mirrors: **AI-powered endotyping for precision medicine**, **personalization & trajectory mapping**, **directed drug-disease mapping**, and **scoring stratification**. It also names **hidden Markov models** as one of three model families they use.

### What this project is for — read before anything else

**This is a portfolio artifact, not a paper.** Its purpose is to make a specific person credible for a specific job: ML Engineer, Computational Biology at a clinical-stage company developing a p38α:MK2 dual signal modulator for ARDS. Success is measured by whether it produces (a) a résumé line that survives scrutiny, and (b) four or five interview answers with real substance behind them.

That reframing, introduced in v1.4, resolves the problem that sank v1.1–v1.3. Those versions chased a novel claim and narrowed each time prior art closed one off. **Novelty is not the objective.** The biology here is established, and the project should say so plainly and cite it — see `prior-art-and-reading-list.md`. What has to be original is the *execution*: the pipeline, the validation discipline, and the judgment on display.

The working question, chosen because it naturally requires the widest span of relevant skills rather than because it is unanswered:

> **Build and validate an end-to-end endotyping pipeline for acute hypoxemic respiratory failure — unsupervised discovery, supervised classification, longitudinal trajectory modelling, molecular characterisation, and trial-enrichment simulation — and stress-test it under the data conditions a real screening deployment would face.**

Framed for a reader in one sentence: *the subphenotypes are known; this asks whether you could actually operationalise them to enrich a trial using data a real ICU captures.*

**Honesty constraint, non-negotiable.** Established findings are cited as established. The README states in its opening paragraph that hyperinflammatory/hypoinflammatory subphenotypes are well-replicated and that this project reproduces rather than discovers them. Overclaiming is the one failure mode that converts a good artifact into a liability — the reader knows this literature.

**Which is itself a selling point.** A candidate who reproduces known work accurately, cites it properly, states limitations before being asked, and stress-tests the result is demonstrating exactly the judgment a small research team needs. A candidate claiming novel discovery from public data is demonstrating the opposite.

**Data-access principle (v1.2):** every dataset in this spec is open access — no credentialing application, no CITI course, no institutional affiliation, no signed DUA, no IRB. This is a deliberate constraint. It costs real capability (§6.1) and buys schedule certainty against a hard September 19 deadline.

Two layers are in scope here:

- **Layer 1 — Molecular endotypes.** Public transcriptomics + plasma proteomics. Biological anchor, biomarker work, p38 arm.
- **Layer 2 — Clinical endotyping and deployment stress-test.** Challenge 2019 hourly ICU data. Discovery, classification, trajectory, missingness robustness, cross-hospital transfer, enrichment simulation, iterative refinement.

Layer 3 (LINCS drug-disease reverser mapping) and Layer 4 (heterogeneous treatment effect / enrichment design) are **out of scope for this spec** and are written up as "next steps" in the final report. Do not start them. See §7.

### 0.2 Non-negotiables

1. **Every dataset is open access and de-identified; cite each one properly.** Challenge 2019 is CC BY 4.0 — attribution is a licence condition, not a courtesy. Put the Reyna et al. citation and the PhysioNet platform citation in the README, and cite the source publication for every GEO series and the Olink cohort.
2. **Still do not commit data.** `.gitignore` excludes `data/`, `*.duckdb`, `*.psv`, `*.csv`, `*.parquet` before the first commit. The licence permits redistribution here, but a repo bloated with re-hosted public data is worse engineering, not better, and the habit is the one you want if you ever touch restricted data. Publish code + aggregate figures + summary statistics; link to the source for the data.
3. **These are still real patients.** De-identified is not anonymous, and the tone of the writeup should reflect that. No "interesting cases," no individual patient walkthroughs, no per-patient trajectory plots with identifiers of any kind.
4. **Every extraction step gets an assertion checkpoint** with an expected row-count range (§8.3). A SQL query that silently returns 0 rows is the single most common way this kind of project dies quietly.

### 0.3 Machine

Acer Predator PTN16-51, Intel Core Ultra 7 155H, NVIDIA GPU, 16GB RAM, Windows 11. Use WSL2 for the whole pipeline — mimic-code SQL and the bioconductor stack are both materially less painful there.

**Disk: under 5GB total. This is no longer a constraint.**

| Item | Size |
|---|---|
| PhysioNet/CinC Challenge 2019 | 42 MB |
| eICU-CRD Demo | ~50 MB |
| MIMIC-IV Clinical Database Demo | ~100 MB |
| GEO series matrices (3 cohorts) | ~1–2 GB |
| MGH Olink NPX + clinical | ~200 MB |
| DuckDB + intermediates | ~2 GB |

The v1.0 figure of 120GB and the v1.1 figure of 40GB both belonged to the MIMIC-IV path and no longer apply. The whole open-data project fits comfortably on a laptop with room to spare, and every dataset is small enough to hold in RAM.

---

### 0.4 Competency coverage matrix — the organizing device

Every component below exists because it evidences a line from the job description. **If a component maps to nothing in this table, cut it.** If a JD line has no component, that is a gap worth closing.

| Job description requirement | Component that evidences it | Where |
|---|---|---|
| Strong Python: pandas, NumPy, scikit-learn, scientific computing | The entire pipeline | throughout |
| Clustering, unsupervised learning | Consensus clustering, GMM/latent class | §3.4, §4.5a |
| Classification, supervised learning | Endotype classifiers, GBM + logistic + tree | §4.5a2 |
| Regression | Time-to-event, trajectory, proteomic scores | §4.5c, §3.7 |
| Feature engineering | Binned longitudinal tensors, mask channels, ratio features | §4.4 |
| Model evaluation, model selection, validation, generalization | Nested CV, BIC/PAC selection, calibration, cross-hospital transfer | §3.4, §4.5d |
| Understanding of statistical modeling, interpreting model outputs | Effect sizes with CIs, SHAP/coefficient interpretation, comparator tests | §3.6, §4.5 |
| **Pipelines where predictions, outcomes, and new data feed back into model refinement** | **Iterative refinement loop** | **§4.7** |
| Clean, transform, merge, analyze large complex datasets | 6 sources, 3 modalities, 4 file formats | §3.1, §4.1 |
| Git/GitHub, reproducible research workflows | `make all`, pinned env, frozen configs, assertion checkpoints | §8 |
| EHR and clinical datasets | Challenge 2019, eICU demo, MIMIC demo | §4.1 |
| ICU / critical-care data | Same — 40,336 ICU patients | §4.1 |
| Longitudinal patient-level data | Hourly ICU trajectories; Olink day 0/3/7 | §4.4, §3.7 |
| Laboratory and physiological measurements | The entire Layer 2 feature set | §4.4 |
| Clinical outcomes and treatment/exposure data | Sepsis onset, eICU mortality, acuity trajectory | §4.4 |
| Familiarity with biology; connect computational results to biological hypotheses | p38/MK2 mechanism, pathway scoring, assay-limits reasoning | §3.6 |
| Read and understand biomedical literature | Prior art review; reproducing a published classifier from its equation | companion doc |
| Drug discovery / ML applied to biomarker data for drug development | Enrichment simulation; p38 target context | §4.6 |
| Independent research experience | The whole thing, self-scoped and self-critiqued | — |

**Deliberately not covered, and why:** predictive enrichment / heterogeneous treatment effect (§6.8 — confounded by indication on observational data), and drug-disease signature mapping (Layer 3, out of scope). Both go in "next steps" with the reason. *Naming what you chose not to do, and why, is itself an interview answer.*

---

## 1. Success criteria

Two independent tests. **Both must pass.** C-criteria are what gets you hired; I-criteria are what stops the artifact backfiring.

### Competency demonstrations (C) — the deliverable

| # | Demonstrates | Passes when |
|---|---|---|
| C1 | Unsupervised learning on real clinical data | Stable clusters recovered in ≥2 independent cohorts; k selected by a stated criterion, not by eye |
| C2 | Supervised classification + external model reuse | A published classifier (Sinha 3-variable) implemented from its equation and applied unmodified; own classifiers benchmarked against it |
| C3 | Longitudinal / trajectory modelling | HMM over 72h ICU trajectories, and proteomic trajectory over day 0/3/7 |
| C4 | Generalization and validation discipline | Frozen-model transfer set A → set B, plus → eICU. Calibration reported, not just discrimination |
| C5 | Deployment realism | Degradation curve under progressive variable ablation; native-sparse vs native-complete comparison |
| C6 | **Feedback-loop pipeline design** | §4.7 loop implemented and its effect measured — the JD bullet nothing else covers |
| C7 | Biology → computation → biology | p38/MK2 pathway reasoning including *why the assay can't settle it*; molecular endotype characterisation |
| C8 | Drug-development relevance | Prognostic enrichment simulation (§4.6) with required-N under stratified vs unstratified enrollment |
| C9 | Multi-source data engineering | ≥5 sources, ≥3 modalities, merged and reconciled with a documented attrition trail |
| C10 | Reproducibility | `make all` regenerates every figure from raw data on a clean machine |

### Integrity constraints (I) — non-negotiable

| # | Constraint |
|---|---|
| I1 | Every established finding is cited as established. The README's opening paragraph states this reproduces known subphenotypes |
| I2 | No claim the data can't support. Specifically: no mortality claim from Challenge 2019 (§6.1); no "ARDS" label on that cohort (§6.3); no causal treatment-effect claim (§6.8) |
| I3 | Limitations stated before a reader finds them, in the section where they apply, not only in a limitations block |
| I4 | Pre-registered analysis choices frozen before results are seen (§1.5) |
| I5 | Negative and null results reported as findings. **S5's null is expected (§3.6) and reporting it is a pass, not a failure** |

**On C-criteria vs novelty:** none of C1–C10 requires an unpublished finding. That is the point. A reviewer assessing whether you can do this work is reading for whether the clustering was selected properly, whether the transfer was frozen, whether calibration was checked, and whether you knew what your assay could not show. Those are all visible in reproduced work.

---

## 1.5 Day 1 — descriptive scoping, before anything is frozen

**Run this before writing the analysis plan. Do not skip it, and do not let it become the analysis.**

You do not yet know whether these datasets support the questions above. Committing to an analysis plan built on assumed data properties is how projects discover on day 15 that the cohort doesn't exist. One day, purely descriptive, no models.

Produce `results/L0_scoping.md` answering:

**Challenge 2019**
- Confirm 40,336 patients; confirm the A/B split and that both parse
- Per-column coverage: what fraction of patients have ≥1 non-null value? Per-hour density?
- **`FiO2` coverage specifically** — the cohort definition depends on it (§4.3). If coverage is very low, §4.3 needs rewriting before it is written
- `HCO3`, `Lactate`, `Creatinine`, `Bilirubin_total`, `Platelets`, `WBC` coverage — these carry P1
- Distribution of `ICULOS`; how many have ≥24h of record?
- Do sets A and B differ in coverage or case mix? **This determines whether P2 is a transfer test or a missingness test** — a crucial distinction

**GEO cohorts**
- `pData` column names and unique values for GSE65682 — **are the MARS endotype labels actually present?** (§3.3 fallback hinges on this)
- Genes surviving probe collapse per cohort; platform overlap across the three

**MGH Olink**
- Exact file schemas; confirm `IL6` and `TNFRSF1A` are on the panel (they carry §3.7's external comparator)
- Fraction below LOD per protein at day 0; how many proteins survive the §3.2 filter
- How many subjects have all three timepoints?

**Then, and only then:** pick the final question wording, write `config/genesets.yaml`, commit both, and start.

### The discipline this protects

Exploration is for learning what the data *is*. It is not for choosing which result to pursue. Clustering four ways and keeping the separation you like is p-hacking with extra steps, and reviewers of this literature are trained to spot it. The rule: **descriptive statistics before the freeze, inferential claims after, never the reverse.** If scoping reveals the plan is infeasible, change the plan explicitly and write down that you did — that is honest, and different from quietly steering toward a nicer answer.

---

## 2. Repo layout

```
ards-endotypes/
├── README.md                      # public-facing; figures + narrative
├── SPEC.md                        # this file
├── environment.yml                # conda
├── .gitignore                     # data/, *.duckdb, *.csv, *.parquet
├── .pre-commit-config.yaml        # PHI guard hook
├── config/
│   ├── genesets.yaml              # PRE-REGISTERED pathway gene sets (§3.6)
│   └── cohort_params.yaml         # Berlin thresholds, time windows
├── layer1_molecular/
│   ├── 01_download.py
│   ├── 02_qc_normalize.R          # GEOquery + limma; R is less painful here
│   ├── 03_anchor_mars.py          # supervised anchor on GSE65682
│   ├── 04_cluster_discovery.py
│   ├── 05_transfer_score.py       # fixed signature -> held-out cohorts
│   ├── 06_olink_longitudinal.py
│   └── 07_pathway_scores.py
├── layer2_clinical/
│   ├── 01_build_duckdb.py         # load Challenge 2019 .psv -> DuckDB
│   ├── 02_missingness_audit.py    # FiO2 coverage; selection check (§6.2)
│   ├── 03_cohort_hypoxemic.sql
│   ├── 04_features_static.sql
│   ├── 05_features_longitudinal.sql
│   ├── 06_reference_assignment.py # (a) cohort-standardised LCA — the baseline
│   ├── 07_portable_classifier.py  # (a2) P1 — no population stats at inference
│   ├── 08_parsimony_curve.py      # (a3) P3
│   ├── 09_transfer_setA_to_setB.py# (d1) P2 — the headline transfer test
│   ├── 10_eicu_demo_mortality.py  # (d2)
│   ├── 11_hmm_trajectory.py       # C3
│   ├── 12_multitask_model.py      # supporting
│   ├── 13_missingness_ablation.py # C5 — degradation curve
│   ├── 14_enrichment_simulation.py# C8 — required-N under stratification
│   └── 15_refinement_loop.py      # C6 — RefinementRound, one full cycle
├── layer0_scoping/
│   └── 01_scope_all_sources.py    # §1.5 — descriptive only, runs day 1
├── shared/
│   ├── assertions.py              # checkpoint helpers (§8.3)
│   ├── no_population_stats.py     # §4.5a2 portability assertion
│   └── plotting.py
├── figures/                       # committed; aggregate only
├── results/                       # committed; summary stats only, no row-level
│   └── refinement_log.md          # C6 deliverable — decisions + triggers
└── logs/
```

---

## 3. LAYER 1 — Molecular endotypes

### 3.1 Datasets

| ID | Cohort | Platform | n | Role | Access |
|---|---|---|---|---|---|
| GSE65682 | MARS (Amsterdam ICU) | Affy HG-U219 (GPL13667) | 802 samples; 760 ICU patients; ~479 with 28-day outcome + published MARS endotype label | **Discovery + supervised anchor** | GEO, open |
| GSE32707 | Dolinay et al., sepsis / sepsis-ARDS | Illumina HT-12 v4 (GPL10558) | ~190, day 0 and day 7 | **Transfer cohort 1** (has a longitudinal axis) | GEO, open |
| GSE66890 | Kangelaris et al. | Affy Human Gene 1.0 ST (GPL6244) | 28 sepsis, 29 sepsis-induced ARDS | **Transfer cohort 2** | GEO, open |
| MGH-COVID | Filbin et al. 2021, MGH ED | Olink Explore 1536 | 384 patients (306 COVID+, 78 controls), days 0/3/7, ~1,470 proteins | **Proteomic + longitudinal cohort** | olink.com/mgh-covid-study, free registration |

**Verify before building on it (hour 1, not day 3):** the MARS endotype labels for GSE65682 may live in the GEO `pData` characteristics fields, or may only be in the Scicluna et al. (Lancet Respir Med 2017) supplement. Print `pData(gse)` column names and unique values immediately and report what is actually there. If labels are absent, the supervised anchor step (§3.3) is dropped and the project proceeds unsupervised-only — this is a planned fallback, not a crisis.

**Verify likewise:** the MGH Olink release ships as a long-format NPX file and a separate clinical file. Do not assume column names. Print `head()` and full column lists for both and write the observed schema into `config/`. Expected content includes per-subject ordinal acuity at days 0/3/7/28, age category, BMI category, comorbidity flags, and selected clinical labs — confirm before writing any join.

### 3.2 Ingest and QC

- GEO via `GEOquery::getGEO`. Check whether the series matrix is already log2 — do not double-log. Rule: if `max(exprs) > 50`, apply `log2(x+1)`.
- Probe → gene symbol collapse: keep the probe with maximum IQR per gene. Document how many genes survive.
- Per-cohort: PCA colored by any available batch/scan-date variable. If PC1 tracks batch rather than phenotype, **say so in the report and proceed with within-cohort analysis only.**
- Olink NPX: pivot long → wide (`SampleID × Assay`). Handle assay warnings and LOD flags per Olink convention — count and report the fraction below LOD per protein; drop proteins >50% below LOD in the day-0 samples.

### 3.3 Supervised anchor (conditional on labels existing)

On GSE65682 only:
- Train a regularized multinomial model (elastic net, `glmnet` / `sklearn` `LogisticRegression(penalty='elasticnet')`) to predict MARS endotype from expression.
- Nested CV: outer 5-fold for performance, inner 5-fold for λ/α.
- Output: per-class AUC, and a **fixed, frozen** signature gene list written to `results/mars_signature.json` with coefficients. This file is never regenerated after §3.5 begins.

Purpose: this gives a published, externally-defined endotype axis to sanity-check the unsupervised clusters against. If the unsupervised clusters recover the MARS labels, that is strong evidence they are real biology.

### 3.4 Unsupervised discovery

On GSE65682, top ~2,000 genes by MAD:
- Consensus clustering (`ConsensusClusterPlus`, 1000 resamples, 80% subsampling, k = 2..6).
- Select k by **proportion of ambiguous clustering (PAC)**, not by eyeballing the consensus heatmap.
- Report silhouette width and cluster sizes. A cluster with n < 30 is not a finding.
- Cross-tabulate against MARS labels (adjusted Rand index) if available.

### 3.5 Transfer — the step that prevents self-deception

**Do not merge cohorts and cluster the merged matrix.** Batch effects will manufacture clusters that perfectly separate by cohort, and this is the single most common failure in this literature.

Instead:
- Freeze the discovery signature (from §3.3 or, if unsupervised-only, the top differentially expressed genes between discovery clusters, capped at 100 up / 100 down).
- Score each held-out cohort (GSE32707, GSE66890) with **single-sample scoring** — `GSVA::gsva(method="ssgsea")` or a fixed-coefficient linear score — with **zero refitting**.
- Test: does the score separate ARDS from sepsis-without-ARDS? Does it associate with 28-day mortality where available? Report effect sizes with CIs, not just p-values.
- On GSE32707 specifically: does an individual patient's score change from day 0 to day 7, and does the direction of change track outcome? This is the trajectory-mapping claim in miniature.

### 3.6 Pre-registered pathway scores — and an honest limit on what they can show

**Write `config/genesets.yaml` before running §3.4. Commit it. Do not edit after.**

#### Read this before interpreting anything in this subsection

**Transcript abundance is a poor readout for p38/MK2 activity, and this is why S5 is a secondary arm rather than the headline.** The p38α:MK2 axis is regulated post-translationally: MK2 is activated by phosphorylation, and its principal effect is phosphorylating tristetraprolin (TTP/ZFP36), which changes the *stability* of TNF and other cytokine mRNAs. None of that is visible in *MAPK14* or *MAPKAPK2* transcript levels, which are largely constitutive. A whole-blood transcriptome cannot measure whether p38 is *active*.

What transcriptomics can see is the **downstream consequence** — elevated TNF, IL6, CXCL8, IL1B. But those are equally downstream of NF-κB and most other inflammatory inputs, so they are not specific to p38. This is exactly why the comparator test below is built in, and it is why a null result is the *expected* outcome rather than a disappointment.

**Consequence for the project:** report the transcript score, report it against comparators, and interpret it as "consistent with, but not evidence of, p38-axis dominance." **Protein space (§3.7) is the primary venue for S5** — measured cytokine output is a closer, though still indirect, readout than transcript abundance. Anything stronger would need phosphoproteomics, which no public cohort here provides. Say that in the limitations section; naming the assay you would need is a stronger signal than pretending the one you have suffices.

Primary — p38α/MK2 axis (mechanism-specific):
```
MAPK14, MAPKAPK2, MAPKAPK3, MAPKAPK5, HSPB1, ZFP36, ELAVL1,
MAP2K3, MAP2K6, MAP3K7, DUSP1, DUSP10, TNF, IL6, CXCL8, IL1B, PTGS2
```
`ZFP36` (tristetraprolin) is the mechanistic centerpiece — MK2 phosphorylates TTP, which stabilizes TNF mRNA. It is the single gene most specific to the axis GEn1E modulates, as opposed to generic inflammation. Report it individually as well as within the set score.

Comparators (to show the p38 score is not just "inflammation"):
- MSigDB Hallmark `TNFA_SIGNALING_VIA_NFKB`
- MSigDB Hallmark `INFLAMMATORY_RESPONSE`
- MSigDB Hallmark `INTERFERON_ALPHA_RESPONSE` (a negative control — should *not* track the p38 endotype specifically)

Scoring: ssGSEA per sample. Compare across endotypes with Wilcoxon + Cliff's delta, and report the p38 score's separation *relative to* the generic comparators. If p38 separates no better than `INFLAMMATORY_RESPONSE` — **the expected outcome, per the caveat above** — the conclusion is "the hyperinflammatory class is broadly inflammatory; blood transcriptome cannot resolve whether p38/MK2 specifically dominates." Write exactly that. It is the correct scientific statement and it demonstrates you understand the assay's limits, which is worth more than a fabricated positive.

### 3.7 Olink longitudinal layer

This is the most GEn1E-shaped part of Layer 1: protein-level, longitudinal, with clinical outcomes.

- Cluster day-0 NPX (COVID+ patients only; controls held out as a reference distribution). Same consensus-clustering protocol as §3.4.
- **Bonus that matters:** Olink Explore 1536 contains IL6 and TNFRSF1A (sTNFR1) — two of the three variables in Calfee's parsimonious hyper/hypo-inflammatory classifier. The third (bicarbonate) is in the clinical file if present; otherwise substitute an available acidosis marker. Implement the 3-variable classifier as a **published external comparator** and cross-tab it against your own clusters. Concordance with a published classifier is the strongest single validation available in Layer 1.
- Trajectory: does day-0 endotype predict acuity at days 3, 7, 28 (ordinal outcome — proportional odds model)? Does endotype membership itself switch between day 0 and day 7?
- p38-axis proxy in protein space: IL6, TNF, CXCL8, IL1B, CCL2, IL10, TNFRSF1A, IL18, HGF, LCN2. Same comparator discipline as §3.6.

### 3.8 Layer 1 deliverables

- `figures/L1_consensus_k_selection.png` — PAC curve + consensus matrices
- `figures/L1_endotype_heatmap.png` — discovery clusters, annotated with MARS labels and outcome
- `figures/L1_transfer_scores.png` — signature score distributions across all held-out cohorts
- `figures/L1_pathway_scores.png` — p38 score vs. comparators, by endotype
- `figures/L1_olink_trajectory.png` — acuity trajectory by day-0 proteomic endotype
- `results/L1_summary.json` — k, silhouette, ARI vs MARS, all effect sizes with CIs

---

## 4. LAYER 2 — Clinical trajectories (fully open data)

### 4.1 Data — no credentialing, no DUA, no application

All three sources below are open access. Nothing here requires a PhysioNet credentialing application, a CITI course, an institutional affiliation, or a signed data use agreement. **This removes the single largest schedule risk in spec v1.1 (old §6.1) and buys back roughly four days.**

| Source | n | Structure | Access | Role |
|---|---|---|---|---|
| **PhysioNet/CinC Challenge 2019** | 40,336 ICU patients (set A 20,336, set B 20,000) | One `.psv` per patient, one row per hour, 40 variables | **Open Access, CC BY 4.0.** 42 MB. `wget -r -N -c -np https://physionet.org/files/challenge-2019/1.0.0/` | **Primary cohort.** Hourly trajectories at real scale |
| **eICU-CRD Demo** | ~2,500 unit stays, ~20 hospitals | Full eICU relational schema | Open access | **Mortality + ventilation cohort** — supplies the two things Challenge 2019 lacks |
| **MIMIC-IV Clinical Database Demo** | 100 patients | Full MIMIC-IV schema | Open access | Schema familiarity; pipeline smoke tests |

The two Challenge 2019 training sets come from **two different hospital systems**. Train on set A, validate on set B. That is a free, built-in, cross-institution external validation — a *better* generalization story than the eICU plan in v1.1, and it costs nothing.

**Set `Unit1`/`Unit2` aside as a covariate, not a feature.** They are administrative MICU/SICU identifiers; clustering on them recovers unit assignment, not biology.

#### Challenge 2019 variables

Vitals (hourly): `HR`, `O2Sat`, `Temp`, `SBP`, `MAP`, `DBP`, `Resp`, `EtCO2`

Labs: `BaseExcess`, `HCO3`, `FiO2`, `pH`, `PaCO2`, `SaO2`, `AST`, `BUN`, `Alkalinephos`, `Calcium`, `Chloride`, `Creatinine`, `Bilirubin_direct`, `Glucose`, `Lactate`, `Magnesium`, `Phosphate`, `Potassium`, `Bilirubin_total`, `TroponinI`, `Hct`, `Hgb`, `PTT`, `WBC`, `Fibrinogen`, `Platelets`

Demographics: `Age`, `Gender`, `Unit1`, `Unit2`, `HospAdmTime`, `ICULOS`

Label: `SepsisLabel` (Sepsis-3, shifted 6h early)

**`HCO3` is present.** That is the third variable in the Calfee parsimonious hyper/hypo-inflammatory classifier, and it means the same external comparator built in Layer 1 (§3.7) can be built here — the strongest cross-layer link available. `FiO2` and `SaO2`/`O2Sat` are present, which is what makes a hypoxemia criterion possible at all.

**What is absent, and it matters:** no mortality outcome, no PEEP, no ventilation flag, no vasopressors, no medication or treatment exposure of any kind. The only label is sepsis onset. See §4.4 and §6.1.

### 4.2 Database

**DuckDB, but the scale problem is gone.** The whole Challenge 2019 corpus is 42 MB — it fits in RAM many times over. The v1.1 warnings about 400M-row tables, filtered ingest, two-pass ordering, and BigQuery are all moot on the open path and have been removed.

Ingest is a directory scan. DuckDB reads the pipe-delimited files directly:

```python
# 01_build_duckdb.py — sketch
con.execute("""
CREATE TABLE hourly AS
SELECT *,
       regexp_extract(filename, 'p(\\d+)\\.psv', 1) AS patient_id,
       CASE WHEN filename LIKE '%training_setA%' THEN 'A' ELSE 'B' END AS hospital_system
FROM read_csv_auto('data/challenge2019/**/*.psv', delim='|',
                   nullstr='NaN', filename=true, union_by_name=true);
""")
```

Confirm the delimiter, the `NaN` null encoding, and the exact directory names from the downloaded tree before running this — do not assume the paths above.

eICU demo and MIMIC demo load as ordinary CSVs. For the eICU demo, the tables needed are `patient` (has `unitdischargestatus`, `hospitaldischargestatus`, `hospitalid`), `vitalPeriodic`, `lab`, `apachePatientResult`, `respiratoryCare`, and `treatment`.

### 4.3 Cohort definition — `03_cohort_hypoxemic.sql`

**Read §6.3 before writing this. The framing matters more than the SQL, and on the open path it matters more than it did in v1.1.**

Without PEEP or a ventilation flag, Berlin criteria cannot be approximated even loosely. What is definable here is **hypoxemic respiratory failure among ICU patients**, and the report must call it exactly that. Do not use the word ARDS for the Challenge 2019 cohort anywhere.

Inclusion:
- Age ≥ 18
- `ICULOS` ≥ 24 hours of record
- At least one hour with **SpO₂/FiO₂ ≤ 235** (a validated surrogate for P/F ≤ 200) using `O2Sat` and `FiO2`, or `SaO2`/`FiO2` where arterial values exist
- A recorded non-null `FiO2` at that hour — itself a weak proxy for being on oxygen support, since `FiO2` is charted mainly for supported patients
- Index time (T0) = first qualifying hour

Exclusion:
- Qualifying hour beyond ICU hour 168
- `O2Sat` ≥ 97 at the qualifying hour — the S/F ratio is uninformative once the oximeter saturates, and including these patients injects noise

**Expected yield:** unknown, and that is the point — do not anchor on a number. Report what the filters actually produce, with a **CONSORT-style attrition table** (`results/L2_attrition.csv`). If the cohort falls below ~1,500 patients, relax the S/F threshold to ≤ 315 (surrogate for P/F ≤ 300, the mild-ARDS boundary) and document the change.

Missingness will be severe — `FiO2` is one of the sparser columns. Report the fraction of patients with any recorded `FiO2` before filtering on it; if it is very low, the cohort is a selected subpopulation and the report must say so.

### 4.4 Features and the outcome problem

**Static (T0 to T0 + 24h):** age, sex, `HospAdmTime`, and first-24h summaries of `WBC`, `Platelets`, `Creatinine`, `Bilirubin_total`, **`HCO3`**, `Lactate`, `pH`, `BaseExcess`, `BUN`, `Glucose`, `PTT`, `Fibrinogen`, plus S/F ratio and `MAP`.

`HCO3` is deliberately called out — third variable in the Calfee classifier, and the anchor for the cross-layer comparator in §5.

**Longitudinal (T0 to T0 + 72h, 6-hour bins → 12 timesteps):** the hourly data aggregates cleanly into bins with min/max/mean/last. **Always emit a parallel binary missingness mask.** Whether a lab was drawn is informative, and forward-filling without a mask fabricates trajectories. On this dataset the labs are sparse enough that the mask may carry more signal than the values — worth reporting as a finding in its own right.

**No vent features, no exposure features.** `PEEP`, tidal volume, plateau pressure, driving pressure, vasopressors, corticosteroids, NMB and RRT are all unavailable. The v1.1 feature schema is cut accordingly.

#### The outcome problem — read carefully

**Challenge 2019 has no mortality label.** The only outcome is `SepsisLabel`. This is the real cost of going open, and it must not be papered over.

Available proxy outcomes, in descending order of defensibility:
1. **Sepsis onset / time-to-sepsis** — a genuine labelled clinical event, Sepsis-3 adjudicated. Use as the primary endpoint.
2. **`ICULOS` at final record** — a weak, censored LOS proxy. Report, but never as a headline.
3. **Trajectory of organ-dysfunction proxies** (creatinine, bilirubin, platelets, lactate over 72h) — not an outcome, but an interpretable severity path.

**Mortality comes from the eICU demo instead.** Re-derive the same static features there, apply the frozen cluster model, and test the mortality gradient on ~2,500 stays across ~20 hospitals. The n is thin and the CIs will be wide — report them wide. A small, honestly-bounded mortality replication is worth more than a large unlabelled claim.

If the eICU demo yields too few qualifying patients for a stable estimate, say so and let the mortality question stand unanswered. **Do not substitute a self-defined "severity outcome" and present it as though it were mortality.**

### 4.5 Models

**(a) Reference assignment — the non-portable baseline.** Gaussian mixture / latent class on log-transformed, **cohort-standardised** static features. k = 2..4, BIC-selected, 100 random restarts. Expected: a hyperinflammatory-like class with lower `HCO3`, higher lactate, higher creatinine and bilirubin, and a higher sepsis rate, matching the published two-class structure. Report class-conditional means as a forest plot.

This is deliberately the *unusable* model — it needs the cohort mean of every variable at inference. It exists to generate the reference labels that (a2) has to reproduce. Do not present it as a result; present it as the baseline whose deployment problem the project addresses.

**(a2) The portable classifier — P1, the headline experiment.** Train a model to reproduce (a)'s assignment using **only per-patient values at inference**: raw measurements, absolute thresholds, and ratios between a patient's own variables. No z-scores against a cohort. No population means, medians, or quantiles anywhere in the inference path.

Candidate forms, in increasing order of deployability:
- Logistic regression on raw (untransformed, unstandardised) values
- A shallow decision tree with absolute cut-points — the most deployable form, and directly readable as a bedside rule
- Gradient boosting on raw values as the performance ceiling

Report κ against the reference, the confusion matrix, and calibration. **Write an explicit assertion that no fitted scaler, cohort statistic, or population parameter appears in the inference function** — a `assert_no_population_stats()` check over the serialized model. This is the technical core of the whole project; make it verifiable rather than asserted.

**(a3) Parsimony curve — P3.** Agreement as a function of variable count, greedily selected. Report the smallest set retaining acceptable κ, and **cross-reference §1.5 coverage numbers**: a variable present in 12% of patients is not deployable regardless of how much it helps. Plot κ against both variable count and variable availability.

**(b) Trajectory HMM.** `hmmlearn.GaussianHMM` over the 12×features tensor, 3–6 hidden states, BIC-selected on a held-out split. Deliverables: transition matrix, state occupancy by endotype, state-at-24h → sepsis-onset mapping. This instantiates the "hidden Markov models" and "trajectory mapping" language from GEn1E's platform page — say so in the README.

**(c) Multi-task encoder.** The port of the Hydra/MAGE architecture; frame it that way. Shared GRU or temporal-conv encoder over the 72h tensor → shared latent → heads for:
- sepsis onset within the window (binary)
- time-to-sepsis (regression, censored)
- endotype membership (categorical, from (a))
- inflammatory-axis score (regression)

Per-head label masking, exactly as in Hydra. Report per-head metrics against single-task baselines — the multi-task claim needs the ablation or it is decoration. **Also report an XGBoost baseline on flattened features.** If XGBoost matches the deep model, say so; GEn1E lists gradient boosting on their own platform, so a well-executed GBM baseline is not a weakness.

**(d) Cross-hospital transfer — P2, the second headline experiment. No longer optional.**

1. **Set A → set B.** Fit (a) and (a2) on hospital system A only. Freeze both. Apply to set B with zero refitting. Report:
   - Agreement of the portable classifier with a set-B-native reference assignment
   - **Calibration drift** — does the portable model's threshold still sit in the right place on B? Reliability diagram for both.
   - The *comparison that matters*: how much does the standardised reference model degrade on B versus the portable one? If the portable model transfers better because it never depended on A's population statistics, **that is the project's central result.** If it transfers worse, that is a real and reportable negative finding.
2. **→ eICU demo.** Same frozen models, different schema entirely, plus the mortality endpoint from §4.4. A harder test: variable definitions and units will not line up cleanly, and that friction is itself informative about deployability.

**Before running either, check the §1.5 scoping result:** if sets A and B differ substantially in coverage or case mix, a performance drop is a missingness artifact, not a transfer failure. Report which one you are measuring. Getting this distinction wrong would invalidate P2.

### 4.6 Enrichment simulation — C8, the drug-development component

The single component that speaks most directly to a Phase 2 company, and it is cheap.

**Prognostic enrichment is a power calculation, not a causal claim.** Given observed event rates in each endotype, compute the sample size a trial would need under (a) unstratified enrollment and (b) enrollment restricted to the higher-event-rate endotype. Standard two-proportion power calculation. Report required N, the screening burden (how many patients you must screen to enroll one), and the net trial-size effect.

**Why this is safe when §6.8's HTE analysis is not.** Prognostic enrichment asks "which patients have more events," which is directly observable. **Predictive** enrichment asks "which patients benefit more from treatment," which requires a treatment-effect estimate and is confounded by indication in observational data. Do the first. Name the second in next steps, with the reason. **Do not blur them** — the distinction is exactly the kind of thing a clinical-development team will probe in an interview, and getting it right is worth more than the analysis itself.

Extend with the §4.5c degradation results: if the classifier degrades under realistic missingness, screening accuracy falls, and enrichment benefit falls with it. **Report required-N as a function of classifier accuracy.** That connects the deployment stress-test to a number a trial designer actually cares about, and it is the strongest single figure in the project.

Caveats to state inline: event rates from a non-ARDS-adjudicated cohort with sepsis onset as the endpoint, not mortality (§6.1, §6.3). The simulation demonstrates the method; the numbers are illustrative.

### 4.7 Iterative refinement loop — C6

The JD asks for the ability to *design computational pipelines where predictions, observed outcomes, and new data feed back into subsequent model development or refinement.* Nothing else in this spec covers it, and almost no portfolio project does. Cheap to build, and distinctive.

Implement one full documented cycle:

1. **Round 1** — train the classifier on set A; predict on set B.
2. **Confront predictions with outcomes** — where does it fail? Stratify errors by missingness pattern, case mix, and endotype. Produce an error taxonomy, not just an aggregate metric.
3. **Derive a refinement from what round 1 revealed** — e.g. add a mask-derived feature for the failing pattern, reweight, adjust the decision threshold for calibration drift, or drop a variable whose availability is too low to deploy. **The refinement must be traceable to a specific round-1 finding.**
4. **Round 2** — retrain with the refinement, using a *held-out* slice never touched in rounds 1–2 diagnosis.
5. **Report the delta honestly**, including if it is null or negative.

**Guard against the obvious failure:** this is a loop, not a licence to tune on the test set. Reserve a final holdout touched exactly once, at the end. Log every refinement decision with its triggering observation in `results/refinement_log.md` — that log *is* the deliverable for C6, more than the performance change is.

Write it as a `RefinementRound` dataclass so the structure is legible in code, not just in prose. An interviewer who opens the repo should see the loop architecture immediately.

### 4.8 Layer 2 deliverables

- `results/L2_attrition.csv` — CONSORT table
- `figures/L2_cohort_table1.png` — Table 1 by endotype
- `figures/L2_endotype_profiles.png` — class-conditional feature means
- `figures/L2_hmm_states.png` — transition matrix + state-outcome mapping
- `figures/L2_model_performance.png` — multi-task vs single-task vs XGBoost
- `figures/L2_eicu_replication.png` — mortality gradient, MIMIC vs eICU, per-hospital forest

---

## 5. Where the layers meet

The one cross-layer claim worth making, and the one to state carefully:

> Layer 1 finds molecular endotypes with differential p38/MK2-axis activity. Layer 2 finds clinical endotypes with differential event rates. **These are different cohorts and cannot be joined at the patient level.**

The defensible bridge is *concordance of structure*: both layers, independently, recover a two-class hyper/hypo-inflammatory split, and both can be scored by the same published parsimonious classifier — Layer 1 in protein space (IL6 and TNFRSF1A from Olink Explore, §3.7), Layer 2 in clinical-lab space (`HCO3`, §4.4). The shared bicarbonate variable is what makes this more than an analogy.

**v1.3 note:** with the thesis now centred on portability, this bridge is supporting material rather than a headline claim. Keep it to a paragraph and one figure. The temptation is to over-invest here because it feels like the intellectually satisfying part; resist it — P1–P3 are what the reader will judge.

**Note the asymmetry introduced in v1.2:** the v1.1 version of this bridge also claimed both layers show higher mortality in the hyperinflammatory class. Layer 2 can no longer support the mortality half (§6.1). The concordance claim now rests on cluster structure and classifier agreement, not on a shared outcome. Weaker, and it must be stated as weaker.

**Do not claim a patient-level link.** Anyone at GEn1E will know it isn't there, and overclaiming it costs more than the finding is worth.

---

## 6. Devil's advocate — failure-mode register

Ordered by expected cost. Each has a **kill criterion**: the point at which you stop and take the fallback rather than keep digging.

### 6.1 🔴 No mortality outcome in the primary Layer 2 cohort
**Likelihood: certain. This is the price of the open-data decision and the most important risk in v1.2.**

Challenge 2019 labels only sepsis onset. Success criterion **S4** as written in v1.0 called for "a mortality gradient," and on the primary cohort that gradient cannot be computed at all. Anyone reading the report who knows the dataset will know this within seconds.

**Mitigation:**
- Restate S4 explicitly as *"endotypes with a differential clinical-event gradient"*, with sepsis onset as the primary endpoint and mortality carried by the eICU demo alone (§4.4).
- Report the eICU mortality replication with wide CIs and an explicit note on the small n. Wide and honest beats narrow and fabricated.
- Say in the README, in one plain sentence, that the primary cohort has no mortality label and why the project accepted that tradeoff.

**What NOT to do:** do not construct a composite "severity outcome" from organ-dysfunction variables and let it read like mortality. That is the single most tempting move here and it is the one that would actually damage credibility.

**Kill criterion:** none — this is a known, accepted constraint rather than a failure. But if the eICU demo cannot support a stable mortality estimate, the honest answer is to leave the mortality question open, not to manufacture a substitute.

### 6.1b 🟢 Credentialing delay — ELIMINATED in v1.2
This was the largest risk in v1.0/v1.1 and no longer exists. No application, no CITI course, no DUA, no waiting. Recorded here so the removal is visible rather than silently dropped.

### 6.2 🟡 Cohort definition — much cheaper now, but sparser
**Downgraded from 🔴. The MIMIC itemid swamp is gone; a missingness problem replaces it.**

Challenge 2019 is a flat hourly table with named columns — no itemid archaeology, no joining ABGs to concurrent vent settings, no P/F ratios of 0.4 and 3,000 in the same table. The v1.1 fallback (S/F ratio) is now the *primary* definition, so the fallback has already been taken.

**What replaces it:** `FiO2` is one of the sparser columns, and the cohort is defined on it. Requiring a recorded `FiO2` silently selects for oxygen-supported and more heavily monitored patients.

**Mitigation:** before filtering, report the fraction of the full 40,336 with any recorded `FiO2`, and compare age/sex/`ICULOS`/lactate between those with and without. If they differ substantially, the cohort is a selected subpopulation and Table 1 must show it. Sanity-plot the S/F distribution before filtering. Drop hours with `O2Sat` ≥ 97 (§4.3) — a saturated oximeter makes S/F meaningless.

**Kill criterion:** 1.5 days, down from 3.

### 6.3 🔴 The cohort is definitely not ARDS — and this got worse in v1.2
**Likelihood: certain. Upgraded from 🟠 because the gap widened.**

Berlin requires bilateral opacities on imaging, a PEEP ≥ 5 criterion, and exclusion of cardiogenic pulmonary edema. On the open path there is no imaging, **no PEEP, and no ventilation flag at all**. What §4.3 defines is hypoxemic respiratory failure by oximetry surrogate — further from ARDS than the v1.1 MIMIC definition was, and that one was already a proxy.

There is also no radiology-note fallback: the MIMIC-IV-Note sensitivity analysis offered in v1.1 is unavailable here.

**Mitigation:** never write "ARDS" about the Challenge 2019 cohort. Not in figure captions, not in the README, not in variable names. Call it hypoxemic respiratory failure. State the gap in the README's opening paragraph and in the Table 1 caption. Frame the project's ARDS relevance as coming from Layer 1's molecular cohorts and from the *structural* argument in §5 — not from the clinical cohort.

**Why this is still a strength:** naming a limitation before a reviewer finds it is the clearest available signal of someone who has done real clinical data work. Someone at a Phase 2 ARDS company will spot an unadjudicated hypoxemia proxy within thirty seconds. Beating them to it converts the weakness into evidence of judgment. Overclaiming here is the fastest way to lose the reader.

### 6.4 🔴 The endotypes are just severity
**Likelihood: high. Cost: the whole thesis.**
If your two clusters separate on "sicker vs. less sick," you have rediscovered SOFA with extra steps. This is the standard and correct criticism of the entire subphenotype literature.
**Mitigation, pre-specified:**
- Re-run clustering with all features standardized *within* SOFA quartile. Do the clusters survive?
- Report cluster separation on *ratios and shapes* (e.g. bicarbonate given creatinine, platelet trajectory slope) not levels alone.
- Report AUC for mortality from cluster membership vs. from SOFA alone vs. both. If cluster membership adds < 0.02 AUC over SOFA, say so plainly.
**Kill criterion:** none — this is a finding either way. But if severity fully explains the clusters, the report's conclusion changes to "clinical-variable endotyping in this cohort is not separable from severity; the molecular layer is where the signal is," which is honest and still interesting.

### 6.5 🟠 Batch effects manufacture fake clusters
**Likelihood: high if the discipline in §3.5 slips.**
Three GEO cohorts on three platforms. Merging and clustering will produce clusters that separate perfectly by cohort.
**Mitigation:** never cluster a merged matrix. Train on one, score fixed signatures on the others. Never run ComBat with the cluster label as a covariate — that is circular and will produce a beautiful, meaningless heatmap.
**Detection:** if any cluster is >90% one cohort, that is the failure mode, not a finding.

### 6.6 🟠 The MARS endotype labels aren't in GEO
**Likelihood: moderate.**
**Mitigation:** verify in hour 1. Fallback is unsupervised-only, which weakens S1's external anchor but does not break the project. Second fallback: extract labels from the Scicluna et al. supplement and join by GSM ID.
**Kill criterion:** 2 hours of searching, then take the fallback.

### 6.7 🟡 COVID-era proteomics ≠ non-COVID ARDS
**Likelihood: certain. Cost: moderate.**
The MGH cohort is COVID-era acute respiratory distress. Its inflammatory biology is related to but not identical to non-COVID ARDS.
**Mitigation:** frame the entire Layer 1 proteomic section as "acute hypoxemic respiratory failure," never as "ARDS." Note that the cohort's 78 COVID-negative patients presented with the *same clinical syndrome*, which lets you test whether the endotype structure is driven by etiology or by host response — genuinely interesting, and worth a figure.

### 6.8 🟠 Confounding by indication if you touch treatment effects
**Likelihood: certain if attempted.**
Corticosteroids in ARDS are given preferentially to sicker patients. A causal forest on observational ICU data will return a confident, wrong, and dangerous-sounding answer about who benefits.
**Mitigation:** this is why HTE is Layer 4 and out of scope. If tempted, it goes in "next steps" as *"a target-trial-emulation design with negative-control outcomes and E-value sensitivity analysis would be the appropriate next step"* — which demonstrates you know why you didn't do it. That sentence is worth more than a badly-executed causal forest.
**Kill criterion:** do not start Layer 4. At all.

### 6.9 🟡 Overfitting 1,470 proteins on n=306
**Likelihood: high.** Feature:sample ratio ~5:1.
**Mitigation:** all Olink modeling inside nested CV. Report permutation-test null distributions for cluster separation. Correct for multiple testing (BH) on all per-protein tests and report the number surviving, not just the top hits.

### 6.10 🟢 Compute and disk — ELIMINATED in v1.2
Every dataset now fits in RAM. Total footprint under 5GB (§0.3). The v1.0 and v1.1 discussions of decompression, filtered ingest, spill space, and BigQuery quota all belonged to MIMIC-IV and no longer apply. Recorded so the removal is visible.

The only compute item left: the multi-task encoder in §4.5c trains on GPU, which the machine has.

### 6.11 🔴 Scope creep against a hard deadline
**Likelihood: very high. This is the most probable actual cause of failure.**
Four layers, three cohorts per layer, a multi-task deep model, an HMM, and external validation — in three and a half weeks, alongside a job search.
**Mitigation:** the cut lines in §7 are load-bearing. A complete, well-documented Layer 1 with three honest figures beats a sprawling half-finished four-layer repo. **The README is the deliverable, not the code.** A hiring manager will spend four minutes on this.

### 6.12 🟡 Claude Code specific failure modes
- All downloads are manual, by James, locally. Claude Code should not attempt to fetch datasets itself.
- It will assume file layouts and column names it has not verified. Every schema in §3.1 and §4.1 is flagged "verify before building on it" — hold it to that. Print `head()` and column lists first, always.
- It will write SQL that returns 0 rows and continue happily to the next step. See §8.3.
- It will suggest running Layer 3 or 4 because they're in the background context. They are out of scope. Say no.
- **It will reach for MIMIC-IV.** MIMIC is heavily represented in its training data and it will drift toward MIMIC-shaped code, MIMIC itemids, and mimic-code concepts out of habit — none of which exist on this path. If MIMIC-IV appears anywhere outside §10, that is drift; correct it.
- It may write "ARDS" about the Challenge 2019 cohort. See §6.3. Grep the repo for the string before shipping.
- **It will drift back to "do endotypes exist."** That framing dominates the literature it was trained on, and it is the v1.0 thesis this spec replaced. If the analysis starts optimizing cluster quality rather than portability and transfer, that is drift. The deliverable is P1–P3.
- It will happily leave a fitted `StandardScaler` in the inference path and call the model portable. See §4.5a2 — make the check executable.

### 6.13 🟠 Derivative — reframed in v1.4
**Downgraded from 🔴.** v1.3 treated "this reproduces known work" as an existential risk. Under the v1.4 objective it is not: reproducing established findings accurately, with citation, is a competency demonstration (§1, C-criteria). The risk was never derivativeness — it was **overclaiming**.

**The actual risk:** presenting reproduced findings as discovery. That converts a good artifact into evidence of poor judgment, in front of readers who know this literature well.

**Mitigation:** I1. The README's opening paragraph states that the subphenotypes are established and cites them. Every reproduced result is labelled as a reproduction in its own figure caption.

**Detection:** if any figure caption or README sentence would imply first discovery to a reader who skims, rewrite it.

### 6.15 🟠 Competency coverage collapses under time pressure
**New in v1.4. Likelihood: moderate. Cost: the deliverable's whole purpose.**

The failure mode is cutting the last evidence for a competency while keeping a third demonstration of one already covered. The multi-task encoder is intellectually appealing and duplicates C-coverage; the enrichment simulation is half a day and is the sole evidence for C8.

**Mitigation:** consult the §0.4 matrix before cutting anything. The question is never "is this interesting" but "is this the only thing evidencing its row."

**Related trap:** over-investing in Layer 1 because the biology is the enjoyable part. Layer 1 carries C7 and C9 and needs two good figures, not eight.

### 6.14 🟠 The p38 measurement can't support the claim
**Likelihood: certain in transcript space. Cost: contained, because S5 is now secondary.**

Fully described in §3.6: transcript abundance cannot measure a post-translationally regulated kinase axis. This risk is listed separately so the demotion is visible in the register rather than buried in a methods subsection.

**Mitigation:** protein space primary, transcript secondary, comparators mandatory, null expected and pre-accepted. Name phosphoproteomics as the assay that would actually answer it.

**What NOT to do:** do not expand the gene set until it separates. If `config/genesets.yaml` is edited after §3.4 has run, the result is uninterpretable — the freeze exists precisely to make this impossible to do accidentally.

---

## 7. Timeline and cut lines

Deadline 2026-09-19. Assume ~24 working days from 2026-08-26, shared with a job search.

| Days | Work | Cut line |
|---|---|---|
| **0** | Repo scaffolded, `.gitignore` in place. All datasets downloaded (~4GB, one afternoon). | Nothing else happens on day 0. `genesets.yaml` is NOT written yet — see day 1. |
| **1** | **§1.5 scoping.** Purely descriptive. Produce `results/L0_scoping.md`. Then finalize the question wording, write and commit `config/genesets.yaml` and the analysis plan. | Hard stop at one day. Scoping that runs long becomes the analysis, which is the failure it exists to prevent. No models on day 1. |
| **2–6** | Layer 1: ingest, QC, discovery clustering, transfer scoring on GEO. | If GEO QC is a swamp by day 5, drop GSE66890 and run with two cohorts. |
| **7–9** | Layer 1: Olink longitudinal + pathway scores (S5 lives here, secondary). Write L1 figures. | **Layer 1 must be figure-complete by day 9.** S5 gets one figure, not three. |
| **10–12** | Layer 2: Challenge 2019 ingest, missingness audit, cohort definition, attrition table. | Cohort timeboxed to 1.5 days (§6.2). No credentialing decision point — that gate is gone. |
| **13–15** | Layer 2: features, clustering (C1), Sinha classifier implemented + own classifiers (C2). | C1/C2. |
| **16–17** | Transfer set A → set B and → eICU (C4). Missingness ablation (C5). | **Never cut.** C4 and C5 are the strongest technical signals in the project. |
| **18** | Enrichment simulation (C8). | Half a day of work, highest interview yield per hour. Do not skip. |
| **19** | Refinement loop, one full cycle (C7 → C6). | The differentiating component. Protect it. |
| **20** | HMM trajectory (C3) + multi-task encoder. | **First to cut.** Supporting evidence only. |
| **21** | Buffer. §6.4 severity checks or Layer 1 depth. | If unused, start the writeup early. |
| **22–24** | README, figures, writeup, limitations section, "next steps" covering Layers 3–4. | Non-negotiable. Reserve these days. |

**Priority order if time collapses:** §1.5 scoping > Layer 2 clustering + classification (C1, C2) > transfer (C4) > missingness ablation (C5) > enrichment simulation (C8) > refinement loop (C6) > Layer 1 endotypes + Olink (C7, C9) > HMM (C3) > multi-task model > eICU mortality > p38 arm.

**Reordered for competency coverage in v1.4.** The rule is different from earlier versions: cut whatever duplicates a competency already demonstrated, keep whatever is the *only* evidence for one. The multi-task encoder is cuttable because C-criteria are already covered elsewhere; the enrichment simulation is not, because nothing else evidences C8. Check the §0.4 matrix before cutting anything.

---

## 8. Engineering standards

### 8.1 Reproducibility
Every script runs from CLI with explicit args, seeds fixed and logged, `environment.yml` pinned. `make all` reproduces every figure from raw data. The JD names "Git/GitHub and reproducible research workflows" — this is being graded.

### 8.2 Logging
Structured logging to `logs/`. Every step logs input row count, output row count, and wall time. When something is wrong three steps downstream, this log is how you find it.

### 8.3 Assertion checkpoints — mandatory

```python
# shared/assertions.py
def checkpoint(df, name, min_rows, max_rows, required_cols=(), no_nulls=()):
    n = len(df)
    assert min_rows <= n <= max_rows, (
        f"{name}: got {n} rows, expected [{min_rows}, {max_rows}]. "
        "Do not proceed — the join is wrong."
    )
    missing = set(required_cols) - set(df.columns)
    assert not missing, f"{name}: missing columns {missing}"
    for c in no_nulls:
        assert df[c].notna().all(), f"{name}: nulls in {c}"
    logging.info(f"✓ {name}: {n} rows")
    return df
```

Required checkpoints with expected ranges:
- Challenge 2019 ingest: exactly 40,336 distinct patients, both hospital systems present, `hospital_system` split ~20,336 / ~20,000
- Patients with any recorded `FiO2`: report the number; no hard bound, but log it prominently (§6.2)
- Cohort after S/F filter: no pre-set range — report the attrition table. Flag if below 1,500 (§4.3 relaxation rule)
- Longitudinal tensor: `n_patients × 12 × n_features`, no all-NaN feature columns, mask channel present and not all-ones
- eICU demo: ~2,500 unit stays, `hospitaldischargestatus` non-null for the large majority
- GEO expression matrix: 10,000–25,000 genes post-collapse
- Olink wide matrix: 300–400 subjects × 900–1,500 proteins

### 8.4 What gets committed
Code, config, figures, `results/*.json` and `results/*.csv` containing **aggregate statistics only**. Never row-level clinical data, never an intermediate parquet — see §0.2 item 2 for why this holds even though the licences would permit it. Include a `CITATIONS.md` listing every dataset with its required citation; CC BY 4.0 makes attribution a licence condition (§0.2 item 1).

---

## 9. Final README structure (the actual deliverable)

1. **One-paragraph thesis.** What was built and what it was tested against. **State in this paragraph that the subphenotypes are established and cited (I1).** Do not open with a discovery claim.
2. **Headline figure** — required trial N as a function of classifier accuracy under realistic missingness (§4.6). It connects modelling to a decision, which is what a clinical-development reader is looking for.
3. **What the data is and isn't** — §6.3 and the §6.1 no-mortality constraint, up front.
4. **Pipeline overview** — one diagram, the whole flow, including the §4.7 refinement loop. This is the figure an interviewer will ask you to walk through; design it for that.
5. **Endotype discovery and classification** (C1, C2) — including the reproduced Sinha classifier.
6. **Generalization** (C4) — transfer, calibration, case-mix vs measurement-practice decomposition.
7. **Deployment stress-test** (C5) — the degradation curve.
8. **Enrichment simulation** (C8) — with the prognostic-vs-predictive distinction stated.
9. **Iterative refinement** (C6) — link to `refinement_log.md`.
10. **Molecular layer** (C7) — 2 figures. p38 arm with the §3.6 assay caveat *in the section*.
11. **Limitations** — specific, unhedged.
12. **What I chose not to do and why** — predictive enrichment (confounding), Layer 3, the credentialed path. **This section is an interview asset; write it deliberately.**
13. **Reproduce it** — `make all`, `CITATIONS.md`, source links.

Keep it under 1,500 words. Every figure captioned so it stands alone.

---

## 10. Optional stretch — the credentialed path

Deliberately placed last, after the deliverable, so it cannot become a distraction. **Do not start this before the report in §9 is finished.**

If time remains after day 24, or if this project continues past the application, MIMIC-IV v3.1 and the full eICU-CRD would materially strengthen it: a real Berlin-criteria proxy with PEEP and a ventilation flag (§6.3), genuine 28-day mortality and ventilator-free days (§6.1), treatment and exposure variables enabling the Layer 4 HTE work (§6.8), and 200-hospital external validation.

The v1.1 spec text covering this path — credentialing steps, table skip-list, filtered ingest, the mimic-code DuckDB build, BigQuery — is preserved in git history at the v1.1 tag. Retrieve it from there rather than reconstructing it.

Two honest notes. First, the credentialing application takes minutes to submit and its review clock runs independently of anything else, so **submitting it on day 0 costs nothing and preserves the option** even under the open-data plan. Second, if it lands mid-project, resist rebuilding Layer 2 on MIMIC — the open pipeline will already be working, and a half-migrated Layer 2 at day 20 is the worst of both.

---

## 11. Interview material — what each component lets you say

The point of §0.4 is that every component maps to a question you will be asked. Drafting these *before* building keeps the work pointed at the goal; revise them with real numbers after.

**"Tell me about a project."**
One pipeline, six data sources, three modalities. Open with the problem framing — established subphenotypes, unclear whether they survive real screening conditions — then the enrichment number. Ninety seconds. Do not narrate the layers in order.

**"How do you validate a model?"** (C4)
Frozen model, held-out hospital system, calibration as well as discrimination, and the decomposition of transfer failure into case mix versus measurement practice. That last distinction is the answer most candidates don't have.

**"How do you handle missing data?"** (C5)
Mask channels rather than silent imputation, because missingness in ICU data encodes ordering behaviour. Then the ablation curve — you can state how performance degrades as variables disappear, which is a quantitative answer to a question usually answered with a technique name.

**"Have you worked with clinical/ICU data?"** (C9)
40,336 ICU patients, hourly. Then go straight to a specific difficulty — S/F ratio being uninformative once the oximeter saturates, or the cohort selection induced by requiring a recorded FiO₂. Specificity is the signal.

**"How do you connect ML results to biology?"** (C7)
The p38 story, including why blood transcriptome cannot resolve pathway activity for a post-translationally regulated kinase, and that phosphoproteomics is the assay that would. **Knowing what your data cannot show is the strongest available signal of scientific maturity** — it will land better than any positive result.

**"Tell me about a time you were wrong."**
The prior-art review that closed off two successive theses. You scoped against the literature, found the question was answered, and narrowed rather than proceeding. That is the answer, and it is true.

**"How would you improve a model in production?"** (C6)
Walk the §4.7 loop: predictions confronted with outcomes, errors stratified into a taxonomy, refinement traceable to a specific finding, final holdout touched once. This maps directly onto their JD language about feedback into subsequent refinement.

**Questions to ask them.** How do they currently select patients for GEn-1124 trials, and what data is available at screening? Given dilmapimod's phase IIa was in unselected at-risk patients, how does their enrichment strategy differ? Does GRID's endotyping run on trial-grade or routine-care data — and how do they handle the gap?

**Résumé line, draft:** *"Built an end-to-end ICU endotyping pipeline (40k patients, 6 public sources) reproducing published ARDS subphenotypes; quantified classifier degradation under realistic EHR missingness and cross-hospital shift, and simulated the trial-enrichment impact."* Tighten with real numbers once you have them.
