# Prior Art Review & Reading List

**Companion to:** `ards-endotype-project-spec.md` v1.3
**Purpose:** establish what is already known before writing code, and give Claude Code a tiered reading list.
**Bottom line up front:** the v1.3 thesis needs narrowing again. Details in §1.

---

## 1. The correction — read this first

**v1.3 claimed portability was the open gap. It is substantially solved.**

The spec's §0.1 argued that LCA subphenotype models can't deploy because they need population-level standardization, citing the prospective-validation literature. That framing is real but incomplete, and I missed two papers that close most of it:

**Sinha, Delucchi, McAuley et al., *Lancet Respir Med* 2020;8:247–257** — "Development and validation of parsimonious algorithms to classify ARDS phenotypes." Derived on three ARDSNet RCTs (ARMA, ALVEOLI, FACTT; n=2,022), validated on SAILS (n=715). Produced classifier models using three or four variables. The published 3-variable form is a plain logistic equation on raw values:

```
P(hyperinflammatory) = -11.9593 + 1.0138·ln(IL-6 + 1) - 0.2436·HCO3 + 1.1913·ln(sTNFR1 + 1)
```

**This is already fully portable in the sense v1.3 proposed as novel.** Fixed coefficients, per-patient values, log transforms, absolute intercept. No cohort mean anywhere. P1 as written is a solved problem.

**Sinha, Churpek, Calfee, *AJRCCM* 2020;202(7):996–1004** — "Machine Learning Classifier Models Can Identify ARDS Phenotypes Using Readily Available Clinical Data." Gradient-boosted trees using only routine clinical variables available at ED/ICU admission — no protein biomarkers at all. AUC 0.95 (95% CI 0.94–0.96) in validation; replicated at AUC ~0.94 in both ALVEOLI and FACTT. It also recovered the differential treatment responses seen in the original LCA analyses.

**Sinha et al., *Lancet Respir Med* 2022** — external validation of those classifier models in two observational cohorts (EARLI n=335, VALID n=452) plus LUNG SAFE, with LCA assignment as reference.

So: parsimonious ✓, biomarker-free ✓, routine-variables-only ✓, externally validated ✓, treatment-effect-recovering ✓.

### What this means

Three revisions in, the literature has eaten the thesis twice. That is not a failure of the revisions — it is what prior art review is *for*, and it is the reason this document should have been written before v1.0 rather than after v1.3. Worth saying plainly in the README as a note on method: the project scoped itself against the literature before committing, and narrowed twice as a result. That is a more honest and more impressive story than a project that never checked.

### The gap that actually survives

Every validation above ran in **curated research cohorts** — ARDSNet RCTs, EARLI, VALID, LUNG SAFE — where variables were measured on protocol, adjudication was performed, and missingness was low and largely by design.

Nobody has publicly stress-tested these classifiers under the conditions a real screening deployment faces:

- **Realistic EHR missingness.** The classifiers assume the variables exist. In routine ICU data they frequently don't, and missingness in EHRs is well established as informative rather than random — driven by ordering behaviour and care-process intensity, not by patient biology alone.
- **Unadjudicated cohorts.** Real screening runs on hypoxemic respiratory failure before anyone has confirmed ARDS.
- **Cross-hospital shift** in routine care data rather than in multi-site research protocols.

**Proposed narrowed thesis:**

> **Do published ARDS subphenotype classifiers survive the conditions of routine, sparsely-measured, unadjudicated ICU data — and how gracefully do they degrade as the variables they need go missing?**

This is defensible because: it takes the published models as given rather than competing with them; it uses Challenge 2019's sparsity as the instrument rather than apologising for it (reversing §6.2); the A/B hospital split tests shift under routine-care conditions; and it answers the question GEn1E faces operationally — *can we actually screen patients for this trial with data we will really have?*

### Why this lands especially well at GEn1E

p38 inhibition in ARDS has been tried essentially once. A phase IIa trial of dilmapimod (SB-681323) in trauma patients at risk for ARDS showed the drug was safe and modestly reduced CRP, but was not powered for ARDS incidence or severity. Losmapimod separately failed its primary endpoint in COPD despite good tolerability. The recurring diagnosis across the p38 field is **unselected patient populations** — the compounds were given to everyone with the syndrome.

That is precisely the problem GEn1E's precision positioning claims to solve. A project asking "can we reliably identify the enrichment population using data a real ICU actually captures" is directly upstream of whether their Phase 2 can succeed. Frame it that way in the README.

---

## 2. Prior art map

| Claim | Status | Key work |
|---|---|---|
| Two inflammatory subphenotypes exist in ARDS | **Settled.** Consistent across 5+ RCT cohorts | Calfee 2014; Famous 2017; Sinha 2018 (SAILS) |
| They carry differential mortality | **Settled** | Calfee 2014 and subsequent |
| They show differential treatment response | **Established**, retrospectively | Simvastatin (Calfee 2018); PEEP (ALVEOLI); fluid strategy (FACTT); corticosteroids in COVID |
| Parsimonious 3–4 variable classifiers work | **Solved** | Sinha 2020 *Lancet Respir Med* |
| Classification without protein biomarkers | **Solved** | Sinha 2020 *AJRCCM* (GBM, AUC 0.95) |
| Validation in observational cohorts | **Done** | Sinha 2022 *Lancet Respir Med* (EARLI, VALID, LUNG SAFE) |
| Prospective bedside identification | **Done, 2026** | PHIND (*Lancet Respir Med*); SPARC (ISPY COVID Network, 17 hospitals) |
| Replication in MIMIC-IV / eICU | **Published 2026** | Frontiers in Medicine 18-variable external validation |
| Subphenotypes extend beyond ARDS to acute respiratory failure generally | **Established** | Multiple; relevant because our cohort is unadjudicated |
| **Robustness to realistic EHR missingness** | **Open** | ← the gap |
| **Degradation under routine-care cross-hospital shift** | **Largely open** | ← the gap |
| Informative missingness in ICU data | Established generally, **not applied to subphenotyping** | See §3 Tier 2 |
| p38/MK2 as ARDS target | Preclinically strong, **one underpowered phase IIa** | Christie 2015 (dilmapimod) |

**Implication for the spec:** S1 (endotypes exist) and P1 (portable classifier) both move to "replicate the published result as a validation step, cite it, spend one figure." The novel contribution concentrates entirely in the missingness/degradation analysis.

---

## 3. Reading list

Access is marked: 🟢 open · 🟡 PMC/preprint available · 🔴 paywalled (read abstract + secondary descriptions)

### Tier 1 — read before writing any Layer 2 code

These define the reference models the project must reproduce and stress-test. Without them the project has no baseline.

1. **Sinha, Delucchi, McAuley et al. (2020).** Development and validation of parsimonious algorithms to classify ARDS phenotypes. *Lancet Respir Med* 8:247–257. 🟡 PMC7543720
   *Why:* the 3-variable equation is reproduced in §1 above. Implement it exactly. It is the reference assignment for everything downstream, and its coefficients are the thing you are testing under missingness.

2. **Sinha, Churpek, Calfee (2020).** Machine Learning Classifier Models Can Identify ARDS Phenotypes Using Readily Available Clinical Data. *AJRCCM* 202(7):996–1004. 🔴 (abstract + the accompanying editorial, which is more accessible)
   *Why:* the biomarker-free GBM. Note the variable list and the admission-time constraint. **Note also:** the authors state a commitment to sharing the model with investigators seeking to validate it. Emailing Pratik Sinha at UCSF to request it is legitimate and costs one email — worth doing on day 0, and a request that goes unanswered costs nothing.

3. **Calfee, Delucchi, Parsons et al. (2014).** Subphenotypes in ARDS: latent class analysis of data from two RCTs. *Lancet Respir Med* 2:611–620. 🔴
   *Why:* the origin. Read for the LCA variable set and the hyper/hypo characterisation.

4. **PhysioNet/CinC Challenge 2019 description page + Reyna et al., *Crit Care Med* 2020;48(2):210–217.** 🟢 The manuscript PDF is on the dataset page.
   *Why:* your primary dataset. Read the sepsis labelling definition carefully — `SepsisLabel` is shifted 6h early, which affects any use as an endpoint.

5. **PHIND (2026).** Bedside identification of subphenotypes in acute respiratory failure. *Lancet Respir Med.* 🟢 CC BY 4.0
   *Why:* current state of the art, open access, and the source of the standardization critique. Read what they concluded about feasibility.

### Tier 2 — read before designing the missingness analysis

This is where the project's actual contribution lives, so read these properly rather than skimming.

6. **Informative missingness in ICU EHR data.** Two accessible entries: the JMIR study on missingness indicators as predictive features in MIMIC-III (🟡 PMC6329436), and the medRxiv paper on sampling-frequency bias and missingness patterns in critical care data (🟢 medRxiv 2024.06.09.24308661).
   *Why:* establishes that missingness encodes care process. Your §4.4 mask channel is grounded in this literature — cite it rather than inventing the justification.

7. **Comparative evaluation of handling missing data and modalities in EHR** (*Int J Med Inform*, 2026). 🔴 abstract accessible
   *Why:* the central methodological finding you will lean on — models that learn directly from incomplete data hold up as missingness increases, while imputation can shift the distribution and degrade separation. This is the design rationale for mask-channel models over imputation.

8. **Sinha et al. (2022).** Validation and utility of ARDS subphenotypes identified by machine-learning models. *Lancet Respir Med.* 🔴
   *Why:* shows what performance looks like in *curated* observational cohorts. Your numbers on sparse routine data should be compared against these. The delta is your result.

### Tier 3 — Layer 1 and the secondary arm

9. **Scicluna, van Vught, Zwinderman et al. (2017).** Classification of patients with sepsis according to blood genomic endotype (MARS). *Lancet Respir Med* 5:816–826. 🔴
   *Why:* the source of GSE65682's endotype labels. Check the supplement for the GSM-level label table (§1.5 scoping question).

10. **Filbin, Mehta, Schneider et al. (2021).** Longitudinal proteomic analysis of severe COVID-19. *Cell Rep Med.* 🟡
    *Why:* the Olink cohort. Read the acuity definitions before using them as outcomes.

11. **Christie, Vaslef, Chang et al. (2015).** A Randomized Dose-Escalation Study of the Safety and Anti-Inflammatory Activity of the p38 MAPK Inhibitor Dilmapimod in Severe Trauma Subjects at Risk for ARDS. *Crit Care Med* 43(9):1859–1869. 🔴
    *Why:* the only p38-in-ARDS clinical trial. Essential context for the README's framing, and short enough to read from the abstract if paywalled.

12. **A p38/MK2 mechanism review** covering MK2 → tristetraprolin → TNF mRNA stability.
    *Why:* grounds the §3.6 caveat about why transcript abundance cannot measure this pathway. Do not write the S5 section without it.

### Tier 4 — background, skim only

13. GEn1E's own site and any published material on GEn-1124 / the GRID platform.
14. Reviews of ARDS precision medicine and treatable traits.
15. Challenge 2019 methods papers (the arXiv analysis of the 2019 Challenge algorithms) — useful for baseline expectations on this dataset.

---

## 4. Instructions for Claude Code

- **Read Tier 1 before writing Layer 2 code.** Items 1, 4, and 5 are open access; fetch and read them properly. For paywalled items, read the abstract and the secondary descriptions — do not fabricate methods detail you have not seen.
- **Implement the Sinha 3-variable equation exactly as printed in §1.** Do not refit it, do not "improve" the coefficients, do not restandardize the inputs. It is a fixed external reference, and the whole design depends on it being external.
- **`HCO3` is the only one of its three variables present in Challenge 2019.** IL-6 and sTNFR1 are not in EHR data. This is the point, not a problem: it is exactly why the biomarker-free GBM (item 2) exists, and exactly the substitution a real deployment must make.
- **Where a paper's model is unavailable, say so in the README** rather than approximating it silently and calling it "the Sinha model."
- **Cite every claim of prior work.** A README that asserts "two subphenotypes have been identified" without a citation reads as though the project discovered it.
- Do not let reading expand indefinitely. Tier 1 is a day. Tier 2 is half a day. Tiers 3–4 are reference material to consult when writing the relevant section, not to read front-to-back.

---

## 5. What this changes in the spec

If the narrowed thesis in §1 is adopted, `ards-endotype-project-spec.md` needs:

- **§0.1** — thesis restated around missingness robustness rather than portability
- **§1** — P1 becomes "reproduce the published classifiers as a validated baseline" (supporting, one figure). New headline: **M1** degradation curve under progressive variable ablation; **M2** performance on natively-sparse vs natively-complete patients; **M3** cross-hospital shift decomposed into case-mix vs measurement-practice components
- **§1.5** — scoping gains a step: quantify how often the Sinha variable sets are fully observed in Challenge 2019. If the answer is "almost never," that single number is arguably the project's headline result
- **§4.5** — a2/a3 reframed as ablation experiments rather than portability experiments
- **§6.2** — sparsity is upgraded from a risk to the instrument
- **§6.13** — derivative risk sharpens: the project must be explicit that it is *stress-testing* published models, not claiming to have invented them
