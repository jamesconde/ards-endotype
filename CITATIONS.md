# Citations

Attribution is a licence condition here, not a courtesy (SPEC.md §0.2 item 1).
Challenge 2019 is CC BY 4.0.

## Datasets

**PhysioNet/CinC Challenge 2019** — primary Layer 2 cohort, 40,336 ICU patients.
> Reyna MA, Josef CS, Jeter R, et al. Early Prediction of Sepsis From Clinical Data:
> The PhysioNet/Computing in Cardiology Challenge 2019. *Crit Care Med.*
> 2020;48(2):210–217. doi:10.1097/CCM.0000000000004145

**eICU Collaborative Research Database Demo v2.0.1** — mortality replication.
> Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU
> Collaborative Research Database, a freely available multi-center database for
> critical care research. *Sci Data.* 2018;5:180178. doi:10.1038/sdata.2018.178

**MIMIC-IV Clinical Database Demo v2.2** — schema familiarity and smoke tests only.
> Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic
> health record dataset. *Sci Data.* 2023;10:1. doi:10.1038/s41597-022-01899-x

**PhysioNet platform** — required alongside every PhysioNet dataset citation.
> Goldberger AL, Amaral LAN, Glass L, et al. PhysioBank, PhysioToolkit, and PhysioNet:
> Components of a New Research Resource for Complex Physiologic Signals.
> *Circulation.* 2000;101(23):e215–e220.

**GSE65682** — MARS cohort, Layer 1 discovery.
> Scicluna BP, van Vught LA, Zwinderman AH, et al. Classification of patients with
> sepsis according to blood genomic endotype: a prospective cohort study.
> *Lancet Respir Med.* 2017;5(10):816–826. doi:10.1016/S2213-2600(17)30294-1

**GSE32707** — Layer 1 transfer cohort 1.
> Dolinay T, Kim YS, Howrylak J, et al. Inflammasome-regulated cytokines are critical
> mediators of acute lung injury. *Am J Respir Crit Care Med.* 2012;185(11):1225–1234.

**GSE66890** — Layer 1 transfer cohort 2.
> Kangelaris KN, Prakash A, Liu KD, et al. Increased expression of neutrophil-related
> genes in patients with early sepsis-induced ARDS. *Am J Physiol Lung Cell Mol Physiol.*
> 2015;308(11):L1102–L1113.

**MGH COVID-19 Olink Explore 1536** — Layer 1 proteomics.
> Filbin MR, Mehta A, Schneider AM, et al. Longitudinal proteomic analysis of severe
> COVID-19 reveals survival-associated signatures, tissue-specific cell death, and
> cell-cell interactions. *Cell Rep Med.* 2021;2(5):100287.

## Reference models reproduced, not invented

**The 3-variable parsimonious classifier** implemented in `shared/sinha2020.py`:

> ⚠️ **PROVISIONAL — do not reproduce this attribution elsewhere until it is confirmed.**
> The implemented equation uses **IL-6, sTNFR1 and bicarbonate**. At least two other
> three-variable classifiers circulate in this literature (IL-8/bicarbonate/protein C, and
> IL-8/bicarbonate/sTNFR1 as used in PHARLAP), and the paper below may publish a different
> variant than the coefficients in hand. The primary source for *these* coefficients is
> being verified. See DECISIONS.md D7.

> Sinha P, Delucchi KL, McAuley DF, et al. Development and validation of parsimonious
> algorithms to classify ARDS phenotypes. *Lancet Respir Med.* 2020;8(3):247–257.
> doi:10.1016/S2213-2600(19)30369-8

**The biomarker-free classifier** this project's own classifiers are benchmarked against:
> Sinha P, Churpek MM, Calfee CS. Machine Learning Classifier Models Can Identify ARDS
> Phenotypes Using Readily Available Clinical Data. *Am J Respir Crit Care Med.*
> 2020;202(7):996–1004. doi:10.1164/rccm.202002-0347OC

## Established findings this project reproduces rather than discovers

> Calfee CS, Delucchi K, Parsons PE, et al. Subphenotypes in acute respiratory distress
> syndrome: latent class analysis of data from two randomised controlled trials.
> *Lancet Respir Med.* 2014;2(8):611–620.

> Famous KR, Delucchi K, Ware LB, et al. ARDS Subphenotypes Respond Differently to
> Randomized Fluid Management Strategy. *Am J Respir Crit Care Med.* 2017;195(3):331–338.

> Sinha P, Delucchi KL, Chen Y, et al. Latent class analysis-derived subphenotypes are
> generalisable to observational cohorts of acute respiratory distress syndrome.
> *Thorax.* 2022;77(1):13–21.

## Clinical context for the target

> Christie JD, Vaslef S, Chang PK, et al. A Randomized Dose-Escalation Study of the
> Safety and Anti-Inflammatory Activity of the p38 MAPK Inhibitor Dilmapimod in Severe
> Trauma Subjects at Risk for Acute Respiratory Distress Syndrome.
> *Crit Care Med.* 2015;43(9):1859–1869.

---

*Entries are added as each source is actually used. A citation here for a dataset that
never made it into the analysis is as much an error as a missing one.*
