# Day 0 — data acquisition

**You run these, not Claude Code** (SPEC.md §6.12: all downloads are manual and local).
Total footprint ≈ 2–4 GB. Nothing here needs credentialing, a CITI course, or a DUA.

Everything lands under `data/`, which is gitignored. Set `ARDS_DATA` if you want it on
another disk; every script reads paths from `shared/paths.py` and nothing hardcodes.

```bash
cd /home/james/ards-endotypes
mkdir -p data
```

Run `python layer0_scoping/01_scope_all_sources.py --check-only` at any point — it reports
which sources are present, which are missing, and their observed sizes, without analysing
anything. That is the acceptance test for this page.

---

## 1. PhysioNet/CinC Challenge 2019 — the primary cohort (42 MB)

40,336 ICU patients, one `.psv` per patient, one row per hour, 40 variables.
Open Access, CC BY 4.0. Sets A and B are **two different hospital systems**, which is the
free cross-institution validation in §4.1.

```bash
cd data
wget -r -N -c -np -nH --cut-dirs=3 \
     https://physionet.org/files/challenge-2019/1.0.0/
mv 1.0.0 challenge2019 2>/dev/null || mv challenge-2019/1.0.0 challenge2019
```

The archive ships `training_setA.zip` / `training_setB.zip` (or already-extracted
`training/` and `training_setB/` directories, depending on mirror). Extract whatever
arrives as zips:

```bash
cd challenge2019 && for z in *.zip; do unzip -q -n "$z"; done
```

**Expected after extraction:** two directories of `p??????.psv` files totalling 40,336.
Do not rename them — the scoping script discovers the layout rather than assuming it
(§6.12), and reports what it finds.

---

## 2. eICU-CRD Demo — the mortality cohort (~50 MB)

~2,500 unit stays across ~20 hospitals. This is the **only** source of a mortality
outcome in Layer 2 (§4.4, §6.1); Challenge 2019 has none.

```bash
cd data
wget -r -N -c -np -nH --cut-dirs=3 \
     https://physionet.org/files/eicu-crd-demo/2.0.1/
mv 2.0.1 eicu_demo 2>/dev/null || mv eicu-crd-demo/2.0.1 eicu_demo
```

Tables needed: `patient`, `vitalPeriodic`, `lab`, `apachePatientResult`,
`respiratoryCare`, `treatment`. Files are `.csv.gz`; leave them compressed, DuckDB reads
them directly.

---

## 3. MIMIC-IV Clinical Database Demo — smoke tests only (~100 MB)

100 patients. Used for schema familiarity and pipeline smoke tests, nothing else.
**This is not a licence to build Layer 2 on MIMIC** — see SPEC.md §6.12 and DECISIONS.md D4.

```bash
cd data
wget -r -N -c -np -nH --cut-dirs=3 \
     https://physionet.org/files/mimic-iv-demo/2.2/
mv 2.2 mimic_iv_demo 2>/dev/null || mv mimic-iv-demo/2.2 mimic_iv_demo
```

---

## 4. GEO transcriptomics — Layer 1 (~1–2 GB)

Three cohorts, three platforms. **Never merge and cluster them** (§3.5, §6.5).

```bash
mkdir -p data/geo && cd data/geo
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/GSE65682_series_matrix.txt.gz
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32707/matrix/GSE32707_series_matrix.txt.gz
wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66890/matrix/GSE66890_series_matrix.txt.gz
```

Platform annotation tables (for probe → gene symbol collapse) download on first use via
GEOparse into `data/geo/`. If you prefer to pre-fetch:

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13667/annot/GPL13667.annot.gz   # GSE65682
wget -c https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10558/annot/GPL10558.annot.gz   # GSE32707
wget -c https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6244/annot/GPL6244.annot.gz     # GSE66890
```

> **Hour-1 question, before any Layer 1 code** (§3.1, §6.6): are the MARS endotype labels
> actually in GSE65682's `pData`? The scoping script prints every characteristics field
> and its unique values. If they are absent, the supervised anchor (§3.3) is dropped and
> Layer 1 runs unsupervised-only — a planned fallback, not a crisis. Timebox the search
> for them at 2 hours, then take the fallback.

---

## 5. MGH COVID Olink Explore 1536 — Layer 1 proteomics (25.8 MB) ✅ DOWNLOADED

**The URL in the spec (`olink.com/mgh-covid-study`) is dead — it 404s.** Olink removed the
landing page but never removed the file. There is **no registration and no DUA**; the old
gate was a marketing form that no longer exists.

```bash
cd data/olink_mgh
wget -O MGH_Olink_COVID_Apr_27_2021.zip \
  "https://f.hubspotusercontent40.net/hubfs/7074596/MGH-Broad%20Covid%20study/MGH_Olink_COVID_Apr_27_2021.zip"
md5sum MGH_Olink_COVID_Apr_27_2021.zip   # c13826c4f9c3ad0501649008abb811f5
unzip -q MGH_Olink_COVID_Apr_27_2021.zip && rm -rf __MACOSX
```

The paper's citable deposit is Mendeley `10.17632/nf853r8xsj` — **but you must use version
`/2`**. The Key Resources Table cites `/1`, whose NPX table is a 255-protein subset that
does not include TNFRSF1A. Prefer the ZIP regardless: it is the April 2021 revision with
per-plate LOD (1,472 assays vs Mendeley's 1,429) and it ships long-format with gene names,
where Mendeley is wide-format with bare OlinkID headers.

**Verified contents** (parsed directly, not assumed):

| File | Shape |
|---|---|
| `MGH_COVID_OLINK_NPX.txt` | 1,271,808 × 16, semicolon-delimited, long format |
| `MGH_COVID_Clinical_Info.txt` | 384 × 44, semicolon-delimited |
| `Variable_descriptions_2021-04-27.xlsx` | data dictionary |

Join on `subject_id`. 792 SampleIDs, of which **8 are QC/plate controls with a blank
`Timepoint` — filter these out first**. Timepoints are `D0` (383), `D3` (219), `D7` (139),
`DE` (41, an event draw such as ICU admission), plus two singleton `DE.1`/`DE.2`.
**132 subjects have all three of D0/D3/D7; 226 have at least two.** One subject (id 180)
appears in the clinical file with no proteomics, so 383 subjects carry NPX.

> ⚠️ **Two version traps.**
> 1. This release uses **`Acuity_0/3/7/28/max` on a 1–5 scale** (1 = death, 5 = discharged).
>    The earlier release used **WHO 0/3/7/28 on a 1–6 scale**. A data dictionary describing
>    WHO columns belongs to the *old* version. Do not mix the two encodings; take the
>    clinical file from this ZIP.
> 2. **TNFRSF1A is on the NEUROLOGY panel (OID21155)**, not Inflammation. Subsetting to the
>    Inflammation panel silently deletes it. IL6 appears on all four panels
>    (OID20101/20563/20911/21276) — pick one and say which.
>
> Match proteins on **UniProt**, never on assay name: this protein is `TNFRSF1A` here,
> `TNF-R1` in Olink Target 96, and `TNF sR-I` in SomaScan.

**There is no bicarbonate, lactate, pH or any blood gas in the clinical file** — see
DECISIONS.md D9. Beware greps for "lactate", which false-match "Lactate dehydrogenase"
(`ldh_*_cat`). LDH is not lactate.

---

## 6. Two things worth doing on day 0 that are not downloads

**a. Submit the PhysioNet credentialing application** (SPEC.md §10). It takes minutes,
its review clock runs independently, and it costs nothing to have the option. If it lands
mid-project, **do not** rebuild Layer 2 on MIMIC-IV — a half-migrated Layer 2 on day 20 is
the worst of both paths. It is a post-deadline stretch only.

**b. Email Pratik Sinha at UCSF** (PRIOR_ART.md item 2). The authors state a commitment to
sharing the biomarker-free GBM with investigators seeking to validate it. One email; an
unanswered request costs nothing. If it arrives, §4.5 gains a second published external
comparator, which is the single cheapest strengthening available to this project.

---

## 7. Tier 1 reading — day 1, before any Layer 2 code

PRIOR_ART.md §3 is the full list. These four are open access and are the ones that gate code:

| # | Paper | Access | Gates |
|---|---|---|---|
| 1 | Sinha et al. 2020, *Lancet Respir Med* 8:247–257 | PMC7543720 | `shared/sinha2020.py` VERIFY-1/2/3 |
| 4 | Reyna et al. 2020, *Crit Care Med* 48(2):210–217 | on the dataset page | `SepsisLabel`'s 6h shift — it is the primary endpoint |
| 5 | PHIND 2026, *Lancet Respir Med* | CC BY 4.0 | current state of the art; the standardisation critique |
| 6 | Informative missingness: PMC6329436; medRxiv 2024.06.09.24308661 | open | the justification for mask channels (§4.4) — cite it, don't invent it |

**All four VERIFY items in `shared/sinha2020.py` were closed on day 0** against the paper's
supplementary appendix, so paper 1 no longer gates Layer 2 code. What was found:

- Units are pg/mL for both cytokines, mmol/L for bicarbonate (Table S2).
- The equation is a linear predictor; the logistic transform is required.
- Two published cut-points, 0.274 (Youden) and 0.5, both now exposed.
- The coefficients come from appendix **Table S4** and are an **ancillary** model, not the
  paper's headline 3-variable model (which uses IL-8 and protein C). Wording matters here.
- The sTNFR1 coefficient in circulation (1.1913) is a typo for **1.1903**.

Still worth reading paper 1 properly for the LCA variable set and the derivation design —
just not as a blocker. Paper 4 (Reyna) remains a genuine gate: `SepsisLabel` is the primary
endpoint and its 6-hour shift changes how it can be used.
