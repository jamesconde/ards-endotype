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

## 5. MGH COVID Olink Explore 1536 — Layer 1 proteomics (~200 MB)

**Manual, requires free registration.** No wget line exists for this one.

1. Go to <https://www.olink.com/mgh-covid-study/>
2. Register (free) and download the NPX data file and the clinical/metadata file
3. Put both in `data/olink_mgh/` under their original filenames

Ships as a long-format NPX file plus a separate clinical file. **Do not assume column
names** (§3.1) — the scoping script prints both schemas and writes the observed columns
into `config/olink_schema_observed.yaml`, which every later script reads.

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

The three VERIFY items in `shared/sinha2020.py` are open questions until paper 1 is read.
They are the highest-value hour in the whole schedule: the module currently infers its
units from an internal consistency argument, and a consistency argument is not a citation.
