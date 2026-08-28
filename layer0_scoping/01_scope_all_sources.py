#!/usr/bin/env python3
"""Day 1 — descriptive scoping across every source. SPEC.md §1.5.

DESCRIPTIVE ONLY. No models, no clustering, no hypothesis tests, no filtering
toward a nicer answer. This script exists to find out what the data *is* before
the analysis plan is frozen, because the alternative is discovering on day 15
that the cohort does not exist.

Writes results/L0_scoping.md. Hard stop at one day (§7): scoping that runs long
becomes the analysis, which is the failure it exists to prevent.

Usage
-----
    python layer0_scoping/01_scope_all_sources.py --check-only   # presence only
    python layer0_scoping/01_scope_all_sources.py                # full scoping
    python layer0_scoping/01_scope_all_sources.py --only challenge

Every source is optional. A missing source is reported as missing, not raised —
the script must run usefully on day 0 with nothing downloaded.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import report_only
from shared.logging_setup import setup

log = setup("L0_scoping")

# Columns whose coverage the analysis plan actually depends on (§1.5).
CRITICAL_COLS = ["FiO2", "O2Sat", "SaO2", "HCO3", "Lactate", "Creatinine",
                 "Bilirubin_total", "Platelets", "WBC", "pH", "BaseExcess"]

SINHA_VARS = {"HCO3": "present in Challenge 2019",
              "IL6": "NOT in EHR data — this is the point, not a problem (PRIOR_ART §4)",
              "sTNFR1": "NOT in EHR data — ditto"}


# ---------------------------------------------------------------------------
# presence check
# ---------------------------------------------------------------------------
def dir_size_mb(p: Path) -> float:
    if not p.exists():
        return 0.0
    if p.is_file():
        return p.stat().st_size / 1e6
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) / 1e6


def check_sources() -> pd.DataFrame:
    rows = []
    for name, path, expect in [
        ("challenge2019", paths.CHALLENGE, "~42 MB, 40,336 .psv files"),
        ("eicu_demo", paths.EICU_DEMO, "~50 MB, eICU relational tables"),
        ("mimic_iv_demo", paths.MIMIC_DEMO, "~100 MB, smoke tests only"),
        ("geo", paths.GEO, "3 series matrices (.txt.gz)"),
        ("olink_mgh", paths.OLINK, "NPX + clinical, manual download"),
    ]:
        present = path.exists()
        n_files = len(list(path.rglob("*"))) if present else 0
        rows.append({"source": name, "present": present, "path": str(path),
                     "size_mb": round(dir_size_mb(path), 1),
                     "n_files": n_files, "expected": expect})
    df = pd.DataFrame(rows)
    log.info("\n%s", df[["source", "present", "size_mb", "n_files"]].to_string(index=False))
    missing = df.loc[~df.present, "source"].tolist()
    if missing:
        log.warning("missing: %s — see DOWNLOADS.md", ", ".join(missing))
    return df


# ---------------------------------------------------------------------------
# Challenge 2019
# ---------------------------------------------------------------------------
def find_psv_files(root: Path) -> dict[str, list[Path]]:
    """Discover the layout rather than assuming it (§6.12, §4.2)."""
    files = sorted(root.rglob("*.psv"))
    by_set: dict[str, list[Path]] = {"A": [], "B": [], "unknown": []}
    for f in files:
        parts = "/".join(f.parts).lower()
        if "seta" in parts or "training_set_a" in parts or parts.endswith("/training/" + f.name.lower()):
            by_set["A"].append(f)
        elif "setb" in parts or "training_set_b" in parts:
            by_set["B"].append(f)
        else:
            by_set["unknown"].append(f)
    # If the split is undetectable, fall back to the top-level directory name.
    if not by_set["A"] and not by_set["B"] and files:
        groups: dict[str, list[Path]] = {}
        for f in files:
            groups.setdefault(f.parent.name, []).append(f)
        log.warning("Could not detect setA/setB from paths. Top-level groups found: %s",
                    {k: len(v) for k, v in groups.items()})
        by_set["unknown"] = files
    return by_set


def scope_challenge(sample_n: int | None = None,
                    sets: tuple[str, ...] | None = None) -> dict:
    """Scope the Challenge 2019 corpus.

    `sets` restricts to named hospital systems. Use it when a set is still
    downloading: a partial set would silently bias the A/B coverage comparison,
    and that comparison decides whether §4.5d measures transfer or missingness.
    """
    root = paths.CHALLENGE
    if not root.exists():
        return {"present": False}

    by_set = find_psv_files(root)
    if sets:
        dropped = {k: len(v) for k, v in by_set.items() if k not in sets and v}
        by_set = {k: v for k, v in by_set.items() if k in sets}
        if dropped:
            log.warning("restricted to sets %s; EXCLUDED %s. The A/B comparison "
                        "is not computed in this run.", sets, dropped)
    counts = {k: len(v) for k, v in by_set.items()}
    total = sum(counts.values())
    log.info("Challenge 2019 layout: %s (total %d)", counts, total)
    report_only("challenge2019_n_patients", total)

    if total == 0:
        return {"present": True, "error": "no .psv files found — extract the zips (DOWNLOADS.md)"}

    if sets is None and total != 40336:
        log.warning("expected 40,336 patient files, found %d — check the extraction", total)
    elif sets:
        log.info("scoping %d patients from set(s) %s", total, ",".join(sets))

    files = [(f, s) for s, fl in by_set.items() for f in fl]
    if sample_n and sample_n < len(files):
        rng = np.random.default_rng(42)
        idx = rng.choice(len(files), size=sample_n, replace=False)
        files = [files[i] for i in idx]
        log.warning("SAMPLING %d of %d patients — scoping only, never for analysis", sample_n, total)

    # Accumulate per-patient coverage without holding every row in memory.
    per_patient = []
    all_cols: set[str] = set()
    hour_counts, col_nonnull_hours, col_total_hours = [], {}, 0

    for i, (f, hosp) in enumerate(files):
        if i % 5000 == 0 and i:
            log.info("  ... %d/%d patients read", i, len(files))
        try:
            df = pd.read_csv(f, sep="|", na_values=["NaN", "nan", ""])
        except Exception as e:  # noqa: BLE001 — report and continue, one bad file must not kill scoping
            log.error("unreadable: %s (%s)", f.name, e)
            continue
        all_cols |= set(df.columns)
        rec = {"patient": f.stem, "hospital_system": hosp, "n_hours": len(df)}
        for c in df.columns:
            rec[f"has_{c}"] = bool(df[c].notna().any())
            col_nonnull_hours[c] = col_nonnull_hours.get(c, 0) + int(df[c].notna().sum())
        col_total_hours += len(df)
        for c in ("Age", "Gender", "ICULOS"):
            if c in df.columns:
                rec[c] = df[c].dropna().iloc[-1] if df[c].notna().any() else np.nan
        if "SepsisLabel" in df.columns:
            rec["sepsis_any"] = bool(df["SepsisLabel"].fillna(0).max() > 0)
        hour_counts.append(len(df))
        per_patient.append(rec)

    pp = pd.DataFrame(per_patient)
    has_cols = [c for c in pp.columns if c.startswith("has_")]

    patient_coverage = (pp[has_cols].mean()
                        .rename(lambda c: c[4:])
                        .sort_values(ascending=False))
    hour_density = pd.Series(
        {c: col_nonnull_hours.get(c, 0) / col_total_hours for c in sorted(all_cols)}
    ).sort_values(ascending=False)

    fio2_cov = float(patient_coverage.get("FiO2", np.nan))
    report_only("challenge2019_fio2_patient_coverage", round(fio2_cov, 4),
                "fraction of patients with >=1 non-null FiO2 (§6.2)")

    # Selection check (§6.2): do FiO2-recorded patients differ from the rest?
    selection = {}
    if "has_FiO2" in pp.columns:
        for c in ("Age", "ICULOS", "n_hours", "sepsis_any"):
            if c in pp.columns:
                a = pd.to_numeric(pp.loc[pp.has_FiO2, c], errors="coerce").dropna()
                b = pd.to_numeric(pp.loc[~pp.has_FiO2, c], errors="coerce").dropna()
                if len(a) and len(b):
                    pooled = np.sqrt((a.var() + b.var()) / 2) or np.nan
                    selection[c] = {"with_fio2_mean": round(float(a.mean()), 3),
                                    "without_fio2_mean": round(float(b.mean()), 3),
                                    "std_diff": round(float((a.mean() - b.mean()) / pooled), 3)
                                    if pooled and np.isfinite(pooled) else None}

    # A vs B: transfer test or missingness test? (§1.5, §4.5d — a crucial distinction)
    ab = {}
    if {"A", "B"} <= set(pp.hospital_system.unique()):
        for c in has_cols:
            a = pp.loc[pp.hospital_system == "A", c].mean()
            b = pp.loc[pp.hospital_system == "B", c].mean()
            ab[c[4:]] = {"A": round(float(a), 4), "B": round(float(b), 4),
                         "abs_diff": round(abs(float(a - b)), 4)}
        worst = sorted(ab.items(), key=lambda kv: -kv[1]["abs_diff"])[:10]
        log.info("largest A/B coverage gaps: %s",
                 {k: v["abs_diff"] for k, v in worst})

    los = pd.to_numeric(pp.get("n_hours"), errors="coerce")
    out = {
        "present": True,
        "n_patients": total,
        "n_scoped": len(pp),
        "set_counts": counts,
        "columns": sorted(all_cols),
        "patient_coverage": patient_coverage.round(4).to_dict(),
        "hour_density": hour_density.round(4).to_dict(),
        "critical_coverage": {c: round(float(patient_coverage.get(c, np.nan)), 4)
                              for c in CRITICAL_COLS},
        "iculos": {"median": float(los.median()), "p25": float(los.quantile(.25)),
                   "p75": float(los.quantile(.75)), "max": float(los.max()),
                   "frac_ge_24h": float((los >= 24).mean())},
        "sepsis_rate": float(pp["sepsis_any"].mean()) if "sepsis_any" in pp else None,
        "fio2_selection_check": selection,
        "ab_coverage_gaps": ab,
    }

    # The number PRIOR_ART.md §5 argues could be the project's headline result.
    if {"FiO2", "HCO3"} <= set(pp.columns.str.replace("has_", "", regex=False)):
        joint = float((pp.get("has_FiO2", False) & pp.get("has_HCO3", False)).mean())
        out["joint_fio2_hco3_coverage"] = round(joint, 4)
        report_only("challenge2019_joint_fio2_hco3", round(joint, 4),
                    "fraction with BOTH — cohort definition and Sinha comparator")
    return out


# ---------------------------------------------------------------------------
# GEO
# ---------------------------------------------------------------------------
def scope_geo() -> dict:
    root = paths.GEO
    if not root.exists():
        return {"present": False}
    out: dict = {"present": True, "series": {}}

    for gse in ("GSE65682", "GSE32707", "GSE66890"):
        matches = list(root.glob(f"{gse}*series_matrix.txt*"))
        if not matches:
            out["series"][gse] = {"present": False}
            continue
        path = matches[0]
        info = {"present": True, "file": path.name}
        try:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", errors="replace") as fh:
                sample_lines, header_lines, n_data = [], [], 0
                in_table = False
                for line in fh:
                    if line.startswith("!series_matrix_table_begin"):
                        in_table = True
                        continue
                    if line.startswith("!series_matrix_table_end"):
                        break
                    if in_table:
                        n_data += 1
                        if n_data <= 3:
                            sample_lines.append(line.rstrip("\n")[:400])
                    elif line.startswith("!Sample_"):
                        header_lines.append(line.rstrip("\n"))
            info["n_probes_approx"] = max(n_data - 1, 0)

            # pData characteristics — the MARS-label question (§3.1, §6.6).
            chars = {}
            for h in header_lines:
                key = h.split("\t", 1)[0].lstrip("!")
                if "characteristics" in key or key in ("!Sample_title", "Sample_title"):
                    vals = [v.strip('"') for v in h.split("\t")[1:]]
                    uniq = sorted({v for v in vals if v})[:12]
                    chars.setdefault(key, []).append(
                        {"n_unique": len(set(vals)), "examples": uniq})
            info["characteristics"] = chars
            info["n_samples"] = max(
                (len(h.split("\t")) - 1 for h in header_lines), default=0)

            if gse == "GSE65682":
                blob = " ".join(header_lines).lower()
                hits = [t for t in ("mars1", "mars2", "mars3", "mars4", "endotype",
                                    "mars_endotype", "class") if t in blob]
                info["mars_label_search"] = {
                    "tokens_found": hits,
                    "verdict": ("MARS labels appear present — confirm the exact field"
                                if any(h.startswith("mars") for h in hits)
                                else "MARS labels NOT visible in pData. Fallback: extract "
                                     "from Scicluna 2017 supplement, join by GSM. "
                                     "Timebox 2h then run unsupervised-only (§6.6)."),
                }
                log.warning("GSE65682 MARS labels: %s", info["mars_label_search"]["verdict"])
        except Exception as e:  # noqa: BLE001
            info["error"] = str(e)
            log.error("%s: %s", gse, e)
        out["series"][gse] = info
    return out


# ---------------------------------------------------------------------------
# Olink
# ---------------------------------------------------------------------------
def _read_any(path: Path, nrows: int | None = None) -> pd.DataFrame:
    for sep in (None, ",", "\t", ";"):
        try:
            return pd.read_csv(path, sep=sep, engine="python", nrows=nrows)
        except Exception:  # noqa: BLE001
            continue
    raise ValueError(f"could not parse {path}")


def scope_olink() -> dict:
    root = paths.OLINK
    if not root.exists():
        return {"present": False}
    files = [f for f in root.rglob("*") if f.is_file()
             and f.suffix.lower() in {".csv", ".txt", ".tsv", ".xlsx", ".gz"}]
    if not files:
        return {"present": True, "error": "no data files found — manual download (DOWNLOADS.md §5)"}

    out: dict = {"present": True, "files": {}}
    observed_schema = {}
    for f in files:
        entry: dict = {"size_mb": round(f.stat().st_size / 1e6, 2)}
        try:
            head = _read_any(f, nrows=2000)
            entry["columns"] = list(head.columns)
            entry["n_cols"] = len(head.columns)
            observed_schema[f.name] = list(head.columns)

            # Detect long-format NPX without assuming column names (§3.1).
            lower = {c.lower(): c for c in head.columns}
            assay_col = next((lower[k] for k in ("assay", "olinkid", "protein", "uniprot")
                              if k in lower), None)
            npx_col = next((lower[k] for k in ("npx", "value", "normalized_npx")
                            if k in lower), None)
            entry["looks_like_long_npx"] = bool(assay_col and npx_col)
            if assay_col:
                full = _read_any(f)
                assays = full[assay_col].astype(str).unique()
                entry["n_assays"] = int(len(assays))
                # The two proteins carrying the Sinha comparator (§3.7).
                for target in ("IL6", "IL-6", "TNFRSF1A", "TNFR1"):
                    entry[f"has_{target}"] = bool(
                        any(target.lower() == a.lower() for a in assays))
                lod_col = next((lower[k] for k in ("lod", "limit_of_detection", "missingfreq")
                                if k in lower), None)
                entry["lod_column"] = lod_col
                if lod_col and npx_col:
                    below = full[npx_col] < full[lod_col]
                    entry["frac_below_lod_overall"] = round(float(below.mean()), 4)
                time_col = next((lower[k] for k in ("timepoint", "day", "visit", "time")
                                 if k in lower), None)
                entry["timepoint_column"] = time_col
                if time_col:
                    entry["timepoints"] = sorted(
                        str(v) for v in full[time_col].dropna().unique()[:12])
        except Exception as e:  # noqa: BLE001
            entry["error"] = str(e)
        out["files"][f.name] = entry

    schema_path = paths.CONFIG / "olink_schema_observed.yaml"
    try:
        import yaml
        schema_path.write_text(yaml.safe_dump(observed_schema, sort_keys=False))
        log.info("observed Olink schema -> %s (every later script reads this, "
                 "no column name is assumed)", schema_path)
    except ImportError:
        schema_path.with_suffix(".json").write_text(json.dumps(observed_schema, indent=2))
    return out


# ---------------------------------------------------------------------------
# eICU demo
# ---------------------------------------------------------------------------
def scope_eicu() -> dict:
    root = paths.EICU_DEMO
    if not root.exists():
        return {"present": False}
    needed = ["patient", "vitalPeriodic", "lab", "apachePatientResult",
              "respiratoryCare", "treatment"]
    out: dict = {"present": True, "tables": {}}
    for t in needed:
        matches = [f for f in root.rglob("*") if f.is_file()
                   and f.name.lower().startswith(t.lower())]
        if not matches:
            out["tables"][t] = {"present": False}
            continue
        f = matches[0]
        entry = {"present": True, "file": f.name, "size_mb": round(f.stat().st_size / 1e6, 2)}
        try:
            head = pd.read_csv(f, nrows=5000)
            entry["columns"] = list(head.columns)
            if t == "patient":
                full = pd.read_csv(f)
                entry["n_rows"] = len(full)
                for c in ("hospitaldischargestatus", "unitdischargestatus", "hospitalid"):
                    if c in full.columns:
                        entry[f"{c}_nonnull_frac"] = round(float(full[c].notna().mean()), 4)
                        if "status" in c:
                            entry[f"{c}_values"] = full[c].value_counts().head(6).to_dict()
                        if c == "hospitalid":
                            entry["n_hospitals"] = int(full[c].nunique())
        except Exception as e:  # noqa: BLE001
            entry["error"] = str(e)
        out["tables"][t] = entry
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def write_report(findings: dict, presence: pd.DataFrame) -> Path:
    L: list[str] = []
    A = L.append
    A("# L0 — descriptive scoping\n")
    A("*Generated by `layer0_scoping/01_scope_all_sources.py`. Descriptive only: "
      "no models, no tests, no filtering. SPEC.md §1.5.*\n")
    A("## Sources present\n")
    A(presence[["source", "present", "size_mb", "n_files", "expected"]]
      .to_markdown(index=False))
    A("")

    ch = findings.get("challenge", {})
    A("## Challenge 2019\n")
    if not ch.get("present"):
        A("Not downloaded. See DOWNLOADS.md §1.\n")
    elif ch.get("error"):
        A(f"**Error:** {ch['error']}\n")
    else:
        A(f"- Patient files: **{ch['n_patients']:,}** "
          f"(expected 40,336); split {ch['set_counts']}")
        A(f"- ICULOS: median {ch['iculos']['median']:.0f}h, "
          f"IQR {ch['iculos']['p25']:.0f}–{ch['iculos']['p75']:.0f}h; "
          f"**{ch['iculos']['frac_ge_24h']:.1%} have ≥24h of record**")
        if ch.get("sepsis_rate") is not None:
            A(f"- Sepsis (any hour): {ch['sepsis_rate']:.1%}")
        A("\n### Coverage of the columns the plan depends on\n")
        A("| Variable | Patients with ≥1 value | Hour-level density |")
        A("|---|---|---|")
        for c in CRITICAL_COLS:
            pc = ch["critical_coverage"].get(c)
            hd = ch["hour_density"].get(c)
            A(f"| `{c}` | {pc:.1%} | {hd:.2%} |" if pc is not None and hd is not None
              else f"| `{c}` | — | — |")
        A("")
        fio2 = ch["critical_coverage"].get("FiO2")
        if fio2 is not None:
            A(f"**FiO2 coverage is {fio2:.1%}.** The cohort definition (§4.3) is built on "
              "it, so requiring a recorded FiO2 selects for oxygen-supported and more "
              "heavily monitored patients. If this is low, §4.3 needs rewriting before "
              "it is written, and Table 1 must show the selection.\n")
        if ch.get("joint_fio2_hco3_coverage") is not None:
            A(f"**Joint FiO2 + HCO3 coverage: {ch['joint_fio2_hco3_coverage']:.1%}.** "
              "PRIOR_ART.md §5 argues that how often the published classifiers' variable "
              "sets are *fully observed* may be the project's headline result. This is "
              "that number, for the two variables this cohort has.\n")
        if ch.get("fio2_selection_check"):
            A("### Selection induced by requiring FiO2 (§6.2)\n")
            A("| Variable | With FiO2 | Without FiO2 | Std. diff |")
            A("|---|---|---|---|")
            for k, v in ch["fio2_selection_check"].items():
                A(f"| {k} | {v['with_fio2_mean']} | {v['without_fio2_mean']} | {v['std_diff']} |")
            A("")
        if ch.get("ab_coverage_gaps"):
            gaps = sorted(ch["ab_coverage_gaps"].items(),
                          key=lambda kv: -kv[1]["abs_diff"])[:12]
            A("### Set A vs set B coverage — *transfer test or missingness test?*\n")
            A("This determines what §4.5d is actually measuring. If A and B differ "
              "substantially in coverage, a performance drop on B is a missingness "
              "artifact, not a transfer failure, and reporting it as transfer would "
              "invalidate the experiment.\n")
            A("| Variable | A | B | |Δ| |")
            A("|---|---|---|---|")
            for k, v in gaps:
                A(f"| `{k}` | {v['A']:.1%} | {v['B']:.1%} | {v['abs_diff']:.1%} |")
            biggest = gaps[0][1]["abs_diff"] if gaps else 0
            A(f"\nLargest coverage gap: **{biggest:.1%}**. "
              + ("Small — A→B is plausibly a genuine case-mix transfer test."
                 if biggest < 0.10 else
                 "Large — A→B confounds transfer with measurement practice. §4.5d must "
                 "decompose the two, and the README must say which is being measured.")
              + "\n")
        A("### Sinha 2020 variables in this cohort\n")
        for var, note in SINHA_VARS.items():
            A(f"- `{var}` — {note}")
        A("")

    for key, title, section in [("geo", "GEO transcriptomics", "§4"),
                                ("olink", "MGH Olink", "§5"),
                                ("eicu", "eICU demo", "§2")]:
        f = findings.get(key, {})
        A(f"## {title}\n")
        if not f.get("present"):
            A(f"Not downloaded. See DOWNLOADS.md {section}.\n")
            continue
        A("```json")
        A(json.dumps(f, indent=2, default=str)[:6000])
        A("```\n")

    A("---\n")
    A("## Decisions this scoping forces, before anything is frozen\n")
    A("1. Does the §4.3 cohort definition survive the observed FiO2 coverage?")
    A("2. Is A→B a transfer test or a missingness test? (§4.5d)")
    A("3. Are the MARS labels in GSE65682's pData, or is Layer 1 unsupervised-only? (§6.6)")
    A("4. Do IL6 and TNFRSF1A survive the Olink LOD filter? (§3.7)")
    A("5. Does the eICU demo have enough qualifying patients for a stable mortality "
      "estimate — or does the mortality question stay open? (§6.1)\n")
    A("**Then, and only then:** finalise the question wording, write and commit "
      "`config/genesets.yaml`, and start. Descriptive statistics before the freeze, "
      "inferential claims after, never the reverse.\n")

    out = paths.RESULTS / "L0_scoping.md"
    out.write_text("\n".join(L))
    (paths.RESULTS / "L0_scoping.json").write_text(
        json.dumps(findings, indent=2, default=str))
    log.info("wrote %s", out)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-only", action="store_true",
                    help="report which sources are present, then stop")
    ap.add_argument("--only", choices=["challenge", "geo", "olink", "eicu"],
                    help="scope a single source")
    ap.add_argument("--sample", type=int, default=None,
                    help="scope only N Challenge patients (development only)")
    ap.add_argument("--sets", default=None,
                    help="restrict Challenge scoping to hospital systems, e.g. --sets A "
                         "(use while a set is still downloading)")
    args = ap.parse_args()

    presence = check_sources()
    if args.check_only:
        n_missing = int((~presence.present).sum())
        print(f"\n{len(presence) - n_missing}/{len(presence)} sources present.")
        if n_missing:
            print("See DOWNLOADS.md — downloads are manual (SPEC.md §6.12).")
        return 0

    findings: dict = {}
    want = {args.only} if args.only else {"challenge", "geo", "olink", "eicu"}
    if "challenge" in want:
        findings["challenge"] = scope_challenge(
            sample_n=args.sample,
            sets=tuple(args.sets.split(",")) if args.sets else None)
    if "geo" in want:
        findings["geo"] = scope_geo()
    if "olink" in want:
        findings["olink"] = scope_olink()
    if "eicu" in want:
        findings["eicu"] = scope_eicu()

    write_report(findings, presence)
    log.info("Scoping complete. Hard stop at one day (§7) — do not let this become "
             "the analysis.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
