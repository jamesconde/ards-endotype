#!/usr/bin/env python3
"""Day 1 — does the §4.3 cohort actually exist? Descriptive only.

SPEC.md §1.5 says the FiO2 coverage number determines whether §4.3 needs
rewriting before it is written. §6.2 says to sanity-plot the S/F distribution
*before* filtering on it. This script does both, and reports the yield at the
primary threshold and at the pre-authorised relaxation.

It answers one question: if we apply the cohort definition as specified, how
many patients come out, and are they a defensible population?

NO MODELS. This is a counting exercise whose output is an attrition table.

Usage
-----
    python layer0_scoping/02_sf_feasibility.py --sets A
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import attrition_step, report_only
from shared.logging_setup import setup

log = setup("L0_sf_feasibility")

sys.path.insert(0, str(paths.REPO / "layer0_scoping"))
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_scope", paths.REPO / "layer0_scoping" / "01_scope_all_sources.py")
_scope = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_scope)


def load_params() -> dict:
    return yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sets", default=None, help="restrict to hospital systems, e.g. A")
    args = ap.parse_args()

    p = load_params()
    inc, exc, rel = p["inclusion"], p["exclusion"], p["relaxation"]
    fio2v = p["fio2_validity"]
    n_dropped = [0]   # FiO2 hours rejected as physiologically impossible (D10)

    by_set = _scope.find_psv_files(paths.CHALLENGE)
    if args.sets:
        keep = tuple(args.sets.split(","))
        by_set = {k: v for k, v in by_set.items() if k in keep}
        log.warning("restricted to sets %s — yields are for this subset only", keep)
    files = [(f, s) for s, fl in by_set.items() for f in fl]
    log.info("examining %d patients", len(files))

    fio2_values: list[float] = []
    rows: list[dict] = []

    for i, (f, hosp) in enumerate(files):
        if i % 5000 == 0 and i:
            log.info("  ... %d/%d", i, len(files))
        df = pd.read_csv(f, sep="|", na_values=["NaN", "nan", ""])
        age = df["Age"].dropna()
        rec = {
            "patient": f.stem, "hospital_system": hosp,
            "age": float(age.iloc[-1]) if len(age) else np.nan,
            "n_hours": len(df),
            "sepsis_any": bool(df["SepsisLabel"].fillna(0).max() > 0),
        }
        fio2_values.extend(df["FiO2"].dropna().tolist())

        # Hours where BOTH the numerator and denominator of S/F exist AND the
        # FiO2 is physiologically possible (D10). An out-of-range FiO2 is not a
        # measurement, and dividing by it manufactures a qualifying hour.
        valid_fio2 = df["FiO2"].between(fio2v["min"], fio2v["max"])
        n_dropped[0] += int((df["FiO2"].notna() & ~valid_fio2).sum())
        both = valid_fio2 & df["O2Sat"].notna()
        rec["n_hours_both"] = int(both.sum())
        rec["has_both_hour"] = bool(both.any())

        if both.any():
            sub = df.loc[both].copy()
            sub["hour"] = sub.index + 1
            # Exclusion applied at the candidate hour, not to the patient (§4.3).
            ok = (sub["O2Sat"] < exc["o2sat_saturation_cutoff"]) & \
                 (sub["hour"] <= exc["max_qualifying_hour"])
            sub = sub.loc[ok]
            if len(sub):
                sf = sub["O2Sat"] / sub["FiO2"]
                rec["min_sf"] = float(sf.min())
                rec["first_qualifying_hour_235"] = (
                    int(sub.loc[sf <= inc["sf_ratio_threshold"], "hour"].min())
                    if (sf <= inc["sf_ratio_threshold"]).any() else None)
                rec["qualifies_235"] = bool((sf <= inc["sf_ratio_threshold"]).any())
                rec["qualifies_315"] = bool((sf <= rel["relaxed_sf_threshold"]).any())
        rows.append(rec)

    d = pd.DataFrame(rows)
    for c in ("qualifies_235", "qualifies_315"):
        d[c] = d.get(c).astype("boolean").fillna(False).astype(bool)
    d["sepsis_any"] = d["sepsis_any"].astype(bool)

    # --- FiO2 scale check: fraction (0.21-1.0) or percent (21-100)? ----------
    fv = pd.Series(fio2_values)
    log.info("FiO2 observed: min=%.3f p50=%.3f max=%.3f n=%d",
             fv.min(), fv.median(), fv.max(), len(fv))
    frac_like = float((fv <= 1.0).mean())
    log.info("fraction of FiO2 values <= 1.0: %.1f%%", frac_like * 100)
    if 0.01 < frac_like < 0.99:
        log.error("FiO2 IS MIXED SCALE — %.1f%% look like fractions and the rest like "
                  "percentages. S/F is meaningless until this is resolved.", frac_like * 100)
    scale = "fraction" if frac_like > 0.99 else ("percent" if frac_like < 0.01 else "MIXED")
    log.info("FiO2 scale: %s", scale)
    report_only("fio2_scale_fraction_like", round(frac_like, 4), f"-> {scale}")
    log.warning("FiO2 hours rejected as out-of-range [%.2f, %.2f]: %d of %d (%.3f%%)",
                fio2v["min"], fio2v["max"], n_dropped[0], len(fv),
                100 * n_dropped[0] / max(len(fv), 1))
    report_only("fio2_hours_rejected_invalid", n_dropped[0], "D10 validity window")

    # --- attrition ----------------------------------------------------------
    at: list[dict] = []
    attrition_step(at, "set A patients", len(d))
    d1 = d[d.age >= inc["min_age"]]
    attrition_step(at, f"age >= {inc['min_age']}", len(d1), "paediatric or missing age")
    d2 = d1[d1.n_hours >= inc["min_iculos_hours"]]
    attrition_step(at, f"ICULOS >= {inc['min_iculos_hours']}h", len(d2), "short record")
    d3 = d2[d2.has_both_hour]
    attrition_step(at, "has >=1 hour with FiO2 AND O2Sat", len(d3), "S/F not computable")
    d4 = d3[d3.min_sf.notna()]
    attrition_step(at, "after O2Sat<97 and hour<=168 at candidate hours", len(d4),
                   "oximeter saturated or too late")
    d5 = d4[d4.qualifies_235]
    attrition_step(at, f"S/F <= {inc['sf_ratio_threshold']} (primary)", len(d5))

    at_df = pd.DataFrame(at)
    out_csv = paths.RESULTS / "L2_attrition_setA_feasibility.csv"
    at_df.to_csv(out_csv, index=False)
    log.info("\n%s", at_df.to_string(index=False))

    n_235, n_315 = len(d5), int(d4.qualifies_315.sum())
    report_only("cohort_yield_sf235_setA", n_235)
    report_only("cohort_yield_sf315_setA", n_315)

    verdict = (f"PRIMARY THRESHOLD VIABLE ({n_235} patients)"
               if n_235 >= rel["min_cohort_n"] else
               f"BELOW {rel['min_cohort_n']} — take the pre-authorised relaxation to "
               f"S/F <= {rel['relaxed_sf_threshold']} ({n_315} patients)")
    log.warning("VERDICT: %s", verdict)

    # --- who gets selected, and how different are they? ---------------------
    sel = {}
    d_elig = d[(d.age >= inc["min_age"]) & (d.n_hours >= inc["min_iculos_hours"])]
    inc_mask = d_elig.index.isin(d5.index)
    for c in ("age", "n_hours", "sepsis_any"):
        a = pd.to_numeric(d_elig.loc[inc_mask, c], errors="coerce").dropna()
        b = pd.to_numeric(d_elig.loc[~inc_mask, c], errors="coerce").dropna()
        if len(a) and len(b):
            pooled = np.sqrt((a.var() + b.var()) / 2)
            sel[c] = {"cohort": round(float(a.mean()), 3),
                      "excluded": round(float(b.mean()), 3),
                      "std_diff": round(float((a.mean() - b.mean()) / pooled), 3)
                      if pooled else None}

    sf_desc = d4.min_sf.describe(percentiles=[.05, .25, .5, .75, .95]).round(1).to_dict()

    lines = [
        "# L0b — does the §4.3 cohort exist?\n",
        "*Descriptive feasibility count on Challenge 2019"
        f"{' set ' + args.sets if args.sets else ''}. No models. SPEC.md §1.5, §4.3, §6.2.*\n",
        f"## Verdict\n\n**{verdict}**\n",
        f"## FiO2 scale and validity\n\nObserved min {fv.min():.3f}, median "
        f"{fv.median():.3f}, max {fv.max():.3f} over {len(fv):,} recorded values. "
        f"{frac_like:.1%} are <= 1.0, so the column is recorded as a **{scale}**. "
        "A mixed scale would make the ratio meaningless; this one is not mixed.\n\n"
        f"But the raw range is impossible: room air is 0.209 and pure oxygen is 1.0. "
        f"**{n_dropped[0]:,} values ({100 * n_dropped[0] / max(len(fv), 1):.2f}%) fall "
        f"outside [{fio2v['min']}, {fio2v['max']}] and are rejected as S/F candidates** "
        "(DECISIONS.md D10). They include zeros, which send S/F to infinity, and values "
        "above 1.0, which *deflate* S/F and would falsely admit a patient as severely "
        "hypoxemic. The patient is not excluded — only the impossible hour is.\n",
        "## Attrition\n", at_df.to_markdown(index=False), "",
        f"\nAt the relaxed threshold (S/F <= {rel['relaxed_sf_threshold']}, "
        f"a surrogate for P/F <= 300): **{n_315}** patients.\n",
        "## Minimum S/F per patient, among those where it is computable\n",
        "| statistic | value |", "|---|---|",
        *[f"| {k} | {v} |" for k, v in sf_desc.items()], "",
        "\n## Who the cohort selects (vs. eligible but excluded)\n",
        "| Variable | In cohort | Excluded | Std. diff |", "|---|---|---|---|",
        *[f"| {k} | {v['cohort']} | {v['excluded']} | {v['std_diff']} |"
          for k, v in sel.items()], "",
        "\nStandardised differences above 0.1 indicate imbalance. This table belongs in "
        "Table 1 and in the README: the cohort is a selected subpopulation, not a random "
        "sample of ICU patients, and the selection is induced by the measurement itself "
        "(§6.2).\n",
    ]
    out = paths.RESULTS / "L0b_cohort_feasibility.md"
    out.write_text("\n".join(lines))
    log.info("wrote %s and %s", out, out_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
