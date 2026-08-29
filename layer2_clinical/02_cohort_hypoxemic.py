#!/usr/bin/env python3
"""Define the hypoxemic respiratory failure cohort. SPEC.md §4.3.

NOT an ARDS cohort, and never called one (§6.3). There is no imaging, no PEEP
and no ventilation flag on this data path, so Berlin criteria cannot be
approximated even loosely. What is defined here is hypoxemic respiratory
failure identified by an SpO2/FiO2 surrogate.

Produces:
  cohort      -- one row per included patient, with T0
  results/L2_attrition.csv  -- CONSORT-style attrition, per hospital system

All thresholds come from config/cohort_params.yaml, frozen day 0 (I4).

Usage
-----
    python layer2_clinical/02_cohort_hypoxemic.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import duckdb
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, step

log = setup("L2_cohort")


def connect() -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(str(paths.DUCKDB))
    con.execute("SET memory_limit='6GB'")
    con.execute(f"SET temp_directory='{paths.SCRATCH}'")
    return con


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--relax", action="store_true",
                    help="use the pre-authorised relaxed S/F threshold (§4.3)")
    args = ap.parse_args()

    p = yaml.safe_load((paths.CONFIG / "cohort_params.yaml").read_text())
    inc, exc, rel = p["inclusion"], p["exclusion"], p["relaxation"]
    fio2v = p["fio2_validity"]
    sf_threshold = rel["relaxed_sf_threshold"] if args.relax else inc["sf_ratio_threshold"]
    log.info("S/F threshold: %d%s", sf_threshold, "  (RELAXED)" if args.relax else "")

    con = connect()

    with step(log, "build candidate hours"):
        # A candidate hour needs BOTH sides of the ratio, a physiologically
        # possible FiO2 (D10), an unsaturated oximeter, and to fall inside the
        # window. ICULOS is the hour, never row position (D14).
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE candidate_hours AS
            SELECT patient_id, hospital_system, ICULOS AS hour,
                   O2Sat, FiO2, O2Sat / FiO2 AS sf_ratio
            FROM hourly
            WHERE FiO2 BETWEEN {fio2v['min']} AND {fio2v['max']}
              AND O2Sat IS NOT NULL
              AND O2Sat < {exc['o2sat_saturation_cutoff']}
              AND ICULOS <= {exc['max_qualifying_hour']}
        """)
        n = con.execute("SELECT COUNT(*) FROM candidate_hours").fetchone()[0]
        log.info("candidate hours: %s", f"{n:,}")
        if n == 0:
            raise CheckpointError("no candidate hours — the filters are wrong")

    with step(log, "assemble cohort"):
        con.execute(f"""
            CREATE OR REPLACE TABLE cohort AS
            WITH patient_level AS (
                SELECT patient_id,
                       ANY_VALUE(hospital_system) AS hospital_system,
                       MAX(Age)                   AS age,
                       MAX(Gender)                AS gender,
                       MAX(HospAdmTime)           AS hosp_adm_time,
                       COUNT(*)                   AS n_hours,
                       MAX(ICULOS)                AS iculos_max,
                       MAX(SepsisLabel)           AS sepsis_any
                FROM hourly GROUP BY patient_id
            ),
            qualifying AS (
                SELECT patient_id, MIN(hour) AS t0, MIN(sf_ratio) AS min_sf
                FROM candidate_hours
                WHERE sf_ratio <= {sf_threshold}
                GROUP BY patient_id
            )
            SELECT pl.*, q.t0, q.min_sf
            FROM patient_level pl
            JOIN qualifying q USING (patient_id)
            WHERE pl.age >= {inc['min_age']}
              AND pl.n_hours >= {inc['min_iculos_hours']}
        """)

    # --- attrition, per hospital system --------------------------------------
    def count_at(where: str) -> pd.DataFrame:
        return con.execute(f"""
            SELECT hospital_system, COUNT(DISTINCT patient_id) AS n
            FROM ({where}) GROUP BY 1 ORDER BY 1""").df()

    steps = [
        ("all ICU patients", "SELECT patient_id, hospital_system FROM hourly"),
        (f"age >= {inc['min_age']}",
         f"SELECT patient_id, ANY_VALUE(hospital_system) hospital_system FROM hourly "
         f"GROUP BY patient_id HAVING MAX(Age) >= {inc['min_age']}"),
        (f"ICULOS >= {inc['min_iculos_hours']}h",
         f"SELECT patient_id, ANY_VALUE(hospital_system) hospital_system FROM hourly "
         f"GROUP BY patient_id HAVING MAX(Age) >= {inc['min_age']} "
         f"AND COUNT(*) >= {inc['min_iculos_hours']}"),
        (">=1 candidate hour (valid FiO2 + O2Sat<97, hour<=168)",
         "SELECT patient_id, hospital_system FROM candidate_hours"),
        (f"S/F <= {sf_threshold}  [COHORT]",
         "SELECT patient_id, hospital_system FROM cohort"),
    ]
    rows = []
    for label, q in steps:
        df = count_at(q).set_index("hospital_system")["n"]
        rows.append({"step": len(rows), "label": label,
                     "A": int(df.get("A", 0)), "B": int(df.get("B", 0)),
                     "total": int(df.sum())})
        log.info("attrition | %-56s A=%6d B=%6d total=%6d",
                 label, rows[-1]["A"], rows[-1]["B"], rows[-1]["total"])

    at = pd.DataFrame(rows)
    at["excluded"] = at["total"].shift(1) - at["total"]
    at.to_csv(paths.RESULTS / "L2_attrition.csv", index=False)

    n_total = int(at.iloc[-1]["total"])
    n_a, n_b = int(at.iloc[-1]["A"]), int(at.iloc[-1]["B"])
    report_only("cohort_n_total", n_total)
    report_only("cohort_n_setA", n_a)
    report_only("cohort_n_setB", n_b)

    if n_total < rel["min_cohort_n"] and not args.relax:
        raise CheckpointError(
            f"cohort is {n_total}, below the {rel['min_cohort_n']} floor. Re-run with "
            f"--relax to take the pre-authorised S/F <= {rel['relaxed_sf_threshold']} "
            "threshold, and record it in DECISIONS.md.")

    # Sanity: T0 must exist, be inside the window, and n must be plausible.
    bad_t0 = con.execute(
        f"SELECT COUNT(*) FROM cohort WHERE t0 IS NULL OR t0 < 1 "
        f"OR t0 > {exc['max_qualifying_hour']}").fetchone()[0]
    if bad_t0:
        raise CheckpointError(f"{bad_t0} cohort rows have an invalid T0")

    dupes = con.execute(
        "SELECT COUNT(*) FROM (SELECT patient_id FROM cohort GROUP BY 1 HAVING COUNT(*)>1)"
    ).fetchone()[0]
    if dupes:
        raise CheckpointError(f"{dupes} duplicate patient_id in cohort")

    log.info("\n%s", at.to_string(index=False))
    log.info("COHORT: %d patients (A=%d, B=%d)", n_total, n_a, n_b)
    log.warning("This cohort is HYPOXEMIC RESPIRATORY FAILURE, not ARDS (§6.3). "
                "The term must not appear in any figure, caption or variable name.")

    # Descriptive Table 1 input, no tests.
    t1 = con.execute("""
        SELECT hospital_system,
               COUNT(*) AS n,
               ROUND(AVG(age), 1) AS age_mean,
               ROUND(AVG(gender), 3) AS male_frac,
               ROUND(AVG(n_hours), 1) AS iculos_mean,
               ROUND(AVG(t0), 1) AS t0_mean,
               ROUND(AVG(min_sf), 1) AS min_sf_mean,
               ROUND(AVG(sepsis_any), 4) AS sepsis_rate
        FROM cohort GROUP BY 1 ORDER BY 1""").df()
    t1.to_csv(paths.RESULTS / "L2_table1_by_system.csv", index=False)
    log.info("\n%s", t1.to_string(index=False))

    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
