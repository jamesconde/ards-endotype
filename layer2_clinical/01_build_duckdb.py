#!/usr/bin/env python3
"""Load Challenge 2019 and the eICU demo into DuckDB. SPEC.md §4.2.

The scale problem the spec worried about does not exist on this path: the whole
Challenge corpus is ~320 MB of text and fits in RAM many times over. DuckDB is
used because it reads the pipe-delimited files directly and gives the later
SQL steps somewhere to live, not because the data is large.

Every stage is checkpointed (§8.3). A silent zero-row load here would poison
everything downstream.

Usage
-----
    python layer2_clinical/01_build_duckdb.py            # build
    python layer2_clinical/01_build_duckdb.py --rebuild  # drop and rebuild
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, step

log = setup("L2_build_duckdb")

EXPECTED_PATIENTS = 40336
EXPECTED_A = 20336
EXPECTED_B = 20000


def build_hourly(con: duckdb.DuckDBPyConnection) -> None:
    """One row per patient-hour, with hospital_system derived from the path.

    The directory layout is discovered rather than assumed (§6.12): the glob is
    verified to match files before the load runs.
    """
    root = paths.CHALLENGE
    globs = {
        "A": str(root / "training" / "training_setA" / "*.psv"),
        "B": str(root / "training" / "training_setB" / "*.psv"),
    }
    for name, g in globs.items():
        n = len(list(Path(g).parent.glob("*.psv")))
        log.info("set %s: %d files at %s", name, n, g)
        if n == 0:
            raise CheckpointError(
                f"set {name}: no .psv files at {g}. The layout differs from the "
                f"expected training/training_set{name}/ — inspect data/challenge2019 "
                f"and fix the glob rather than guessing."
            )

    # No window function here. ROW_NUMBER() OVER (PARTITION BY filename) forces a
    # full materialisation of 1.5M rows and OOMs on a 16 GB laptop — and it is
    # redundant, because ICULOS already numbers the hours 1..n within a patient.
    # That is asserted below rather than assumed.
    con.execute("DROP TABLE IF EXISTS hourly")
    for setname, g in globs.items():
        con.execute(f"""
            {'CREATE TABLE hourly AS' if setname == 'A' else 'INSERT INTO hourly BY NAME'}
            SELECT
                CAST(regexp_extract(filename, 'p(\\d+)\\.psv', 1) AS INTEGER) AS patient_id,
                '{setname}' AS hospital_system,
                * EXCLUDE (filename)
            FROM read_csv_auto('{g}', delim='|', nullstr='NaN',
                               filename=true, header=true)
        """)
        log.info("  loaded set %s", setname)

    # ICULOS is the hour index everywhere downstream, so verify its structure
    # before anything depends on it.
    #
    # It is CONTIGUOUS (no gaps, no duplicates) for every patient, but it does
    # NOT always start at 1: 8,793 of 40,336 patients begin at ICULOS 2 or later,
    # because their record starts some hours after ICU admission. ICULOS is
    # therefore hours-since-ICU-admission, not row position, and the two differ
    # for 22% of patients. Anything that indexes hours by row order is wrong for
    # those patients — which is exactly the kind of off-by-a-variable-amount error
    # that produces plausible output and no error message.
    bad = con.execute("""
        SELECT COUNT(*) FROM (
            SELECT patient_id
            FROM hourly
            GROUP BY patient_id
            HAVING (MAX(ICULOS) - MIN(ICULOS) + 1) <> COUNT(*)
                OR COUNT(DISTINCT ICULOS) <> COUNT(*)
        )""").fetchone()[0]
    if bad:
        raise CheckpointError(
            f"{bad} patients have a non-contiguous or duplicated ICULOS. The hour "
            "index cannot be taken from ICULOS; add an explicit row number per patient.")

    late = con.execute(
        "SELECT COUNT(*) FROM (SELECT patient_id FROM hourly GROUP BY 1 HAVING MIN(ICULOS) > 1)"
    ).fetchone()[0]
    log.info("✓ ICULOS contiguous for all patients; %d (%.1f%%) start after hour 1 — "
             "use ICULOS, never row position", late, 100 * late / EXPECTED_PATIENTS)
    report_only("patients_starting_after_hour_1", late,
                "ICULOS != row position for these (D14)")

    n_rows, n_pat = con.execute(
        "SELECT COUNT(*), COUNT(DISTINCT patient_id) FROM hourly").fetchone()
    by_set = con.execute(
        "SELECT hospital_system, COUNT(DISTINCT patient_id) "
        "FROM hourly GROUP BY 1 ORDER BY 1").fetchall()
    log.info("hourly: %s rows, %s patients, split %s", f"{n_rows:,}", f"{n_pat:,}", by_set)

    split = dict(by_set)
    if n_pat != EXPECTED_PATIENTS:
        raise CheckpointError(
            f"expected {EXPECTED_PATIENTS} distinct patients, got {n_pat}. "
            "Do not proceed — the load is incomplete or patient_id parsing is wrong.")
    if split.get("A") != EXPECTED_A or split.get("B") != EXPECTED_B:
        raise CheckpointError(f"expected A={EXPECTED_A} B={EXPECTED_B}, got {split}")

    # patient_id must be unique across the two systems, or every later join is wrong.
    overlap = con.execute("""
        SELECT COUNT(*) FROM (
            SELECT patient_id FROM hourly GROUP BY 1 HAVING COUNT(DISTINCT hospital_system) > 1
        )""").fetchone()[0]
    if overlap:
        raise CheckpointError(
            f"{overlap} patient_ids appear in BOTH hospital systems. The id is not a "
            "global key; joins would silently merge patients from different hospitals.")

    report_only("duckdb_hourly_rows", n_rows)
    report_only("duckdb_patients", n_pat)


def build_eicu(con: duckdb.DuckDBPyConnection) -> None:
    """eICU demo tables needed for the mortality replication (§4.4) and C2 (D13)."""
    root = paths.EICU_DEMO
    wanted = {
        "eicu_patient": "patient.csv.gz",
        "eicu_apache": "apachePatientResult.csv.gz",
        "eicu_lab": "lab.csv.gz",
        "eicu_vital": "vitalPeriodic.csv.gz",
        "eicu_resp": "respiratoryCharting.csv.gz",
        "eicu_apachevar": "apacheApsVar.csv.gz",
    }
    for table, fname in wanted.items():
        f = root / fname
        if not f.exists():
            log.warning("eICU %s missing (%s) — skipping", table, fname)
            continue
        con.execute(f"DROP TABLE IF EXISTS {table}")
        con.execute(f"CREATE TABLE {table} AS SELECT * FROM read_csv_auto('{f}', "
                    f"header=true, union_by_name=true, sample_size=-1)")
        n = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        log.info("  %-16s %8s rows  (%s)", table, f"{n:,}", fname)
        if n == 0:
            raise CheckpointError(f"{table}: loaded 0 rows from {fname}")

    n_stay = con.execute(
        "SELECT COUNT(DISTINCT patientunitstayid) FROM eicu_patient").fetchone()[0]
    n_hosp = con.execute("SELECT COUNT(DISTINCT hospitalid) FROM eicu_patient").fetchone()[0]
    log.info("eICU: %d unit stays across %d hospitals", n_stay, n_hosp)
    if not (2000 <= n_stay <= 3000):
        raise CheckpointError(f"eICU demo: expected ~2,500 stays, got {n_stay}")

    # The external published model whose predictions anchor C2 (D13).
    apache = con.execute("""
        SELECT apacheversion, COUNT(*) n,
               SUM(CASE WHEN predictedhospitalmortality IS NOT NULL THEN 1 ELSE 0 END) with_pred
        FROM eicu_apache GROUP BY 1 ORDER BY 1""").fetchall()
    log.info("APACHE predictions available: %s", apache)

    mort = con.execute("""
        SELECT hospitaldischargestatus, COUNT(*) FROM eicu_patient GROUP BY 1 ORDER BY 2 DESC
        """).fetchall()
    log.info("eICU hospital discharge status: %s", mort)
    died = dict((str(k), v) for k, v in mort).get("Expired", 0)
    report_only("eicu_deaths", died, "the only mortality outcome in Layer 2 (§6.1)")
    if died < 50:
        raise CheckpointError(
            f"only {died} deaths in the eICU demo — too few for a stable mortality "
            "estimate. §6.1 says leave the mortality question open rather than "
            "manufacture a substitute.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rebuild", action="store_true", help="drop the database first")
    args = ap.parse_args()

    if args.rebuild and paths.DUCKDB.exists():
        paths.DUCKDB.unlink()
        log.warning("dropped existing %s", paths.DUCKDB)

    paths.DUCKDB.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(paths.DUCKDB))
    # Leave headroom on a 16 GB laptop and give DuckDB somewhere to spill rather
    # than dying. The scratch dir is gitignored (§0.2).
    con.execute("SET memory_limit='6GB'")
    con.execute(f"SET temp_directory='{paths.SCRATCH}'")
    con.execute("SET preserve_insertion_order=false")

    with step(log, "load Challenge 2019 hourly") as r:
        build_hourly(con)
        r["n_out"] = con.execute("SELECT COUNT(*) FROM hourly").fetchone()[0]

    with step(log, "load eICU demo") as r:
        build_eicu(con)
        r["n_out"] = con.execute("SELECT COUNT(*) FROM eicu_patient").fetchone()[0]

    tables = con.execute("SHOW TABLES").fetchall()
    log.info("database %s contains: %s", paths.DUCKDB, [t[0] for t in tables])
    con.close()
    log.info("✓ build complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
