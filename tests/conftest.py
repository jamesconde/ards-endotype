"""Synthetic fixtures shaped like the real sources.

These let the pipeline be exercised on day 0, before anything is downloaded, and
they are the only data that ever lives in this repo. They are deliberately
*structurally* faithful (column names, pipe delimiter, NaN encoding, hourly
rows, informative missingness) and biologically meaningless — they must never
be mistaken for, or substituted into, an analysis.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

CHALLENGE_COLS = [
    "HR", "O2Sat", "Temp", "SBP", "MAP", "DBP", "Resp", "EtCO2",
    "BaseExcess", "HCO3", "FiO2", "pH", "PaCO2", "SaO2", "AST", "BUN",
    "Alkalinephos", "Calcium", "Chloride", "Creatinine", "Bilirubin_direct",
    "Glucose", "Lactate", "Magnesium", "Phosphate", "Potassium",
    "Bilirubin_total", "TroponinI", "Hct", "Hgb", "PTT", "WBC", "Fibrinogen",
    "Platelets", "Age", "Gender", "Unit1", "Unit2", "HospAdmTime", "ICULOS",
    "SepsisLabel",
]

# Roughly the real sparsity ordering: vitals dense, labs sparse, FiO2 sparser.
HOURLY_PRESENCE = {c: 0.95 for c in ["HR", "O2Sat", "SBP", "MAP", "DBP", "Resp"]}
HOURLY_PRESENCE.update({c: 0.05 for c in ["HCO3", "pH", "BaseExcess", "PaCO2", "Lactate"]})
HOURLY_PRESENCE.update({c: 0.02 for c in ["FiO2", "SaO2", "Fibrinogen", "TroponinI"]})


def _make_patient(rng: np.random.Generator, n_hours: int, sick: bool) -> pd.DataFrame:
    df = pd.DataFrame(index=range(n_hours))
    for c in CHALLENGE_COLS:
        df[c] = np.nan
    df["HR"] = rng.normal(95 if sick else 78, 12, n_hours)
    df["O2Sat"] = np.clip(rng.normal(93 if sick else 97, 3, n_hours), 60, 100)
    df["MAP"] = rng.normal(68 if sick else 82, 10, n_hours)
    df["SBP"] = df["MAP"] + rng.normal(35, 8, n_hours)
    df["DBP"] = df["MAP"] - rng.normal(12, 5, n_hours)
    df["Resp"] = rng.normal(24 if sick else 17, 4, n_hours)
    df["HCO3"] = rng.normal(19 if sick else 25, 3, n_hours)
    df["Lactate"] = np.abs(rng.normal(3.5 if sick else 1.2, 1.0, n_hours))
    df["FiO2"] = rng.choice([0.21, 0.4, 0.6, 0.8, 1.0], n_hours)
    df["pH"] = rng.normal(7.31 if sick else 7.40, 0.06, n_hours)
    df["Creatinine"] = np.abs(rng.normal(1.9 if sick else 0.9, 0.5, n_hours))
    df["Platelets"] = np.abs(rng.normal(140 if sick else 240, 60, n_hours))
    df["WBC"] = np.abs(rng.normal(16 if sick else 8, 4, n_hours))
    df["Bilirubin_total"] = np.abs(rng.normal(1.8 if sick else 0.7, 0.6, n_hours))

    # Informative missingness: sicker patients get drawn more often (§4.4).
    for c, base in HOURLY_PRESENCE.items():
        p = min(base * (2.2 if sick else 1.0), 0.98)
        df.loc[rng.random(n_hours) > p, c] = np.nan

    df["Age"] = float(rng.integers(18, 92))
    df["Gender"] = float(rng.integers(0, 2))
    df["Unit1"] = float(rng.integers(0, 2))
    df["Unit2"] = 1.0 - df["Unit1"]
    df["HospAdmTime"] = float(-rng.exponential(20))
    df["ICULOS"] = np.arange(1, n_hours + 1)
    df["SepsisLabel"] = 0
    if sick and n_hours > 12 and rng.random() < 0.6:
        df.loc[rng.integers(8, n_hours):, "SepsisLabel"] = 1
    return df[CHALLENGE_COLS]


@pytest.fixture
def fake_challenge(tmp_path):
    """A miniature Challenge 2019 tree: two hospital systems, real file layout.

    Set B is given systematically worse lab coverage so the A-vs-B section of
    the scoping report has something to detect — the exact confound §4.5d has
    to distinguish from genuine transfer failure.
    """
    rng = np.random.default_rng(42)
    root = tmp_path / "challenge2019"
    for setname, n_pat, cov_mult in [("training_setA", 40, 1.0), ("training_setB", 30, 0.4)]:
        d = root / setname
        d.mkdir(parents=True)
        for i in range(n_pat):
            sick = rng.random() < 0.4
            df = _make_patient(rng, int(rng.integers(8, 100)), sick)
            if cov_mult < 1.0:
                for c in ("FiO2", "HCO3", "Lactate"):
                    drop = rng.random(len(df)) > cov_mult
                    df.loc[drop, c] = np.nan
            df.to_csv(d / f"p{100000 + i + (0 if setname.endswith('A') else 500000)}.psv",
                      sep="|", index=False, na_rep="NaN")
    return root
