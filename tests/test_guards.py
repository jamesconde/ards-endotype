"""Tests for the machinery that is supposed to stop this project going wrong quietly."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier

from shared import sinha2020
from shared.assertions import CheckpointError, attrition_step, checkpoint
from shared.no_population_stats import (
    PortabilityError,
    assert_batch_invariant,
    assert_no_population_stats,
)


# --- checkpoints -------------------------------------------------------------
def test_checkpoint_passes_and_returns_frame():
    df = pd.DataFrame({"a": [1, 2, 3]})
    assert checkpoint(df, "t_ok", 1, 10, required_cols=["a"]) is df


def test_checkpoint_catches_the_silent_zero_row_join():
    """The single most common way this kind of project dies quietly (§8.3)."""
    with pytest.raises(CheckpointError, match="got 0 rows"):
        checkpoint(pd.DataFrame({"a": []}), "t_empty", 1, 10)


def test_checkpoint_catches_missing_column_and_nulls_and_dupes():
    df = pd.DataFrame({"a": [1, 1, None]})
    with pytest.raises(CheckpointError, match="missing columns"):
        checkpoint(df, "t", 1, 10, required_cols=["b"])
    with pytest.raises(CheckpointError, match="nulls in"):
        checkpoint(df, "t", 1, 10, no_nulls=["a"])
    with pytest.raises(CheckpointError, match="duplicate"):
        checkpoint(pd.DataFrame({"id": [1, 1]}), "t", 1, 10, unique_on="id")


def test_attrition_trail_records_exclusions():
    rows: list[dict] = []
    attrition_step(rows, "all ICU patients", 40336)
    attrition_step(rows, "age >= 18", 39000, "paediatric")
    attrition_step(rows, "ICULOS >= 24h", 21000, "short stay")
    assert [r["excluded"] for r in rows] == [None, 1336, 18000]


# --- portability -------------------------------------------------------------
@pytest.fixture
def xy():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 4)) * [10, 1, 100, 0.1]
    y = (X[:, 0] + X[:, 2] / 10 > 0).astype(int)
    return X, y


def test_scaler_in_pipeline_is_rejected(xy):
    """The exact failure §6.12 predicts: a fitted StandardScaler called portable."""
    X, y = xy
    pipe = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression())]).fit(X, y)
    with pytest.raises(PortabilityError, match="StandardScaler"):
        assert_no_population_stats(pipe, "scaled_lr")


def test_raw_value_models_pass(xy):
    X, y = xy
    for name, model in [("lr", LogisticRegression(max_iter=2000)),
                        ("tree", DecisionTreeClassifier(max_depth=3)),
                        ("gbm", GradientBoostingClassifier(n_estimators=10))]:
        assert_no_population_stats(model.fit(X, y), name)


def test_imputer_counts_as_a_population_statistic(xy):
    """Cohort-median imputation is a population statistic, not a convenience."""
    from sklearn.impute import SimpleImputer
    X, y = xy
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                     ("lr", LogisticRegression())]).fit(X, y)
    with pytest.raises(PortabilityError, match="SimpleImputer|imputation"):
        assert_no_population_stats(pipe, "imputed_lr")


def test_batch_invariance_catches_what_introspection_cannot(xy):
    """A custom predict that z-scores inside itself passes introspection."""
    X, y = xy
    lr = LogisticRegression(max_iter=2000).fit(X, y)

    def sneaky(Z):                      # standardises against whoever is in the batch
        Z = np.asarray(Z, dtype=float)
        return lr.predict_proba((Z - Z.mean(0)) / (Z.std(0) + 1e-9))[:, 1]

    def honest(Z):
        return lr.predict_proba(np.asarray(Z, dtype=float))[:, 1]

    assert_no_population_stats(lr, "lr")            # introspection says fine
    assert_batch_invariant(honest, X, label="honest")
    with pytest.raises(PortabilityError, match="depend on the batch"):
        assert_batch_invariant(sneaky, X, label="sneaky")


# --- the external reference model -------------------------------------------
def test_sinha_self_test():
    sinha2020.self_test()


def test_sinha_rejects_unit_errors():
    """sTNFR1 in ng/mL instead of pg/mL must fail loudly, not silently reclassify."""
    with pytest.raises(ValueError, match="unit error"):
        sinha2020.classify(il6_pg_ml=[100.0], hco3_mmol_l=[22.0], stnfr1_pg_ml=[5.0e6])
    with pytest.raises(ValueError, match="unit error"):
        sinha2020.classify(il6_pg_ml=[100.0], hco3_mmol_l=[2200.0], stnfr1_pg_ml=[5000.0])


def test_sinha_never_imputes_a_missing_input():
    r = sinha2020.classify([np.nan, 1500.0], [22.0, 18.0], [5000.0, 15000.0])
    assert r.n_incomplete == 1
    assert r.label[0] == "unclassifiable"
    assert r.label[1] == "hyperinflammatory"


def test_sinha_is_itself_batch_invariant():
    """The reference model must satisfy the same portability standard it anchors."""
    rng = np.random.default_rng(1)
    X = np.column_stack([rng.uniform(10, 5000, 100),
                         rng.uniform(12, 34, 100),
                         rng.uniform(500, 40000, 100)])

    def predict(Z):
        return sinha2020.classify(Z[:, 0], Z[:, 1], Z[:, 2]).probability

    assert_batch_invariant(predict, X, label="sinha2020")


def test_rank_surrogate_is_not_sold_as_the_classifier():
    """NPX must not reach classify(); the surrogate exists so it cannot (D3)."""
    df = pd.DataFrame({"IL6": [1.2, 4.5, 3.0], "HCO3": [22, 18, 25],
                       "TNFRSF1A": [2.1, 5.5, 3.3]})
    s = sinha2020.rank_surrogate(df, "IL6", "HCO3", "TNFRSF1A")
    assert len(s) == 3 and s.notna().all()
    # NPX values are far below the pg/mL plausible range and must be refused.
    with pytest.raises(ValueError, match="unit error"):
        sinha2020.classify(df.IL6, df.HCO3, df.TNFRSF1A)
