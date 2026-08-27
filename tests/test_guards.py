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


def _fake_npx(n: int = 200, seed: int = 3) -> pd.DataFrame:
    """NPX-like values: relative, log2, roughly 0-15. Not concentrations."""
    rng = np.random.default_rng(seed)
    hyper = rng.random(n) < 0.30
    return pd.DataFrame({
        "IL6": rng.normal(np.where(hyper, 8.0, 3.5), 1.0),
        "HCO3": rng.normal(np.where(hyper, 19.0, 25.0), 2.5),
        "TNFRSF1A": rng.normal(np.where(hyper, 6.5, 3.0), 1.0),
        "truth": np.where(hyper, "hyperinflammatory", "hypoinflammatory"),
    })


def test_npx_cannot_reach_the_published_equation():
    """The silent failure D3 exists to prevent: NPX sits inside the pg/mL range."""
    df = _fake_npx()
    with pytest.raises(ValueError, match="unit error"):
        sinha2020.classify(df.IL6, df.HCO3, df.TNFRSF1A)


def test_rank_surrogate_recovers_structure_without_claiming_probabilities():
    df = _fake_npx()
    labels = sinha2020.rank_surrogate_classes(df, "IL6", "HCO3", "TNFRSF1A")
    # Split at the assumed prevalence, so the class sizes are set by construction.
    assert abs((labels == "hyperinflammatory").mean() - 0.30) < 0.02
    # It should still track the underlying structure well above chance.
    agree = (labels == df.truth).mean()
    assert agree > 0.75, f"rank surrogate recovered only {agree:.2%} agreement"


def test_rank_surrogate_is_deliberately_not_portable():
    """It uses within-cohort ranks, so it must fail the portability check.

    Stated as a test so nobody later mistakes it for a deployable classifier.
    """
    df = _fake_npx()

    def predict(Z):
        d = pd.DataFrame(Z, columns=["IL6", "HCO3", "TNFRSF1A"])
        return sinha2020.rank_surrogate(d, "IL6", "HCO3", "TNFRSF1A").to_numpy()

    X = df[["IL6", "HCO3", "TNFRSF1A"]].to_numpy()
    with pytest.raises(PortabilityError, match="depend on the batch"):
        assert_batch_invariant(predict, X, label="rank_surrogate")


def test_model_variant_is_named_explicitly():
    """VERIFY-4: several parsimonious models circulate; the code must not be vague."""
    assert sinha2020.MODEL_VARIANT == "IL6_sTNFR1_HCO3"
    assert sinha2020.MODEL_VARIANT in sinha2020.VARIANTS_IN_CIRCULATION
    if sinha2020.CITATION_IS_PROVISIONAL:
        assert "PROVISIONAL" in sinha2020.CITATION
    else:
        # Resolved: it is an *ancillary* model from the appendix, not the
        # paper's headline 3-variable model. The citation must say so, or the
        # README will eventually claim the wrong thing.
        assert "ancillary" in sinha2020.CITATION.lower()
        assert "Table S4" in sinha2020.CITATION


def test_stnfr1_coefficient_is_the_primary_source_value():
    """1.1903 from two Sinha-authored appendices, not the 1.1913 in circulation.

    Pinned as a test because the difference is ~0.009 log-odds and would never
    be visible in any output — exactly why it needs a guard rather than a review.
    """
    assert sinha2020.BETA_LN_STNFR1 == 1.1903


def test_both_published_cutpoints_are_available_and_ordered():
    assert sinha2020.YOUDEN_THRESHOLD == 0.274
    assert sinha2020.THRESHOLD == 0.5
    assert sinha2020.YOUDEN_THRESHOLD < sinha2020.THRESHOLD


def test_published_class_medians_land_on_the_right_side_of_both_cutpoints():
    """The strongest available units check: the paper's own reference distributions."""
    lo = sinha2020.REFERENCE_MEDIANS["hypoinflammatory"]
    hi = sinha2020.REFERENCE_MEDIANS["hyperinflammatory"]
    p_lo = sinha2020.classify([lo["il6_pg_ml"]], [lo["hco3_mmol_l"]],
                              [lo["stnfr1_pg_ml"]]).probability[0]
    p_hi = sinha2020.classify([hi["il6_pg_ml"]], [hi["hco3_mmol_l"]],
                              [hi["stnfr1_pg_ml"]]).probability[0]
    assert p_lo < sinha2020.YOUDEN_THRESHOLD
    assert p_hi > sinha2020.THRESHOLD


def test_il8_variant_agrees_with_il6_variant_on_clear_cases():
    """Two ancillary models from the same table should not disagree on archetypes."""
    hyper_il6 = sinha2020.classify([933.0], [20.0], [7452.0])
    hyper_il8 = sinha2020.classify_il8([200.0], [20.0], [7452.0])
    hypo_il8 = sinha2020.classify_il8([10.0], [26.0], [2000.0])
    assert hyper_il6.label[0] == hyper_il8.label[0] == "hyperinflammatory"
    assert hypo_il8.label[0] == "hypoinflammatory"


def test_il8_variant_also_refuses_npx():
    df = _fake_npx()
    with pytest.raises(ValueError, match="unit error"):
        sinha2020.classify_il8(df.IL6, df.HCO3, df.TNFRSF1A)
