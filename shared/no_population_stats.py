"""Make the portability claim executable. SPEC.md §4.5a2, §6.12.

A "portable" classifier is one whose prediction for a patient depends only on
that patient's own values. The failure mode the spec names explicitly is
leaving a fitted `StandardScaler` in the inference path and calling the model
portable anyway.

Two independent checks, because either alone is weak:

  1. `assert_no_population_stats(model)` — static introspection. Walks a fitted
     sklearn estimator/Pipeline and rejects any component that carries a
     cohort-derived statistic used at transform time.

  2. `assert_batch_invariant(predict_fn, X)` — functional. Predicts each row
     alone and again inside the full batch; the two must agree exactly. Any
     dependence on cohort statistics computed at inference shows up here even
     if it hides in a custom function that no introspection would catch.

Check 2 is the one to trust. Check 1 explains *why* when check 2 fails.

Distinction that matters: a parameter *learned at fit time and frozen*
(a logistic coefficient, a tree's split threshold on a raw value) is fine —
it is a fixed number, like the Sinha 2020 intercept. A statistic *of the
cohort being scored* (a z-score, a median imputation, a cohort quantile) is
not, because it makes one patient's prediction depend on who else was admitted.
"""
from __future__ import annotations

import logging
from typing import Callable, Sequence

import numpy as np

log = logging.getLogger(__name__)


class PortabilityError(AssertionError):
    """The inference path depends on the cohort. It is not deployable."""


# Attributes that exist only because a transformer memorised a cohort moment
# and will re-apply it at transform time.
BANNED_ATTRS = {
    "mean_": "cohort mean (StandardScaler / PCA centring)",
    "var_": "cohort variance",
    "scale_": "cohort scale factor",
    "center_": "cohort centre (RobustScaler)",
    "quantiles_": "cohort quantiles (QuantileTransformer)",
    "data_min_": "cohort minimum (MinMaxScaler)",
    "data_max_": "cohort maximum (MinMaxScaler)",
    "data_range_": "cohort range (MinMaxScaler)",
    "statistics_": "cohort imputation statistic (SimpleImputer)",
    "lambdas_": "cohort Box-Cox/Yeo-Johnson lambdas (PowerTransformer)",
    "bin_edges_": "cohort-derived bin edges (KBinsDiscretizer)",
    "components_": "cohort covariance structure (PCA/factor rotation)",
    "n_features_outlier_": "cohort outlier structure",
}

# Estimators whose *whole purpose* is a population transform.
BANNED_CLASSES = {
    "StandardScaler", "MinMaxScaler", "RobustScaler", "MaxAbsScaler",
    "QuantileTransformer", "PowerTransformer", "Normalizer",
    "KBinsDiscretizer", "SimpleImputer", "IterativeImputer", "KNNImputer",
    "PCA", "TruncatedSVD", "FactorAnalysis", "KernelPCA", "IncrementalPCA",
}


def _walk(obj, path: str, hits: list[str], seen: set[int]) -> None:
    if id(obj) in seen:
        return
    seen.add(id(obj))

    cls = type(obj).__name__
    if cls in BANNED_CLASSES:
        hits.append(f"{path}: {cls} — a population transform in the inference path")

    for attr, why in BANNED_ATTRS.items():
        if hasattr(obj, attr) and getattr(obj, attr) is not None:
            hits.append(f"{path}.{attr}: {why}")

    # Descend into pipelines, column transformers, ensembles, containers.
    for name in ("steps", "transformers_", "transformers", "named_steps",
                 "estimators_", "estimators", "transformer_list", "_final_estimator"):
        child = getattr(obj, name, None)
        if child is None:
            continue
        items = child.items() if isinstance(child, dict) else child
        if not isinstance(items, (list, tuple)) and not hasattr(items, "__iter__"):
            _walk(child, f"{path}.{name}", hits, seen)
            continue
        for i, item in enumerate(items):
            sub = item
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                sub = item[1]
            if hasattr(sub, "__dict__") or hasattr(sub, "get_params"):
                _walk(sub, f"{path}.{name}[{i}]", hits, seen)


def assert_no_population_stats(model, label: str = "model") -> None:
    """Reject any fitted cohort statistic reachable from `model`."""
    hits: list[str] = []
    _walk(model, label, hits, set())
    if hits:
        raise PortabilityError(
            f"{label} is NOT portable — the inference path carries cohort statistics:\n  "
            + "\n  ".join(sorted(set(hits)))
            + "\n\nA deployed screening model sees one patient at a time and has no "
              "cohort to standardise against. Refit on raw values (§4.5a2)."
        )
    log.info("✓ %s: no population statistics in the inference path", label)


def assert_batch_invariant(
    predict_fn: Callable[[np.ndarray], np.ndarray],
    X: np.ndarray,
    n_probe: int = 64,
    rtol: float = 0.0,
    atol: float = 1e-10,
    label: str = "model",
    seed: int = 0,
) -> None:
    """Each row predicted alone must equal that row predicted in the batch.

    This is the check that actually binds. A model can pass introspection and
    still normalise inside a custom predict function; it cannot pass this.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"expected 2-D X, got shape {X.shape}")

    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), size=min(n_probe, len(X)), replace=False)

    batch = np.asarray(predict_fn(X), dtype=float)
    if batch.ndim == 1:
        batch = batch[:, None]

    bad: list[str] = []
    for i in idx:
        alone = np.asarray(predict_fn(X[[i]]), dtype=float).reshape(1, -1)
        if not np.allclose(alone, batch[[i]], rtol=rtol, atol=atol, equal_nan=True):
            bad.append(f"row {i}: alone={alone.ravel()[:4]} batch={batch[i][:4]}")

    if bad:
        raise PortabilityError(
            f"{label} is NOT portable — predictions depend on the batch "
            f"({len(bad)}/{len(idx)} probe rows disagree):\n  "
            + "\n  ".join(bad[:8])
            + "\n\nSomething in the inference path is computed across rows."
        )
    log.info("✓ %s: batch-invariant over %d probe rows", label, len(idx))


def assert_portable(model, predict_fn, X, label: str = "model") -> None:
    """Both checks. This is what the Layer 2 scripts call."""
    assert_no_population_stats(model, label)
    assert_batch_invariant(predict_fn, X, label=label)


def describe_inference_path(model, label: str = "model") -> list[str]:
    """Enumerate the fitted parameters a deployment would have to ship.

    Not a check — reporting. Goes in the README so a reader can see exactly
    what the portable model consists of.
    """
    out: list[str] = []

    def rec(obj, path: str, seen: set[int]) -> None:
        if id(obj) in seen:
            return
        seen.add(id(obj))
        cls = type(obj).__name__
        for attr in ("coef_", "intercept_", "feature_names_in_", "tree_",
                     "classes_", "n_features_in_"):
            if hasattr(obj, attr):
                val = getattr(obj, attr)
                shape = getattr(val, "shape", None)
                out.append(f"{path} [{cls}] {attr}"
                           + (f" shape={shape}" if shape is not None else ""))
        for name in ("steps", "estimators_", "transformers_"):
            child = getattr(obj, name, None)
            if child is None:
                continue
            for i, item in enumerate(child):
                sub = item[1] if isinstance(item, (list, tuple)) and len(item) >= 2 else item
                rec(sub, f"{path}.{name}[{i}]", seen)

    rec(model, label, set())
    return out
