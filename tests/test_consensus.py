"""Tests for the consensus-clustering implementation.

DECISIONS.md D1 states the cost of replacing ConsensusClusterPlus: the code is
ours, so it needs its own evidence against data with known structure. These are
that evidence. Without them, a PAC curve could be wrong in a way that looks
entirely plausible.
"""
from __future__ import annotations

import numpy as np
import pytest

from shared.consensus import consensus_cluster, consensus_matrix, pac_score


def blobs(n_per: int, k: int, dim: int = 20, sep: float = 8.0, seed: int = 0):
    """k well-separated Gaussian blobs with known ground-truth labels."""
    rng = np.random.default_rng(seed)
    centres = rng.normal(0, sep, size=(k, dim))
    X = np.vstack([centres[c] + rng.normal(0, 1.0, size=(n_per, dim))
                   for c in range(k)])
    y = np.repeat(np.arange(k), n_per)
    return X, y


def test_recovers_known_clusters():
    X, y = blobs(40, 3)
    _, lab = consensus_cluster(X, 3, n_resample=60, seed=1, metric='euclidean')
    from sklearn.metrics import adjusted_rand_score
    assert adjusted_rand_score(y, lab) > 0.95


def test_pac_ties_on_perfectly_separated_data():
    """A real limitation of PAC, pinned so it cannot surprise us later.

    With well-separated blobs, PAC is 0.0 for k=2 AND k=3: merging two clusters
    that are far apart still yields a perfectly STABLE assignment, and stability
    is all PAC measures. So PAC cannot distinguish k here, and any code that
    takes `min(pac)` silently returns the smallest tied k.

    This is why `select_k` reports ties rather than hiding them.
    """
    X, _ = blobs(40, 3, sep=8.0)
    pac = {k: pac_score(consensus_cluster(X, k, n_resample=60, seed=1, metric='euclidean')[0])
           for k in (2, 3, 4, 5)}
    assert pac[2] == pytest.approx(pac[3], abs=1e-9), (
        f"expected a tie between k=2 and k=3 on separated data: {pac}")


def test_pac_is_biased_toward_smaller_k():
    """Documents a real limitation rather than engineering around it.

    Across a wide range of cluster separations, PAC does NOT recover the true
    k=3: it scores k=2 and k=3 identically at 0.0 and rises for larger k. The
    reason is structural — average-linkage agglomerative clustering stably
    merges two of the three blobs, and a stable merge is unambiguous, which is
    the only thing PAC measures.

    CONSEQUENCE FOR §3.4: PAC alone cannot be the sole arbiter of k. It rules
    OUT unstable solutions, which is genuinely useful, but among stable ones it
    prefers the smallest. That is why the analysis reports the whole PAC curve,
    reports ties, and cross-checks the chosen k against the published MARS
    labels by adjusted Rand index — a second, external criterion.

    Tuning the blob separation until this test passed would have hidden a
    property that will affect the real result.
    """
    for sep in (1.5, 2.5, 3.5):
        X, _ = blobs(50, 3, sep=sep, seed=5)
        pac = {k: pac_score(consensus_cluster(X, k, n_resample=60, seed=1, metric='euclidean')[0])
               for k in (2, 3, 4, 5)}
        assert pac[2] <= pac[4], f"PAC should not favour k=4 over k=2 here: {pac}"
        assert min(pac, key=pac.get) == 2, (
            f"expected PAC's small-k bias at sep={sep}, got {pac}")


def test_select_k_flags_ties():
    from shared.consensus import select_k
    k, info = select_k({2: 0.0, 3: 0.0, 4: 0.05})
    assert k == 2 and info["tied"] == [2, 3]
    k, info = select_k({2: 0.30, 3: 0.02, 4: 0.19})
    assert k == 3 and info["tied"] == [3]


def test_pac_higher_on_structureless_data():
    """Unstructured data must look ambiguous, or PAC means nothing."""
    rng = np.random.default_rng(3)
    noise = rng.normal(size=(120, 20))
    X, _ = blobs(40, 3)
    cm_noise, _ = consensus_cluster(noise, 3, n_resample=60, seed=1, metric='euclidean')
    cm_real, _ = consensus_cluster(X, 3, n_resample=60, seed=1, metric='euclidean')
    assert pac_score(cm_noise) > pac_score(cm_real)


def test_consensus_matrix_is_a_valid_similarity():
    X, _ = blobs(25, 2)
    cm = consensus_matrix(X, 2, n_resample=40, seed=2, metric='euclidean')
    assert cm.shape == (50, 50)
    assert np.allclose(np.diag(cm), 1.0)
    assert np.allclose(cm, cm.T), "consensus must be symmetric"
    assert cm.min() >= 0.0 and cm.max() <= 1.0


def test_pac_ignores_diagonal():
    """A diagonal of ones must not be counted as unambiguous agreement."""
    cm = np.full((10, 10), 0.5)
    np.fill_diagonal(cm, 1.0)
    assert pac_score(cm) == pytest.approx(1.0)
