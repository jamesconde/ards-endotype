"""Consensus clustering with PAC-based k selection. SPEC.md §3.4.

Replaces ConsensusClusterPlus, which is R-only (DECISIONS.md D1). The algorithm
is the standard one (Monti et al. 2003):

  repeat n_resample times:
      subsample a fraction of the samples
      cluster them into k groups
      record which pairs landed together
  consensus[i,j] = (times i and j co-clustered) / (times both were sampled)

k is then chosen by the **proportion of ambiguous clustering** (PAC): the
fraction of consensus values falling strictly between two bounds. A perfect
clustering puts every pair at 0 or 1, so a low PAC means a stable solution.
§3.4 requires k be selected by a stated criterion rather than by eyeballing a
heatmap, and PAC is that criterion.

D1 records the cost of this decision: the implementation is ours, not the
canonical package, so it carries its own test against synthetic data with known
structure (tests/test_consensus.py).
"""
from __future__ import annotations

import numpy as np
from sklearn.cluster import AgglomerativeClustering


def _base_cluster(sub: np.ndarray, k: int, algo: str, metric: str, seed: int):
    """One base clustering. Defaults follow transcriptomics practice.

    Euclidean distance with average linkage CHAINS badly on high-dimensional
    expression data: it returns one giant cluster plus singletons, which is a
    degenerate partition rather than an unfavourable one. ConsensusClusterPlus
    defaults to a correlation distance for exactly this reason, so `pearson` is
    the default here too.
    """
    if algo == "km":
        from sklearn.cluster import KMeans
        return KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(sub)
    if metric == "pearson":
        c = np.corrcoef(sub)
        c = np.nan_to_num(c, nan=0.0)
        d = 1.0 - c
        np.fill_diagonal(d, 0.0)
        return AgglomerativeClustering(
            n_clusters=k, metric="precomputed", linkage="average").fit_predict(d)
    return AgglomerativeClustering(n_clusters=k, linkage="average").fit_predict(sub)


def consensus_matrix(X: np.ndarray, k: int, n_resample: int = 200,
                     frac: float = 0.8, seed: int = 42,
                     linkage: str = "average", algo: str = "hc",
                     metric: str = "pearson") -> np.ndarray:
    """Pairwise co-clustering frequency over subsamples."""
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    co = np.zeros((n, n), dtype=np.float64)     # times clustered together
    both = np.zeros((n, n), dtype=np.float64)   # times both were sampled

    size = max(int(round(frac * n)), k + 1)
    for _ in range(n_resample):
        idx = rng.choice(n, size=size, replace=False)
        sub = X[idx]
        lab = _base_cluster(sub, k, algo, metric, seed)
        both[np.ix_(idx, idx)] += 1
        for c in range(k):
            members = idx[lab == c]
            if len(members) > 1:
                co[np.ix_(members, members)] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        cm = np.where(both > 0, co / both, 0.0)
    np.fill_diagonal(cm, 1.0)
    return cm


def pac_score(cm: np.ndarray, lower: float = 0.1, upper: float = 0.9) -> float:
    """Proportion of ambiguous clustering — lower is more stable.

    Only the strict upper triangle is used, so the diagonal of ones and the
    duplicated lower triangle cannot dilute the estimate.
    """
    iu = np.triu_indices_from(cm, k=1)
    v = cm[iu]
    if v.size == 0:
        return float("nan")
    return float(((v > lower) & (v < upper)).mean())


def consensus_cluster(X: np.ndarray, k: int, n_resample: int = 200,
                      frac: float = 0.8, seed: int = 42,
                      algo: str = "hc", metric: str = "pearson"
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Return (consensus matrix, final labels from the consensus distance)."""
    cm = consensus_matrix(X, k, n_resample, frac, seed, algo=algo, metric=metric)
    dist = 1.0 - cm
    np.fill_diagonal(dist, 0.0)
    labels = AgglomerativeClustering(
        n_clusters=k, metric="precomputed", linkage="average").fit_predict(dist)
    return cm, labels


def select_k(pac_by_k: dict[int, float], tol: float = 0.01
             ) -> tuple[int, dict]:
    """Choose k by minimum PAC, reporting ties instead of hiding them.

    PAC measures STABILITY, not correctness. On well-separated data several k
    can be equally stable — merging two distant clusters is still an unambiguous
    assignment — so `min(pac)` can silently return an arbitrary member of a tie.
    tests/test_consensus.py pins that behaviour with a case where k=2 and k=3
    both score exactly 0.0.

    Convention adopted: among k within `tol` of the minimum, take the SMALLEST
    (parsimony). The tie is returned so it can be reported rather than assumed
    away, which §3.4 requires — k must be selected by a stated criterion, and
    "the criterion could not distinguish these" is part of stating it.
    """
    best = min(pac_by_k.values())
    tied = sorted(k for k, v in pac_by_k.items() if v <= best + tol)
    return tied[0], {"min_pac": best, "tied": tied, "tolerance": tol,
                     "ambiguous": len(tied) > 1}
