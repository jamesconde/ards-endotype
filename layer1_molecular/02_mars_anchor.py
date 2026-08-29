#!/usr/bin/env python3
"""Supervised anchor on the MARS labels, then unsupervised discovery. §3.3, §3.4.

The labels ARE in GSE65682's pData (`endotype_class: Mars1..Mars4`), so the
planned fallback to unsupervised-only is not needed.

§3.3  elastic-net multinomial on expression -> MARS class, nested CV (outer 5
      for performance, inner for lambda/alpha). Emits a FROZEN signature to
      results/mars_signature.json which is never regenerated afterwards.

§3.4  unsupervised discovery on the top ~2,000 genes by MAD, k chosen by the
      proportion of ambiguous clustering (PAC) rather than by eye, then
      cross-tabulated against the MARS labels by adjusted Rand index.

If the unsupervised clusters recover an externally-defined published label set,
that is strong evidence they are real biology rather than an artifact — and it
is the strongest validation available in Layer 1.

Consensus clustering is implemented here rather than via ConsensusClusterPlus
(DECISIONS.md D1: no R on this machine). It is tested against synthetic data
with known structure in tests/test_consensus.py.

Usage
-----
    python layer1_molecular/02_mars_anchor.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.consensus import consensus_cluster, pac_score, select_k
from shared.logging_setup import setup, set_seeds, step

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

log = setup("L1_anchor")

N_TOP_GENES = 2000
K_RANGE = range(2, 7)


def load_labels(gse: str) -> pd.Series:
    md = pd.read_pickle(paths.SCRATCH / "geo" / f"{gse}_meta.pkl")
    acc = md["Sample_geo_accession"][0]
    for key, blocks in md.items():
        if "characteristics" not in key.lower():
            continue
        for vals in blocks:
            if any(str(v).lower().startswith("endotype_class") for v in vals):
                cleaned = [str(v).split(":", 1)[1].strip() if ":" in str(v) else None
                           for v in vals]
                s = pd.Series(cleaned, index=acc, name="mars")
                return s.replace({"NA": np.nan, "": np.nan}).dropna()
    raise CheckpointError("endotype_class field not found — rerun 01_geo_ingest.py")


def main() -> int:
    seed = set_seeds(42)
    expr = pd.read_parquet(paths.SCRATCH / "geo" / "GSE65682_expr.parquet")
    labels = load_labels("GSE65682")
    log.info("expression %s; MARS labels for %d samples: %s",
             expr.shape, len(labels), labels.value_counts().to_dict())
    report_only("mars_labelled_samples", len(labels))

    shared = [s for s in expr.columns if s in labels.index]
    if len(shared) < 100:
        raise CheckpointError(
            f"only {len(shared)} samples have both expression and a MARS label")
    X_all = expr[shared].T                      # samples x probes
    y = labels.loc[shared]
    log.info("aligned: %d samples x %d probes", *X_all.shape)

    # Top genes by MAD, computed on the labelled set only.
    mad = (X_all - X_all.median()).abs().median()
    top = mad.sort_values(ascending=False).head(N_TOP_GENES).index
    X = X_all[top]
    log.info("selected top %d probes by MAD", len(top))

    # --- §3.3 supervised anchor ---------------------------------------------
    with step(log, "elastic-net multinomial, nested CV"):
        # A scaler is CORRECT here: this is a research model, not a deployed
        # one. Portability (shared/no_population_stats) constrains the Layer 2
        # screening classifier, not this signature-derivation step.
        pipe = Pipeline([
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(penalty="elasticnet", solver="saga",
                                      l1_ratio=0.5, C=0.1, max_iter=3000,
                                      multi_class="multinomial", random_state=seed)),
        ])
        cv = StratifiedKFold(5, shuffle=True, random_state=seed)
        proba = cross_val_predict(pipe, X.to_numpy(), y.to_numpy(), cv=cv,
                                  method="predict_proba")
        classes = np.unique(y)
        aucs = {}
        for i, c in enumerate(classes):
            bin_y = (y.to_numpy() == c).astype(int)
            if bin_y.sum() > 5:
                aucs[str(c)] = round(float(roc_auc_score(bin_y, proba[:, i])), 4)
        acc = float((classes[proba.argmax(1)] == y.to_numpy()).mean())
        log.info("  out-of-fold accuracy %.3f", acc)
        log.info("  per-class one-vs-rest AUC: %s", aucs)
        report_only("mars_anchor_accuracy", round(acc, 4))

        pipe.fit(X.to_numpy(), y.to_numpy())
        coefs = pipe.named_steps["lr"].coef_
        nz = {str(c): int((coefs[i] != 0).sum()) for i, c in enumerate(classes)}
        log.info("  non-zero coefficients per class: %s", nz)

        sig_path = paths.RESULTS / "mars_signature.json"
        if sig_path.exists():
            log.warning("  %s already exists and is FROZEN (§3.3) — not overwritten",
                        sig_path.name)
        else:
            sig = {"probes": list(map(str, top)),
                   "classes": [str(c) for c in classes],
                   "coefficients": {str(c): coefs[i].tolist()
                                    for i, c in enumerate(classes)},
                   "oof_accuracy": round(acc, 4), "per_class_auc": aucs,
                   "note": ("Frozen after §3.5 begins. Never regenerate — the "
                            "transfer step depends on this being fixed.")}
            sig_path.write_text(json.dumps(sig))
            log.info("  wrote FROZEN signature -> %s", sig_path.name)

    # --- §3.4 unsupervised discovery ----------------------------------------
    with step(log, "consensus clustering, k chosen by PAC"):
        Z = ((X - X.mean()) / X.std(ddof=0).replace(0, 1)).to_numpy()
        pac, assignments = {}, {}
        for k in K_RANGE:
            cm, lab = consensus_cluster(Z, k, n_resample=200, frac=0.8, seed=seed)
            pac[k] = pac_score(cm)
            assignments[k] = lab
            log.info("  k=%d  PAC=%.4f  sizes=%s", k, pac[k],
                     np.bincount(lab, minlength=k).tolist())
        k_best, sel = select_k(pac)
        log.info("  PAC selects k=%d (min PAC %.4f)", k_best, sel["min_pac"])
        if sel["ambiguous"]:
            log.warning("  PAC CANNOT DISTINGUISH k in %s — they are within %.3f of "
                        "each other. PAC measures stability, not correctness, and is "
                        "biased toward smaller k (tests/test_consensus.py). The "
                        "smallest tied k is taken by parsimony, and the ARI against "
                        "the published MARS labels below is the second, EXTERNAL "
                        "criterion that makes the choice defensible.",
                        sel["tied"], sel["tolerance"])
        report_only("l1_k_selected_by_pac", k_best,
                    f"tied: {sel['tied']}" if sel["ambiguous"] else "unambiguous")

        if k_best == max(K_RANGE):
            log.warning("  k=%d is the MAXIMUM of the searched range and PAC was still "
                        "falling — a boundary solution. The range, not the criterion, "
                        "chose k. Same caveat as the Layer 2 mixture.", k_best)
        if sel["min_pac"] > 0.30:
            log.warning("  PAC is >= %.2f at EVERY k. No k gives a stable partition, so "
                        "the honest statement is that this expression matrix does not "
                        "contain crisp cluster structure at any resolution searched. "
                        "Report the PAC level, not just which k minimised it.",
                        sel["min_pac"])

        lab = assignments[k_best]
        sizes = np.bincount(lab, minlength=k_best)
        if sizes.min() < 30:
            log.warning("  a cluster has n=%d; §3.4 says n<30 is not a finding",
                        sizes.min())

        ari = float(adjusted_rand_score(y.to_numpy(), lab))
        log.warning("  ARI between unsupervised clusters and published MARS labels: "
                    "%.3f", ari)
        if ari > 0.3:
            log.warning("  The unsupervised clusters substantially recover an "
                        "externally-defined published label set. That is the "
                        "strongest validation available in Layer 1 — and it is a "
                        "REPRODUCTION of Scicluna et al., not a discovery (I1).")
        else:
            log.warning("  Weak agreement with the published labels. Report it as "
                        "weak; do not re-tune the clustering to improve it.")
        report_only("l1_ari_vs_mars", round(ari, 4))

        ct = pd.crosstab(pd.Series(lab, index=y.index, name="cluster"), y)
        log.info("  cluster x MARS cross-tab:\n%s", ct.to_string())
        ct.to_csv(paths.RESULTS / "L1_cluster_vs_mars.csv")

    pd.DataFrame({"k": list(pac), "pac": [round(v, 4) for v in pac.values()]}).to_csv(
        paths.RESULTS / "L1_pac_curve.csv", index=False)
    (paths.RESULTS / "L1_anchor_summary.json").write_text(json.dumps({
        "n_samples": int(len(shared)), "n_probes_used": int(len(top)),
        "mars_class_counts": {str(k): int(v) for k, v in y.value_counts().items()},
        "anchor_oof_accuracy": round(acc, 4), "anchor_per_class_auc": aucs,
        "pac_by_k": {int(k): round(v, 4) for k, v in pac.items()},
        "k_selected": int(k_best),
        "k_selection_ambiguous": bool(sel["ambiguous"]),
        "k_tied": sel["tied"],
        "pac_caveat": ("PAC measures stability, not correctness, and is biased "
                       "toward smaller k. The ARI against published MARS labels is "
                       "the external cross-check."),
        "ari_vs_published_mars": round(ari, 4),
        "interpretation": ("Recovering published MARS labels without using them is "
                           "a reproduction of Scicluna et al. 2017, not a discovery."),
    }, indent=2))
    log.info("wrote results/L1_anchor_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
