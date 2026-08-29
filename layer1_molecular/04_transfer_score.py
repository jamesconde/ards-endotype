#!/usr/bin/env python3
"""Score a FROZEN signature on held-out GEO cohorts. §3.5, C4.

"The step that prevents self-deception."

**The cohorts are never merged.** Three platforms, three batches: merging and
clustering would manufacture clusters that separate perfectly by cohort, which
§6.5 names as the single most common failure in this literature. Instead the
signature is derived on GSE65682 only, frozen, and applied to the others with
ZERO refitting.

Cross-platform transfer needs gene symbols, since probe IDs are
platform-specific. Each platform's probe->symbol map is extracted from its GEO
annotation, and probes are collapsed per gene by maximum IQR (§3.2).

Tests, per §3.5:
  GSE66890  does the score separate sepsis-with-ARDS from sepsis alone?
            (`disease state` and `berlin_ards` are both recorded)
  GSE32707  does an individual's score change from day 0 to day 7, and does the
            direction track phenotype? (`SIRS/Sepsis/se-ARDS Day 0|7`)

Single-sample scoring: each sample's score depends only on its own expression,
so no cohort statistic enters — the same discipline as the Layer 2 portable
classifier, for the same reason.

Usage
-----
    python layer1_molecular/04_transfer_score.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.logging_setup import setup, set_seeds, step

log = setup("L1_transfer")

PLATFORM = {"GSE65682": "GPL13667", "GSE32707": "GPL10558", "GSE66890": "GPL6244"}


def probe_map(gpl: str) -> pd.Series:
    f = paths.GEO / f"{gpl}_map.tsv"
    if not f.exists():
        raise CheckpointError(f"{f} missing — extract the platform annotation first")
    m = pd.read_csv(f, sep="\t", dtype=str).dropna()
    m = m[(m.symbol != "") & (m.symbol != "---")]
    # Multi-symbol entries ("A///B") are ambiguous; keep the first.
    m["symbol"] = m.symbol.str.split("///").str[0].str.strip()
    return m.set_index("probe").symbol


def to_genes(expr: pd.DataFrame, gpl: str) -> pd.DataFrame:
    """Collapse probes to genes, keeping the probe with maximum IQR (§3.2)."""
    mp = probe_map(gpl)
    idx = expr.index.astype(str)
    sym = mp.reindex(idx)
    keep = sym.notna().to_numpy()
    e = expr.loc[keep]
    s = sym[keep]
    iqr = (e.quantile(.75, axis=1) - e.quantile(.25, axis=1))
    order = iqr.groupby(s.to_numpy()).idxmax()
    out = e.loc[order.to_numpy()]
    out.index = order.index
    log.info("    %s: %d probes -> %d genes", gpl, expr.shape[0], out.shape[0])
    return out


def sample_labels(gse: str) -> pd.Series:
    md = pd.read_pickle(paths.SCRATCH / "geo" / f"{gse}_meta.pkl")
    acc = md["Sample_geo_accession"][0]
    if gse == "GSE32707":
        return pd.Series(md["Sample_source_name_ch1"][0], index=acc, name="group")
    for key, blocks in md.items():
        if "characteristics" not in key.lower():
            continue
        for vals in blocks:
            if any(str(v).lower().startswith("disease state") for v in vals):
                return pd.Series(
                    [str(v).split(":", 1)[1].strip() for v in vals],
                    index=acc, name="group")
    raise CheckpointError(f"{gse}: no phenotype field found")


def main() -> int:
    set_seeds(42)
    sig = json.loads((paths.RESULTS / "mars_signature.json").read_text())
    log.info("frozen signature: %d probes, classes %s",
             len(sig["probes"]), sig["classes"])

    with step(log, "build a gene-level signature from the frozen probe model"):
        disc = pd.read_parquet(paths.SCRATCH / "geo" / "GSE65682_expr.parquet")
        mp = probe_map(PLATFORM["GSE65682"])
        # Mars1 is the class most consistently reported as the severe/
        # hyperinflammatory-like endotype; its coefficient vector defines the axis.
        target = "Mars1" if "Mars1" in sig["classes"] else sig["classes"][0]
        coef = np.asarray(sig["coefficients"][target], dtype=float)
        probes = pd.Index(sig["probes"])
        w = pd.Series(coef, index=probes)
        w = w[w != 0]
        genes = mp.reindex(w.index.astype(str))
        gw = (pd.DataFrame({"symbol": genes.to_numpy(), "w": w.to_numpy()})
              .dropna().groupby("symbol").w.mean())
        gw = gw[gw.index.notna() & (gw.index != "---")]
        log.info("  %s axis: %d non-zero probes -> %d genes", target, len(w), len(gw))
        report_only("l1_signature_genes", len(gw))
        if len(gw) < 15:
            raise CheckpointError(f"only {len(gw)} signature genes survived mapping")

    results = {}
    for gse in ("GSE32707", "GSE66890"):
        with step(log, f"transfer to {gse}"):
            expr = pd.read_parquet(paths.SCRATCH / "geo" / f"{gse}_expr.parquet")
            g = to_genes(expr, PLATFORM[gse])
            shared = [s for s in gw.index if s in g.index]
            log.info("    signature genes present: %d of %d", len(shared), len(gw))
            if len(shared) < 10:
                log.error("    too few shared genes; skipping %s", gse)
                continue

            # Single-sample score: z within GENE across that sample's own profile
            # is not possible for one sample, so use the fixed weights on the
            # sample's own rank-normalised expression. No cohort statistic.
            sub = g.loc[shared]
            ranks = sub.rank(axis=0, pct=True)      # within each SAMPLE, across genes
            score = ranks.mul(gw[shared], axis=0).sum() / np.abs(gw[shared]).sum()

            lab = sample_labels(gse).reindex(score.index)
            df = pd.DataFrame({"score": score, "group": lab}).dropna()
            log.info("    scored %d samples; groups %s", len(df),
                     df.group.value_counts().to_dict())

            entry = {"n": int(len(df)), "genes_used": int(len(shared)),
                     "group_means": df.groupby("group").score.mean().round(4).to_dict()}

            if gse == "GSE66890":
                # EXACT match only. The exclusion label literally contains the
                # phrase "sepsis with ARDS" ("not meeting the Berlin definition
                # of sepsis with ARDS and that they were excluded from
                # analysis"), so a substring match silently folded 5 explicitly
                # EXCLUDED patients into the ARDS group.
                ards = df[df.group.str.fullmatch("sepsis with ARDS", case=False, na=False)]
                sep = df[df.group.str.fullmatch("sepsis alone", case=False, na=False)]
                excluded = df[df.group.str.contains("excluded", case=False, na=False)]
                if len(excluded):
                    log.info("    excluded by the authors and left out: %d", len(excluded))
                if len(ards) > 5 and len(sep) > 5:
                    u = stats.mannwhitneyu(ards.score, sep.score)
                    cliff = 2 * u.statistic / (len(ards) * len(sep)) - 1
                    log.warning("    sepsis+ARDS (n=%d) vs sepsis alone (n=%d): "
                                "Cliff's delta %+.3f, p=%.4g",
                                len(ards), len(sep), cliff, u.pvalue)
                    entry.update({"cliffs_delta": round(float(cliff), 4),
                                  "p": float(f"{u.pvalue:.4g}"),
                                  "n_ards": int(len(ards)), "n_sepsis": int(len(sep))})
                    if u.pvalue >= 0.05:
                        log.warning("    NOT significant. The frozen signature does not "
                                    "separate ARDS from sepsis alone in this cohort. "
                                    "Reported as a null (I5) — the signature was "
                                    "derived to predict MARS endotype, not ARDS status, "
                                    "so this is a demanding test rather than a failure.")

            if gse == "GSE32707":
                df["day"] = df.group.str.extract(r"Day (\d)")[0]
                df["pheno"] = df.group.str.replace(r"\s*Day \d", "", regex=True)
                by = df.groupby(["pheno", "day"]).score.agg(["mean", "count"]).round(4)
                log.info("    score by phenotype and day:\n%s", by.to_string())
                entry["by_pheno_day"] = by.reset_index().to_dict("records")
                for ph, gg in df.groupby("pheno"):
                    d0 = gg[gg.day == "0"].score
                    d7 = gg[gg.day == "7"].score
                    if len(d0) > 4 and len(d7) > 4:
                        u = stats.mannwhitneyu(d7, d0)
                        log.info("      %-14s day0 %.3f -> day7 %.3f  p=%.3g",
                                 ph, d0.mean(), d7.mean(), u.pvalue)

            results[gse] = entry
            df.groupby("group").score.describe().round(4).to_csv(
                paths.RESULTS / f"L1_transfer_{gse}.csv")

    (paths.RESULTS / "L1_transfer_summary.json").write_text(json.dumps({
        "signature_class": target,
        "signature_genes": int(len(gw)),
        "cohorts": results,
        "design": ("Signature frozen on GSE65682 and applied with zero refitting. "
                   "Cohorts never merged (§3.5, §6.5). Single-sample scoring, so no "
                   "cohort statistic enters the score."),
    }, indent=2, default=str))
    log.info("wrote results/L1_transfer_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
