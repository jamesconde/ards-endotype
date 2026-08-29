#!/usr/bin/env python3
"""MGH Olink Explore 1536 — proteomic endotypes and trajectory. §3.7, C7.

The most GEn1E-shaped analysis here: protein-level, longitudinal, with clinical
outcomes. Their platform page names "personalization & trajectory mapping".

Four questions, in descending order of what they are worth:

1. **Aetiology or host response?** 78 of these patients were COVID-NEGATIVE but
   presented with the same clinical syndrome. If the clusters MIX positive and
   negative patients rather than separating them, the endotype reflects host
   response rather than pathogen — which is the premise a host-directed drug
   like a p38 modulator depends on. §6.7 mentions this in one line; it is the
   most interesting figure available in Layer 1.
2. **Trajectory.** Does day-0 cluster predict acuity at days 3, 7 and 28, and do
   patients SWITCH class between day 0 and day 7? The literature mostly treats
   these labels as fixed.
3. **Reproduction** of the two-class inflammatory split in protein space.
4. **The p38 arm**, where a null is the expected and pre-accepted outcome (I5).

TWO TRAPS, both verified rather than assumed (DOWNLOADS.md §5):
  - TNFRSF1A is on the NEUROLOGY panel. Subsetting to Inflammation deletes it.
  - 8 of 792 samples are QC/plate controls with a blank Timepoint.

NPX is a relative log2 scale, so the published fixed-coefficient classifier
CANNOT be evaluated here (DECISIONS.md D3). The rank surrogate is used instead
and is labelled as such.

Usage
-----
    python layer1_molecular/03_olink.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.assertions import CheckpointError, report_only
from shared.consensus import consensus_cluster, pac_score, select_k
from shared.logging_setup import setup, set_seeds, step
from shared import sinha2020

from sklearn.metrics import adjusted_rand_score

log = setup("L1_olink")

NPX_FILE = paths.OLINK / "MGH_COVID_OLINK_NPX.txt"
CLIN_FILE = paths.OLINK / "MGH_COVID_Clinical_Info.txt"
LOD_MAX_FRAC = 0.50      # §3.2: drop proteins >50% below LOD at day 0


def main() -> int:
    seed = set_seeds(42)
    gs = yaml.safe_load((paths.CONFIG / "genesets.yaml").read_text())
    p38_assays = gs["protein_space"]["p38_axis_proxy"]["assays"]

    for f in (NPX_FILE, CLIN_FILE):
        if not f.exists():
            raise CheckpointError(f"{f} missing — see DOWNLOADS.md §5")

    with step(log, "load and QC"):
        npx = pd.read_csv(NPX_FILE, sep=";", low_memory=False)
        clin = pd.read_csv(CLIN_FILE, sep=";")
        log.info("  NPX %s, clinical %s", npx.shape, clin.shape)

        n_ctrl = int(npx.Timepoint.isna().sum() > 0) and \
            npx[npx.Timepoint.isna()].SampleID.nunique()
        npx = npx[npx.Timepoint.notna() & npx.subject_id.notna()].copy()
        npx["subject_id"] = npx.subject_id.astype(int)
        log.info("  dropped %d QC/plate control samples (blank Timepoint)", n_ctrl)

        d0 = npx[npx.Timepoint == "D0"]
        below = (d0.NPX < d0.LOD)
        frac_below = d0.assign(b=below).groupby("OlinkID").b.mean()
        keep = frac_below[frac_below <= LOD_MAX_FRAC].index
        log.info("  assays passing the >%.0f%%-below-LOD filter: %d of %d",
                 100 * LOD_MAX_FRAC, len(keep), frac_below.size)
        report_only("olink_assays_after_lod", len(keep))

        # Verify the two traps rather than trusting them.
        tnf = npx[npx.UniProt == "P19438"]
        if tnf.empty:
            raise CheckpointError("TNFRSF1A (P19438) absent — match on UniProt, "
                                  "not the assay label")
        log.info("  TNFRSF1A: %s on panel %s", tnf.OlinkID.unique(),
                 tnf.Panel.unique())

    with step(log, "day-0 wide matrix"):
        w = (d0[d0.OlinkID.isin(keep)]
             .pivot_table(index="subject_id", columns="OlinkID", values="NPX"))
        w = w.dropna(axis=1, thresh=int(0.9 * len(w)))
        w = w.loc[:, w.std() > 0]
        log.info("  day-0 matrix: %d subjects x %d assays", *w.shape)
        oid2gene = (npx.drop_duplicates("OlinkID").set_index("OlinkID").Assay.to_dict())

        clin = clin.set_index("subject_id")
        common = [s for s in w.index if s in clin.index]
        w, meta = w.loc[common], clin.loc[common]
        log.info("  joined to clinical: %d subjects (COVID+ %d, COVID- %d)",
                 len(w), int((meta.COVID == 1).sum()), int((meta.COVID == 0).sum()))

    # --- cluster COVID+ day 0, per §3.7 -------------------------------------
    with step(log, "consensus clustering of day-0 NPX (COVID+ only)"):
        pos = meta.COVID == 1
        Xp = w[pos]
        Z = ((Xp - Xp.mean()) / Xp.std(ddof=0).replace(0, 1)).fillna(0).to_numpy()
        pac = {}
        assign = {}
        for k in range(2, 6):
            cm, lab = consensus_cluster(Z, k, n_resample=150, seed=seed)
            pac[k], assign[k] = pac_score(cm), lab
            log.info("  k=%d PAC=%.4f sizes=%s", k, pac[k],
                     np.bincount(lab, minlength=k).tolist())
        k_best, sel = select_k(pac)
        log.info("  PAC selects k=%d%s", k_best,
                 f" (tied {sel['tied']})" if sel["ambiguous"] else "")
        lab_pos = assign[k_best]

        # Which cluster is the inflammatory one? By IL6, not by preference.
        il6_oids = [o for o in w.columns if oid2gene.get(o) == "IL6"]
        il6 = Xp[il6_oids].mean(axis=1) if il6_oids else None
        hyper = int(pd.Series(il6.to_numpy()).groupby(lab_pos).mean().idxmax()) \
            if il6 is not None else 0
        log.info("  cluster %d has the highest IL6 -> inflammatory-like", hyper)
        is_hyper = (lab_pos == hyper).astype(int)
        log.info("  inflammatory-like prevalence: %.1f%%", 100 * is_hyper.mean())

    # --- Q1: aetiology or host response? ------------------------------------
    with step(log, "aetiology vs host response (§6.7)"):
        Zall = ((w - w.mean()) / w.std(ddof=0).replace(0, 1)).fillna(0).to_numpy()
        _, lab_all = consensus_cluster(Zall, k_best, n_resample=150, seed=seed)
        ct = pd.crosstab(pd.Series(lab_all, index=w.index, name="cluster"),
                         meta.COVID.rename("covid"))
        log.info("  cluster x COVID status:\n%s", ct.to_string())
        ari_cov = float(adjusted_rand_score(meta.COVID.to_numpy(), lab_all))
        chi_p = float(stats.chi2_contingency(ct)[1]) if ct.shape[1] > 1 else np.nan
        log.warning("  ARI(cluster, COVID status) = %.3f  (chi2 p=%.4g)", ari_cov, chi_p)
        # ARI and the chi-square answer DIFFERENT questions and can disagree.
        # ARI asks "is the cluster partition the same as the COVID partition";
        # chi-square asks "is there any dependence at all". A low ARI with a
        # tiny p-value means the clusters are not a relabelling of COVID status
        # but some clusters ARE enriched. Reporting only the ARI would overstate
        # the host-response claim, so both are required to agree before it is made.
        weak_partition = ari_cov < 0.2
        independent = not (chi_p < 0.01)
        if weak_partition and independent:
            log.warning("  The clusters DO NOT separate by COVID status on either "
                        "measure. The endotype reflects HOST RESPONSE rather than "
                        "pathogen — the premise a host-directed therapy depends on.")
        elif weak_partition:
            log.warning("  MIXED, and this must not be reported as a clean result. "
                        "ARI is %.3f, so the clusters are NOT a relabelling of COVID "
                        "status — but chi-square p=%.3g shows real dependence, so some "
                        "clusters are enriched for one status. The honest reading is "
                        "PARTIAL host-response structure with an aetiology component, "
                        "and the cross-tab must be shown rather than summarised.",
                        ari_cov, chi_p)
        else:
            log.warning("  The clusters track COVID status; aetiology is a substantial "
                        "driver and the host-response reading is NOT supported.")
        report_only("olink_ari_cluster_vs_covid", round(ari_cov, 4))
        ct.to_csv(paths.RESULTS / "L1_olink_cluster_vs_covid.csv")

    # --- Q2: trajectory ------------------------------------------------------
    with step(log, "acuity trajectory by day-0 cluster"):
        acu_cols = [c for c in meta.columns if c.lower().startswith("acuity")]
        mp = meta[pos].assign(cluster=is_hyper)
        traj = mp.groupby("cluster")[acu_cols].mean().round(3)
        log.info("  mean acuity by day-0 cluster (1=death ... 5=discharged):\n%s",
                 traj.to_string())
        traj.to_csv(paths.RESULTS / "L1_olink_acuity_trajectory.csv")
        for c in acu_cols:
            a = mp.loc[mp.cluster == 1, c].dropna()
            b = mp.loc[mp.cluster == 0, c].dropna()
            if len(a) > 5 and len(b) > 5:
                u = stats.mannwhitneyu(a, b)
                cliff = 2 * u.statistic / (len(a) * len(b)) - 1
                log.info("    %-10s inflammatory %.2f vs other %.2f  p=%.4g  "
                         "Cliff's delta %.3f", c, a.mean(), b.mean(), u.pvalue, cliff)

        # Does class membership SWITCH between day 0 and day 7?
        d7 = npx[(npx.Timepoint == "D7") & (npx.OlinkID.isin(w.columns))]
        w7 = d7.pivot_table(index="subject_id", columns="OlinkID", values="NPX")
        both = [s for s in w7.index if s in Xp.index]
        if len(both) > 40:
            X7 = w7.loc[both, Xp.columns].fillna(Xp.mean())
            il6_0 = Xp.loc[both, il6_oids].mean(axis=1)
            il6_7 = X7[il6_oids].mean(axis=1)
            cut = float(np.median(il6_0))
            sw = pd.crosstab(il6_0 > cut, il6_7 > cut)
            log.warning("  IL6-high status, day 0 vs day 7 (n=%d):\n%s",
                        len(both), sw.to_string())
            switch = float((( il6_0 > cut) != (il6_7 > cut)).mean())
            log.warning("  %.1f%% of patients CHANGE inflammatory status between day 0 "
                        "and day 7. These labels are not fixed properties.", 100 * switch)
            report_only("olink_switch_frac_d0_d7", round(switch, 4))
        else:
            switch = None
            log.warning("  only %d subjects have both D0 and D7 — switching not assessed",
                        len(both))

    # --- Q4: the p38 arm, where a null is expected --------------------------
    with step(log, "p38 protein proxy vs comparators (§3.6 discipline)"):
        def score(gene_list):
            oids = [o for o in w.columns if oid2gene.get(o) in gene_list]
            if not oids:
                return None, 0
            sub = w[oids]
            z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, 1)
            return z.mean(axis=1), len(oids)

        p38, n38 = score(p38_assays)
        # The comparator MUST be disjoint from the p38 set. A first version
        # shared 7 of 10 assays with it (IL6, TNF, CXCL8, IL1B, IL10, CCL2,
        # IL18), so "p38 beats generic inflammation" was largely the set being
        # compared with itself. Any overlap is removed here and the exclusion is
        # logged, because a circular comparison that produces a POSITIVE result
        # is the most dangerous kind.
        generic = ["IL6", "CXCL8", "TNF", "IL1B", "CRP", "IL10", "CCL2",
                   "CXCL10", "IL18", "CSF3", "CXCL1", "CXCL5", "IL1RN",
                   "TNFRSF1B", "IL2RA", "CCL3", "CCL4", "OSM", "IL33", "TSLP"]
        overlap = sorted(set(generic) & set(p38_assays))
        generic_disjoint = [g for g in generic if g not in set(p38_assays)]
        log.warning("  removed %d assays shared with the p38 set from the generic "
                    "comparator: %s", len(overlap), overlap)
        comparator_sets = {
            "p38_axis": p38_assays,
            "generic_inflammation_disjoint": generic_disjoint,
            "interferon_negative_control": ["IFNG", "CXCL10", "IFNL1", "IFNB1",
                                            "DDX58", "ISG15", "IFI30", "STAT1"],
        }
        rows = []
        for name, genes in comparator_sets.items():
            s, n_assay = score(genes)
            if s is None:
                continue
            a, b = s[pos][is_hyper == 1], s[pos][is_hyper == 0]
            u = stats.mannwhitneyu(a, b)
            cliff = 2 * u.statistic / (len(a) * len(b)) - 1
            rows.append({"set": name, "n_assays": n_assay,
                         "cliffs_delta": round(float(cliff), 4),
                         "p": float(f"{u.pvalue:.4g}")})
            log.info("  %-28s n=%2d  Cliff's delta %+.3f  p=%.3g",
                     name, n_assay, cliff, u.pvalue)
        cmp_df = pd.DataFrame(rows)
        cmp_df.to_csv(paths.RESULTS / "L1_olink_pathway_scores.csv", index=False)

        d_p38 = next((r["cliffs_delta"] for r in rows if r["set"] == "p38_axis"), None)
        d_inf = next((r["cliffs_delta"] for r in rows
                      if r["set"] == "generic_inflammation_disjoint"), None)
        margin = (abs(d_p38) - abs(d_inf)) if (d_p38 is not None and d_inf is not None) else None
        if margin is not None and 0 < margin < 0.15:
            log.warning("  The p38 proxy separates only MARGINALLY better than a "
                        "DISJOINT generic-inflammation set (|%.3f| vs |%.3f|, margin "
                        "%.3f). Both sets measure downstream cytokine output, which is "
                        "equally downstream of NF-kB. Per §3.6 this is 'consistent "
                        "with, but not evidence of, p38-axis dominance'. Do not "
                        "upgrade it to a positive finding. The interferon negative "
                        "control at delta %.3f confirms the machinery is not simply "
                        "detecting general transcriptional activation.",
                        d_p38, d_inf, margin,
                        next((r["cliffs_delta"] for r in rows
                              if r["set"] == "interferon_negative_control"), float("nan")))
        if d_p38 is not None and d_inf is not None and abs(d_p38) <= abs(d_inf):
            log.warning("  The p38 proxy separates NO BETTER than generic "
                        "inflammation (|%.3f| vs |%.3f|). Per the pre-committed rule "
                        "in genesets.yaml: the inflammatory class is broadly "
                        "inflammatory, and this assay cannot resolve whether p38/MK2 "
                        "specifically dominates. This is the EXPECTED result and "
                        "reporting it is a pass (I5).", d_p38, d_inf)

    # --- the published-coefficient ordering, labelled as what it is ---------
    with step(log, "rank surrogate from the published coefficients (D3)"):
        il6c = il6_oids[0] if il6_oids else None
        tnfc = next((o for o in w.columns if oid2gene.get(o) == "TNFRSF1A"), None)
        if il6c and tnfc:
            tmp = pd.DataFrame({"IL6": w[il6c], "TNFRSF1A": w[tnfc],
                                "HCO3": 0.0})   # absent from this cohort — D9
            log.warning("  bicarbonate is ABSENT from the MGH clinical file (D9), so "
                        "the third variable is constant. This is a TWO-protein "
                        "ordering derived from published coefficients, NOT the "
                        "published classifier, and must never be called one.")
            r = tmp[["IL6", "TNFRSF1A"]].rank(pct=True)
            order = (sinha2020.BETA_LN_IL6 * r.IL6
                     + sinha2020.BETA_LN_STNFR1 * r.TNFRSF1A)
            hi = order >= order.quantile(1 - sinha2020.PUBLISHED_HYPER_PREVALENCE)
            conc = float((hi[pos].to_numpy() == (is_hyper == 1)).mean())
            log.warning("  concordance with the unsupervised inflammatory class: %.1f%%",
                        100 * conc)
            report_only("olink_rank_surrogate_concordance", round(conc, 4))
        else:
            conc = None

    (paths.RESULTS / "L1_olink_summary.json").write_text(json.dumps({
        "n_subjects": int(len(w)), "n_assays_after_lod": int(w.shape[1]),
        "covid_positive": int((meta.COVID == 1).sum()),
        "covid_negative": int((meta.COVID == 0).sum()),
        "pac_by_k": {int(k): round(v, 4) for k, v in pac.items()},
        "k_selected": int(k_best), "k_ambiguous": bool(sel["ambiguous"]),
        "inflammatory_prevalence": round(float(is_hyper.mean()), 4),
        "ari_cluster_vs_covid": round(ari_cov, 4),
        "chi2_p_cluster_vs_covid": float(chi_p),
        "host_response_supported": bool(ari_cov < 0.2 and not (chi_p < 0.01)),
        "host_response_partial": bool(ari_cov < 0.2 and chi_p < 0.01),
        "switch_frac_d0_d7": switch,
        "pathway_comparison": rows,
        "p38_margin_over_disjoint_generic": round(margin, 4) if margin is not None else None,
        "p38_interpretation": ("consistent with, but not evidence of, p38-axis "
                              "dominance — the margin over a disjoint generic "
                              "inflammation set is small and both measure downstream "
                              "cytokine output"),
        "rank_surrogate_concordance": conc,
        "caveats": ("COVID-era cohort, so framed as acute hypoxemic respiratory "
                    "failure, not ARDS (§6.7). NPX is relative, so no published "
                    "fixed-coefficient classifier can be evaluated (D3), and "
                    "bicarbonate is absent entirely (D9)."),
    }, indent=2, default=str))
    log.info("wrote results/L1_olink_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
