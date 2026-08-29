#!/usr/bin/env python3
"""Figures. SPEC.md §9 — every figure must stand alone from its caption.

Aggregate statistics only; no row-level data reaches a figure (§8.4).

Usage
-----
    python layer2_clinical/12_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import paths
from shared.logging_setup import setup

log = setup("L2_figures")

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold",
    "figure.facecolor": "white", "savefig.bbox": "tight",
})
INK, ACCENT, WARN = "#2b2b2b", "#1f6feb", "#c0392b"


def fig_headline() -> None:
    """The §9 headline: screening burden depends on the site, not the model."""
    site = pd.read_csv(paths.RESULTS / "L2_enrichment_by_site.csv")
    curve = pd.read_csv(paths.RESULTS / "L2_enrichment_curve.csv")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.0))

    ax1.plot(curve.classifier_accuracy, curve.n_per_arm, "o-", color=ACCENT, lw=2)
    ax1.set_xlabel("Classifier accuracy")
    ax1.set_ylabel("Patients required per arm")
    ax1.set_title("Required trial size vs classifier accuracy")
    ax1.invert_xaxis()
    ax1.grid(alpha=.25, ls=":")
    n_unstrat = float(curve.n_per_arm.max())
    ax1.axhline(n_unstrat, color=INK, ls="--", lw=1)
    ax1.annotate("unstratified enrollment", (0.5, n_unstrat), xycoords=("axes fraction", "data"),
                 textcoords="offset points", xytext=(0, 6), fontsize=8, color=INK,
                 ha="center")

    lab = {"restricted_setB": "Coverage-matched\nmodel, hospital B",
           "naive_setB": "Blood-gas model,\nhospital B"}
    sub = site[site.scenario.isin(lab)].copy()
    sub["label"] = sub.scenario.map(lab)
    bars = ax2.bar(sub.label, sub.patients_to_screen,
                   color=[ACCENT, WARN], width=.55)
    ax2.set_yscale("log")
    ax2.set_ylabel("Patients screened to enrol 3,596  (log scale)")
    ax2.set_title("Screening burden at the second hospital")
    for b, v in zip(bars, sub.patients_to_screen):
        ax2.annotate(f"{int(v):,}", (b.get_x() + b.get_width() / 2, v),
                     ha="center", va="bottom", fontsize=9, fontweight="bold")
    ratio = sub.patients_to_screen.max() / sub.patients_to_screen.min()
    ax2.annotate(f"{ratio:.0f}x", (0.5, sub.patients_to_screen.max() * .35),
                 ha="center", fontsize=13, fontweight="bold", color=WARN)
    ax2.grid(alpha=.25, ls=":", axis="y")

    fig.suptitle("Enrichment is limited by what the hospital measures, "
                 "not by the classifier", y=1.02, fontsize=11, fontweight="bold")
    fig.savefig(paths.FIGURES / "headline_enrichment.png")
    plt.close(fig)
    log.info("wrote headline_enrichment.png")


def fig_coverage() -> None:
    """The coverage cliff between the two hospital systems."""
    d = json.load(open(paths.RESULTS / "L0_scoping.json"))
    ab = d["challenge"]["ab_coverage_gaps"]
    rows = sorted(((k, v["A"], v["B"]) for k, v in ab.items()),
                  key=lambda r: -(abs(r[1] - r[2])))[:14]
    names = [r[0] for r in rows]
    y = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    ax.barh(y - .2, [r[1] for r in rows], height=.38, color=ACCENT, label="Hospital A")
    ax.barh(y + .2, [r[2] for r in rows], height=.38, color="#8a8a8a", label="Hospital B")
    ax.set_yticks(y, names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of patients with the variable recorded")
    ax.set_title("Same variables, two hospitals: measurement practice differs")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.grid(alpha=.25, ls=":", axis="x")
    ax.annotate("bicarbonate: 97.4% vs 2.1%\nthe variable the published\nclassifier needs",
                (0.55, 0.4), fontsize=8, color=WARN, fontweight="bold")
    fig.savefig(paths.FIGURES / "coverage_by_hospital.png")
    plt.close(fig)
    log.info("wrote coverage_by_hospital.png")


def fig_degradation() -> None:
    """The C5 degradation curve."""
    c = pd.read_csv(paths.RESULTS / "L2_ablation_curve.csv")
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(c.n_features, c.kappa, "o-", color=ACCENT, lw=2, label="agreement (kappa)")
    ax.plot(c.n_features, c.auc, "s--", color="#8a8a8a", lw=1.5, label="AUC")
    ax.set_xlabel("Variables available to the model")
    ax.set_ylabel("Agreement with the reference assignment")
    ax.set_title("Degradation as variables go missing")
    ax.invert_xaxis()
    ax.grid(alpha=.25, ls=":")
    ax.legend(frameon=False, fontsize=8)
    half = float(c.kappa.max()) / 2
    ax.axhline(half, color=WARN, ls=":", lw=1)
    ax.annotate("half of full-feature agreement", (c.n_features.min(), half),
                textcoords="offset points", xytext=(6, 4), fontsize=8, color=WARN)
    fig.savefig(paths.FIGURES / "degradation_curve.png")
    plt.close(fig)
    log.info("wrote degradation_curve.png")


def fig_mortality() -> None:
    """eICU mortality replication with intervals reported wide."""
    m = pd.read_csv(paths.RESULTS / "L2_eicu_mortality.csv")
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    x = np.arange(len(m))
    lo = m.mortality - m.ci_low
    hi = m.ci_high - m.mortality
    ax.errorbar(x, m.mortality, yerr=[lo, hi], fmt="o", capsize=6,
                color=INK, markersize=9, lw=1.6)
    ax.set_xticks(x, [f"hypo-like\n(n={m.n[0]})", f"hyper-like\n(n={m.n[1]})"])
    ax.set_ylabel("Hospital mortality")
    ax.set_ylim(0, max(m.ci_high) * 1.25)
    ax.set_title("Mortality replicates in eICU (186 hospitals)")
    ax.grid(alpha=.25, ls=":", axis="y")
    for xi, r in zip(x, m.itertuples()):
        ax.annotate(f"{r.mortality:.1%}", (xi, r.ci_high),
                    textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=9, fontweight="bold")
    ax.annotate("frozen model trained on a\ndifferent dataset and schema\np < 0.0001",
                (0.02, 0.86), xycoords="axes fraction", fontsize=8, color=ACCENT)
    fig.savefig(paths.FIGURES / "eicu_mortality.png")
    plt.close(fig)
    log.info("wrote eicu_mortality.png")


if __name__ == "__main__":
    paths.FIGURES.mkdir(parents=True, exist_ok=True)
    for fn in (fig_headline, fig_coverage, fig_degradation, fig_mortality):
        try:
            fn()
        except Exception as e:  # noqa: BLE001
            log.error("%s failed: %s", fn.__name__, e)
    log.info("figures in %s", paths.FIGURES)
