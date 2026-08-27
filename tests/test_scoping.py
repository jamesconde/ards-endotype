"""Exercise the scoping script against a synthetic Challenge 2019 tree.

The point is to validate the code path on day 0, so that day 1 is spent reading
what the real data says rather than debugging a file walker.
"""
from __future__ import annotations

import pandas as pd
import pytest

from shared import paths


@pytest.fixture
def scope_module(fake_challenge, monkeypatch):
    # The module name starts with a digit, so it is imported by path.
    import importlib.util
    monkeypatch.setattr(paths, "CHALLENGE", fake_challenge)
    spec = importlib.util.spec_from_file_location(
        "scope_all_sources", paths.REPO / "layer0_scoping" / "01_scope_all_sources.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod.paths, "CHALLENGE", fake_challenge)
    return mod


def test_finds_both_hospital_systems(scope_module, fake_challenge):
    by_set = scope_module.find_psv_files(fake_challenge)
    assert len(by_set["A"]) == 40
    assert len(by_set["B"]) == 30
    assert not by_set["unknown"], "every file should be assigned to a hospital system"


def test_scoping_reports_coverage_without_filtering(scope_module):
    out = scope_module.scope_challenge()
    assert out["present"] and out["n_patients"] == 70

    cov = out["critical_coverage"]
    # Vitals dense, labs sparse — the ordering the real data has.
    assert cov["O2Sat"] > cov["HCO3"] >= cov["FiO2"] - 0.5
    assert 0.0 <= cov["FiO2"] <= 1.0
    assert out["iculos"]["frac_ge_24h"] > 0


def test_detects_the_ab_coverage_gap(scope_module):
    """The fixture builds set B with worse lab coverage on purpose.

    If A and B differ in *measurement practice*, a drop in performance on B is a
    missingness artifact rather than transfer failure (§4.5d). Scoping has to
    surface that before the transfer experiment is interpreted, not after.
    """
    out = scope_module.scope_challenge()
    gaps = out["ab_coverage_gaps"]
    assert gaps["FiO2"]["A"] > gaps["FiO2"]["B"]
    assert gaps["FiO2"]["abs_diff"] > 0.10


def test_selection_check_runs(scope_module):
    out = scope_module.scope_challenge()
    sel = out["fio2_selection_check"]
    assert sel, "requiring FiO2 induces selection; the check must report on it (§6.2)"
    assert all("std_diff" in v for v in sel.values())


def test_report_writes_and_names_the_cohort_correctly(scope_module, tmp_path, monkeypatch):
    monkeypatch.setattr(scope_module.paths, "RESULTS", tmp_path)
    findings = {"challenge": scope_module.scope_challenge()}
    presence = pd.DataFrame([{"source": "challenge2019", "present": True,
                              "size_mb": 1.0, "n_files": 70, "expected": "test"}])
    out = scope_module.write_report(findings, presence)
    text = out.read_text()
    assert "transfer test or missingness test" in text
    # §6.3: the word ARDS must never be applied to this cohort.
    assert "ARDS" not in text.replace("mild-ARDS", "")
