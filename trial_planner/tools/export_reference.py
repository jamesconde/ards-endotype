"""Run the REAL ards-endotypes enrichment simulation unmodified and record its unrounded values.

Produces two stamped files:
  reference/python_reference.json  every value the R port must reproduce (parity test)
  app/defaults.json                the cohort inputs the app starts from

Instrumentation only, per CLAUDE.md verification rule 1 (instrument the real path, never reimplement it):
  - recording wrappers around the module's n_per_arm() and ppv_from_accuracy(); main() looks both up at call time;
  - a sys.setprofile snapshot of main()'s locals at return (prevalence, event rates, scoreable fractions, ...);
  - side effects redirected into scratch/ards_run/: paths.RESULTS (main reads L2_transfer_results.csv there and
    writes its outputs there), logging_setup.LOGS and assertions.CHECKPOINT_LOG; DuckDB opened read-only.
Afterwards the files main() wrote are compared byte-for-byte with the committed ards-endotypes results. A mismatch
means data or code drifted since that commit, and the script refuses to write the reference.

Run from trial_planner/ with the ards env:  ~/miniconda/envs/ards/bin/python tools/export_reference.py
"""
from __future__ import annotations

import importlib.util
import inspect
import shutil
import subprocess
import sys
from pathlib import Path

ARDS = Path(__file__).resolve().parents[2]          # this repo (trial_planner/ sits at its root)
SCRIPT = ARDS / "layer2_clinical" / "09_enrichment_simulation.py"
HERE = Path(__file__).resolve().parent.parent
RUN = HERE / "scratch" / "ards_run"
OUTPUTS = ["L2_enrichment_curve.csv", "L2_enrichment_by_site.csv", "L2_enrichment_summary.json"]

sys.path.insert(0, "/home/james/molmage/src")
from molmage.util.provenance import CHECKED, dump_json, sha_of  # noqa: E402


def git(repo: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=repo, stderr=subprocess.DEVNULL).decode().strip()
    except subprocess.CalledProcessError:
        return ""


def commit_of(repo: Path) -> str:
    c = git(repo, "rev-parse", "--short", "HEAD") or "uncommitted"
    return c + ("+dirty" if git(repo, "status", "--porcelain") else "")


def main() -> int:
    (RUN / "results").mkdir(parents=True, exist_ok=True)
    (RUN / "logs").mkdir(parents=True, exist_ok=True)
    shutil.copy2(ARDS / "results" / "L2_transfer_results.csv", RUN / "results")

    sys.path.insert(0, str(ARDS))
    import duckdb
    import shared.assertions as assertions
    import shared.logging_setup as logging_setup
    from shared import paths

    paths.RESULTS = RUN / "results"
    logging_setup.LOGS = RUN / "logs"
    assertions.CHECKPOINT_LOG = RUN / "results" / "checkpoints.jsonl"
    _connect = duckdb.connect
    duckdb.connect = lambda database=":memory:", *a, **k: _connect(database, read_only=True)

    spec = importlib.util.spec_from_file_location("enrichment09", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    calls: dict[str, list] = {"n_per_arm": [], "ppv_from_accuracy": []}

    def recorder(name):
        fn = getattr(mod, name)
        sig = inspect.signature(fn)

        def wrapped(*args, **kwargs):
            out = fn(*args, **kwargs)
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            calls[name].append({"args": dict(bound.arguments), "out": out})
            return out
        return wrapped

    mod.n_per_arm = recorder("n_per_arm")
    mod.ppv_from_accuracy = recorder("ppv_from_accuracy")

    snap: dict = {}

    def profiler(frame, event, arg):
        if event == "return" and frame.f_code is mod.main.__code__:
            snap.update(frame.f_locals)

    sys.setprofile(profiler)
    try:
        rc = mod.main()
    finally:
        sys.setprofile(None)
    assert rc == 0, f"main() returned {rc}"

    same = {f: (RUN / "results" / f).read_bytes() == (ARDS / "results" / f).read_bytes() for f in OUTPUTS}
    if not all(same.values()):
        print("REFUSING: the real path no longer reproduces the committed outputs:", same)
        return 1

    inputs = {
        "09_enrichment_simulation.py": sha_of(str(SCRIPT)),
        "ards.duckdb": sha_of(str(paths.DUCKDB)),
        "scratch/features/static.parquet": sha_of(str(ARDS / "scratch" / "features" / "static.parquet")),
        "results/L2_transfer_results.csv": sha_of(str(ARDS / "results" / "L2_transfer_results.csv")),
    }
    # The provenance tool hard-codes project="MolMAGE" and reads git from the molmage repo; stamp() applies
    # **extra last, so these overrides give THIS repo's identity. #Consult with James: worth a repo= parameter?
    ident = {"project": "ards-trial-planner", "git_commit": commit_of(HERE),
             "git_branch": git(HERE, "rev-parse", "--abbrev-ref", "HEAD") or "?",
             "source_repo": str(ARDS), "source_commit": commit_of(ARDS), "inputs": inputs, "tier": CHECKED}
    observed = {k: snap[k] for k in ("prevalence", "p_hyper", "p_hypo", "p_all")}
    constants = {"alpha": mod.ALPHA, "power": mod.POWER, "rrr": mod.RRR, "acc_assumed": snap["acc_assumed"]}

    reference = {
        "constants": constants,
        "observed": observed,
        "scoreable": snap["scoreable"],
        "n_unstrat": snap["n_unstrat"],
        "n_perfect": snap["n_perfect"],
        "break_even_accuracy": snap["be"],
        "site": {"ppv": snap["ppv"], "p_enrolled": snap["p_enrolled"], "n_per_arm": snap["n"],
                 "frac_positive": snap["frac_positive"]},
        "calls": calls,
        "curve_rows": snap["rows"],
        "site_rows": snap["site_rows"],
        "reproduces_committed_outputs": same,
    }
    dump_json(reference, str(HERE / "reference" / "python_reference.json"),
              "Unrounded values from the real 09_enrichment_simulation.main(), for the R port's parity test",
              basis="float64 values captured inside main(); endpoint sepsis onset; curve uses sens = spec = accuracy",
              **ident)

    defaults = {
        "constants": constants,
        "observed": observed,
        "scoreable": snap["scoreable"],
        "labels": {"restricted_setA": "System A, coverage-matched model",
                   "naive_setA": "System A, blood-gas model",
                   "restricted_setB": "System B, coverage-matched model",
                   "naive_setB": "System B, blood-gas model"},
        "endpoint": "sepsis onset (NOT mortality, NOT ARDS-adjudicated)",
    }
    dump_json(defaults, str(HERE / "app" / "defaults.json"),
              "Cohort inputs for the trial planner, from the real 09_enrichment_simulation.main()",
              basis="unrounded float64; hypoxemic respiratory failure cohort (oximetry surrogate)", **ident)
    print("reproduces committed outputs:", same)
    print("observed:", observed)
    print("n_unstrat %.6f  n_perfect %.6f  break-even %s" % (snap["n_unstrat"], snap["n_perfect"], snap["be"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
