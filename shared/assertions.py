"""Assertion checkpoints. SPEC.md §8.3.

The failure this exists to prevent: a SQL query silently returns 0 rows and the
pipeline continues happily for three more steps before anything looks wrong.

Every extraction step calls `checkpoint`. Every checkpoint also appends to
`results/checkpoints.jsonl`, which becomes the audit trail for the attrition
table and the writeup.
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from .paths import RESULTS

log = logging.getLogger(__name__)

CHECKPOINT_LOG = RESULTS / "checkpoints.jsonl"


class CheckpointError(AssertionError):
    """Raised when a pipeline stage produces something outside its expected shape."""


def _record(name: str, n_rows: int, ok: bool, detail: str = "") -> None:
    CHECKPOINT_LOG.parent.mkdir(parents=True, exist_ok=True)
    with CHECKPOINT_LOG.open("a") as fh:
        fh.write(json.dumps({
            "ts": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "checkpoint": name,
            "n_rows": n_rows,
            "ok": ok,
            "detail": detail,
        }) + "\n")


def checkpoint(
    df: pd.DataFrame,
    name: str,
    min_rows: int,
    max_rows: int | float = float("inf"),
    required_cols: Sequence[str] = (),
    no_nulls: Sequence[str] = (),
    unique_on: str | Sequence[str] | None = None,
) -> pd.DataFrame:
    """Assert a frame's shape, log it, and return it unchanged.

    Raises CheckpointError rather than continuing. That is the entire point:
    a wrong join must stop the pipeline, not decorate the next figure.
    """
    n = len(df)
    problems: list[str] = []

    if not (min_rows <= n <= max_rows):
        problems.append(f"got {n} rows, expected [{min_rows}, {max_rows}]")

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        problems.append(f"missing columns {missing}")

    for c in no_nulls:
        if c not in df.columns:
            problems.append(f"no_nulls column {c!r} absent")
        elif df[c].isna().any():
            problems.append(f"{df[c].isna().sum()} nulls in {c!r}")

    if unique_on is not None:
        keys = [unique_on] if isinstance(unique_on, str) else list(unique_on)
        present = [k for k in keys if k in df.columns]
        if len(present) != len(keys):
            problems.append(f"unique_on columns absent: {set(keys) - set(present)}")
        else:
            dupes = int(df.duplicated(subset=present).sum())
            if dupes:
                problems.append(f"{dupes} duplicate rows on {keys}")

    if problems:
        detail = "; ".join(problems)
        _record(name, n, ok=False, detail=detail)
        raise CheckpointError(
            f"{name}: {detail}. Do not proceed — the upstream step is wrong."
        )

    _record(name, n, ok=True)
    log.info("✓ checkpoint %s: %d rows", name, n)
    return df


def report_only(name: str, value: float | int, note: str = "") -> None:
    """A quantity with no pre-set bound that must still be logged prominently.

    Used where the spec says 'report the number, do not anchor on one' —
    FiO2 coverage (§6.2), cohort yield after the S/F filter (§4.3).
    """
    _record(name, int(value) if float(value).is_integer() else -1,
            ok=True, detail=f"value={value} {note}".strip())
    log.info("▸ %s = %s %s", name, value, note)


def attrition_step(rows: list[dict], label: str, n: int, reason: str = "") -> list[dict]:
    """Append one row to a CONSORT-style attrition trail (§4.3, results/L2_attrition.csv)."""
    prev = rows[-1]["n"] if rows else None
    rows.append({
        "step": len(rows),
        "label": label,
        "n": int(n),
        "excluded": None if prev is None else int(prev - n),
        "reason": reason,
    })
    log.info("attrition | %-44s n=%7d %s", label, n,
             "" if prev is None else f"(-{prev - n})")
    return rows


def expect_columns(df: pd.DataFrame, cols: Iterable[str], where: str) -> None:
    """Fail loudly on an assumed schema. SPEC.md §6.12: never assume column names."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise CheckpointError(
            f"{where}: expected columns {missing} not found. "
            f"Observed: {sorted(df.columns)[:40]}"
        )


def not_empty(path: Path, where: str) -> Path:
    if not path.exists():
        raise CheckpointError(f"{where}: {path} does not exist. Download it first — see DOWNLOADS.md")
    if path.is_file() and path.stat().st_size == 0:
        raise CheckpointError(f"{where}: {path} is empty")
    return path
