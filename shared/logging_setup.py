"""Structured logging. SPEC.md §8.2: every step logs in/out row counts and wall time."""
from __future__ import annotations

import logging
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .paths import LOGS

FMT = "%(asctime)s %(levelname)-7s %(name)-28s %(message)s"


def setup(script_name: str, level: int = logging.INFO) -> logging.Logger:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    logfile = LOGS / f"{script_name}_{stamp}.log"
    logging.basicConfig(
        level=level,
        format=FMT,
        handlers=[logging.FileHandler(logfile), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    log = logging.getLogger(script_name)
    log.info("=" * 78)
    log.info("start %s  ->  %s", script_name, logfile)
    return log


@contextmanager
def step(log: logging.Logger, name: str, n_in: int | None = None):
    """Time a stage and report rows in / rows out."""
    t0 = time.perf_counter()
    log.info("→ %s%s", name, f" (in: {n_in:,} rows)" if n_in is not None else "")
    result: dict = {}
    try:
        yield result
    finally:
        dt = time.perf_counter() - t0
        n_out = result.get("n_out")
        log.info("← %s  %.2fs%s", name, dt,
                 f"  (out: {n_out:,} rows)" if n_out is not None else "")


def set_seeds(seed: int = 42) -> int:
    """Fix and log every seed. SPEC.md §8.1."""
    import random
    import numpy as np
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass
    logging.getLogger(__name__).info("seed = %d", seed)
    return seed
