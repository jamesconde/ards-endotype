#!/usr/bin/env python3
"""Pre-commit guards for the three ways this repo can quietly go wrong.

1. PHI / row-level data guard  - SPEC.md §0.2 item 2, §8.4
2. "ARDS" vocabulary guard     - SPEC.md §6.3, §6.12
3. Discovery-claim guard       - SPEC.md §6.13, I1

Run standalone:  python tools/repo_guard.py [--all]
Pre-commit runs it on staged files only.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# --- Guard 1: data files -----------------------------------------------------
DATA_SUFFIXES = {
    ".psv", ".parquet", ".feather", ".duckdb", ".h5", ".hdf5",
    ".rds", ".RData", ".npz", ".npy", ".pkl", ".joblib",
}
# results/ may hold aggregate csv/json; anything else may not hold csv at all.
CSV_ALLOWED_DIRS = {"results", "config"}
MAX_COMMITTED_BYTES = 2_000_000


# --- Guard 2: the word "ARDS" applied to the Challenge 2019 cohort -----------
# The cohort is hypoxemic respiratory failure, not ARDS (§6.3). The word is
# legitimate when discussing the literature or the Layer 1 molecular cohorts,
# so this flags proximity to Challenge-2019 tokens rather than the bare word.
CHALLENGE_TOKENS = re.compile(
    r"challenge[\s_-]?2019|physionet[\s_-]?2019|cinc|setA|setB|set_a|set_b|"
    r"training_set[AB]|SepsisLabel|ICULOS|hypoxemic_cohort",
    re.IGNORECASE,
)
ARDS_TOKEN = re.compile(r"\bARDS\b", re.IGNORECASE)

# --- Guard 3: discovery language ---------------------------------------------
DISCOVERY_PATTERNS = [
    (re.compile(r"\bwe (?:discovered|identified|found) (?:a |two |novel )?"
                r"(?:new |novel )?(?:sub)?(?:phenotypes?|endotypes?)", re.I),
     "reads as a discovery claim; these subphenotypes are established (I1)"),
    (re.compile(r"\bnovel (?:sub)?(?:phenotypes?|endotypes?)\b", re.I),
     "'novel endotype' - this project reproduces, it does not discover (I1)"),
    (re.compile(r"\bfor the first time\b", re.I),
     "priority claim (I1)"),
    (re.compile(r"\b(?:our|this) (?:discovery|novel finding)\b", re.I),
     "discovery framing (I1)"),
]

# --- Guard 4: secrets ---------------------------------------------------------
# The Discord bridge puts a bot token and a webhook URL on disk. Both are
# credentials: anyone with the webhook can post to the channel, and anyone with
# the token controls the bot. .gitignore covers secrets/, but a gitignore is one
# `git add -f` away from being bypassed, so the content is matched too.
SECRET_PATTERNS = [
    (re.compile(r"https://discord(?:app)?\.com/api/webhooks/\d+/[\w-]+"),
     "Discord webhook URL — a credential, never commit it"),
    (re.compile(r"\b[MNO][A-Za-z\d_-]{23,27}\.[A-Za-z\d_-]{6}\.[A-Za-z\d_-]{27,}\b"),
     "Discord bot token"),
    (re.compile(r"\b(?:ghp|gho|ghs|ghu)_[A-Za-z0-9]{36}\b"), "GitHub token"),
    (re.compile(r"\bsk-(?:ant-)?[A-Za-z0-9_-]{20,}\b"), "API key"),
    (re.compile(r"AKIA[0-9A-Z]{16}"), "AWS access key id"),
]

# Files exempt from the vocabulary guards: they quote the spec verbatim.
VOCAB_EXEMPT = {"SPEC.md", "PRIOR_ART.md", "tools/repo_guard.py", "CITATIONS.md"}

TEXT_SUFFIXES = {".md", ".py", ".sql", ".yaml", ".yml", ".R", ".txt", ".ipynb"}


def staged_files() -> list[Path]:
    out = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACM"],
        cwd=REPO, capture_output=True, text=True, check=False,
    ).stdout.split()
    return [REPO / f for f in out]


def all_files() -> list[Path]:
    out = subprocess.run(
        ["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=False
    ).stdout.split()
    return [REPO / f for f in out]


def check_data(path: Path, rel: str, fail: list[str]) -> None:
    if path.suffix in DATA_SUFFIXES:
        fail.append(f"{rel}: data file ({path.suffix}) must not be committed (§0.2)")
        return
    if path.suffix == ".csv" and rel.split("/")[0] not in CSV_ALLOWED_DIRS:
        fail.append(f"{rel}: .csv outside results/ or config/ (§8.4)")
    if path.exists() and path.stat().st_size > MAX_COMMITTED_BYTES:
        mb = path.stat().st_size / 1e6
        fail.append(f"{rel}: {mb:.1f} MB exceeds the 2 MB committed-file limit")


def check_vocabulary(rel: str, text: str, fail: list[str]) -> None:
    for i, line in enumerate(text.splitlines(), 1):
        if ARDS_TOKEN.search(line) and CHALLENGE_TOKENS.search(line):
            fail.append(
                f"{rel}:{i}: 'ARDS' on the same line as a Challenge-2019 token. "
                f"That cohort is hypoxemic respiratory failure (§6.3).\n    {line.strip()[:120]}"
            )


def check_secrets(rel: str, text: str, fail: list[str]) -> None:
    """Credentials must never be committed, exempt file or not."""
    for i, line in enumerate(text.splitlines(), 1):
        for pat, why in SECRET_PATTERNS:
            if pat.search(line):
                fail.append(f"{rel}:{i}: {why}. Move it to secrets/ and keep it gitignored.")


def check_discovery(rel: str, text: str, fail: list[str]) -> None:
    for i, line in enumerate(text.splitlines(), 1):
        for pat, why in DISCOVERY_PATTERNS:
            if pat.search(line):
                fail.append(f"{rel}:{i}: {why}\n    {line.strip()[:120]}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="check every tracked file")
    ap.add_argument("files", nargs="*")
    args = ap.parse_args()

    if args.files:
        paths = [Path(f).resolve() for f in args.files]
    else:
        paths = all_files() if args.all else staged_files()

    fail: list[str] = []
    for p in paths:
        try:
            rel = str(p.relative_to(REPO))
        except ValueError:
            rel = str(p)
        check_data(p, rel, fail)
        if p.suffix in TEXT_SUFFIXES and p.exists():
            text = p.read_text(errors="replace")
            check_secrets(rel, text, fail)          # no exemptions for credentials
            if rel not in VOCAB_EXEMPT:
                check_vocabulary(rel, text, fail)
                check_discovery(rel, text, fail)

    if fail:
        print("repo_guard: BLOCKED\n", file=sys.stderr)
        for f in fail:
            print(f"  ✗ {f}", file=sys.stderr)
        print(
            "\nIf a hit is a false positive, fix the wording rather than the guard.\n"
            "Genuine exemptions go in VOCAB_EXEMPT with a comment saying why.",
            file=sys.stderr,
        )
        return 1
    print(f"repo_guard: ok ({len(paths)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
