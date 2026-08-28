#!/usr/bin/env bash
# Fetch the four open-access sources. DOWNLOADS.md §1-4.
# The MGH Olink data (§5) is NOT here — it needs a human at a registration form.
#
# Safe to re-run: -N -c make it resumable and idempotent.
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
mkdir -p data logs
LOG=logs/download.log
: > "$LOG"

say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# --- 1. Challenge 2019 (primary cohort) -------------------------------------
# 40,336 individual .psv files; no archive is offered. Skip the ancillary
# directories, keep the Reyna et al. manuscript (Tier 1 reading).
say "Challenge 2019: starting (40,336 files, expect this to be the slow one)"
wget -r -N -c -np -nH --cut-dirs=3 -nv -P data/challenge2019 \
     -X /files/challenge-2019/1.0.0/img \
     -X /files/challenge-2019/1.0.0/hackathon \
     -X /files/challenge-2019/1.0.0/sources \
     https://physionet.org/files/challenge-2019/1.0.0/ >> "$LOG" 2>&1
say "Challenge 2019: done ($(find data/challenge2019 -name '*.psv' | wc -l) .psv files)"

# --- 2. eICU-CRD Demo (the only mortality outcome in Layer 2) ---------------
say "eICU demo: starting"
wget -r -N -c -np -nH --cut-dirs=3 -nv -P data/eicu_demo \
     https://physionet.org/files/eicu-crd-demo/2.0.1/ >> "$LOG" 2>&1
say "eICU demo: done ($(find data/eicu_demo -type f | wc -l) files)"

# --- 3. MIMIC-IV Demo (smoke tests only — see DECISIONS.md D4) -------------
say "MIMIC-IV demo: starting"
wget -r -N -c -np -nH --cut-dirs=3 -nv -P data/mimic_iv_demo \
     https://physionet.org/files/mimic-iv-demo/2.2/ >> "$LOG" 2>&1
say "MIMIC-IV demo: done ($(find data/mimic_iv_demo -type f | wc -l) files)"

# --- 4. GEO series matrices (Layer 1) ---------------------------------------
say "GEO: starting (3 series matrices, ~120 MB)"
mkdir -p data/geo
for u in \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/GSE65682_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32707/matrix/GSE32707_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66890/matrix/GSE66890_series_matrix.txt.gz
do
  wget -c -nv -P data/geo "$u" >> "$LOG" 2>&1
done
say "GEO: done ($(ls data/geo | wc -l) files)"

say "ALL OPEN SOURCES COMPLETE"
du -sh data/* 2>/dev/null | tee -a "$LOG"
