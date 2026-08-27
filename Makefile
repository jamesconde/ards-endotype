# `make all` regenerates every figure from raw data. SPEC.md C10, §8.1.
# Targets are added as each stage lands; a target that does not exist yet is
# better than a target that silently does nothing.

PY := python
export PYTHONPATH := $(CURDIR)

.PHONY: help all check scope test guard clean-figures env

help:
	@echo "make env      - create the conda environment"
	@echo "make check    - which data sources are present (no analysis)"
	@echo "make scope    - day 1: descriptive scoping -> results/L0_scoping.md"
	@echo "make test     - unit tests + repo guards"
	@echo "make guard    - integrity guards over every tracked file"
	@echo "make all      - full pipeline, raw data -> every figure"

env:
	conda env create -f environment.yml || conda env update -f environment.yml

check:
	$(PY) layer0_scoping/01_scope_all_sources.py --check-only

scope:
	$(PY) layer0_scoping/01_scope_all_sources.py

test:
	pytest -q tests/
	$(PY) tools/repo_guard.py --all

guard:
	$(PY) tools/repo_guard.py --all

# Stage targets land here as they are built:
#   layer1: L1 ingest -> QC -> discovery -> transfer -> pathway -> olink
#   layer2: duckdb -> missingness audit -> cohort -> features -> models -> transfer
#           -> ablation -> enrichment -> refinement
all: check scope
	@echo
	@echo "Pipeline stages beyond scoping are not built yet."
	@echo "See SPEC.md §7 for the day-by-day order and the cut lines."
	@exit 1

clean-figures:
	rm -f figures/*.png results/*.json
