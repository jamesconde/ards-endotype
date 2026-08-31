# `make all` regenerates every figure from raw data. SPEC.md C10, §8.1.
# Targets are added as each stage lands; a target that does not exist yet is
# better than a target that silently does nothing.

PY := python
export PYTHONPATH := $(CURDIR)

.PHONY: help all check scope test guard clean-figures env figures paper

help:
	@echo "make env      - create the conda environment"
	@echo "make check    - which data sources are present (no analysis)"
	@echo "make scope    - day 1: descriptive scoping -> results/L0_scoping.md"
	@echo "make test     - unit tests + repo guards"
	@echo "make guard    - integrity guards over every tracked file"
	@echo "make all      - full pipeline, raw data -> every figure"
	@echo "make figures  - regenerate every figure at 300 dpi"
	@echo "make paper    - build paper.pdf from paper.tex (two pdflatex passes)"

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

figures:
	$(PY) layer2_clinical/12_figures.py
	$(PY) layer2_clinical/13_pipeline_diagram.py

# TinyTeX lives outside the default PATH. Two passes resolve refs; no bibtex,
# because the bibliography is a manual thebibliography environment.
paper: figures
	PATH=/home/james/.local/bin:$$PATH pdflatex -interaction=nonstopmode -halt-on-error -file-line-error paper.tex
	PATH=/home/james/.local/bin:$$PATH pdflatex -interaction=nonstopmode -halt-on-error -file-line-error paper.tex
	@echo "built paper.pdf"

clean-figures:
	rm -f figures/*.png results/*.json
