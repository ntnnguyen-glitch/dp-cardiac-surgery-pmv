#!/usr/bin/env bash
# Usage: MIMIC_DIR=/path/mimic-iv-3.1 EICU_DIR=/path/eicu-crd-2.0 ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p output; cd output
export DUCKDB="${DUCKDB:-$PWD/w.duckdb}"
[ -f "$DUCKDB" ] || Rscript ../R/00_extract.R
for s in 01_cohort 02_analysis 03a_extra_data 03b_extra_sensitivity 03c_splines 04_manuscript_summaries 05_figures 06_tables; do Rscript ../R/$s.R; done
