# Early postoperative driving pressure and prolonged mechanical ventilation after cardiac surgery

R code for a retrospective cohort study using MIMIC-IV v3.1 (development) and eICU-CRD v2.0 (external validation). No machine-learning models are used.

## Data access
Both databases require PhysioNet credentialed access and a signed data use agreement. No patient-level data are included in this repository.

## Requirements
R ≥ 4.3 (analyses were run with R 4.6.1) with: data.table, DBI, duckdb, survival, splines, cmprsk, WeightIt, MatchIt, mice, pROC, EValue, sandwich, logistf, officer, flextable, ggplot2, scales.

## Run
```bash
MIMIC_DIR=/path/to/mimic-iv-3.1 EICU_DIR=/path/to/eicu-crd-2.0 ./run_all.sh
```

| Script | Content |
|---|---|
| `R/00_extract.R` | Extracts the required tables from raw csv.gz files into a local DuckDB file |
| `R/01_cohort.R` | Cohort: adults, first ICU stay directly after open CABG and/or valve surgery, invasively ventilated on ICU arrival, plateau pressure and PEEP recorded 0–6 h; excludes ECMO, heart transplant, VAD and death ≤24 h. Exposure: median driving pressure (plateau pressure − PEEP) 0–6 h |
| `R/02_analysis.R` | Covariates with >20% missing values excluded; multiple imputation by chained equations (20 datasets) and Rubin's rules. Primary outcome ventilation >24 h: logistic models 1–3, continuous and categorical driving pressure, IPTW, overlap weights, doubly robust, propensity score matching, E-value, Firth and parsimonious models, hospital-clustered SE (eICU-CRD); secondary outcomes (reintubation, Fine–Gray time to extubation, ventilator-free days, ICU stay, mortality); sensitivity and subgroup analyses; prespecified early-stop rule on the primary analysis |
| `R/03a_extra_data.R` | MIMIC-IV ventilator mode at each plateau pressure measurement, alternative outcome definition (ventilator settings charted 25–30 h) and patients ventilated without plateau pressure |
| `R/03b_extra_sensitivity.R` | Controlled-mode exposure, harmonised outcome definition and selection-weighted analyses (pooled across imputations) |
| `R/03c_splines.R`, `R/rcs.R` | Restricted cubic splines (4 knots) with pooled Wald tests |
| `R/04_manuscript_summaries.R` | Table 1 and event counts |
| `R/05_figures.R` | Figures 1–4 (TIFF, 300 dpi, uncompressed) |
| `R/06_tables.R` | Tables 1–4 (docx) |

All outputs are written to `output/`.

Reporting follows STROBE.
