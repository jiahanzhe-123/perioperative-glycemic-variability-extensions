# Perioperative glycemic variability: JAHA v5 extension release v1.1.0

> This repository combines the tracked code snapshot from the original public release with the final-lineage JAHA v5 measurement-context extension package.

> **Repository release v1.1.0 (2026-08-22).** This release adds the Phase 3B
> manuscript-production evidence package: aggregate source-dependence closure,
> the 18-row SD-GV analytic-context landscape, public aggregate figure code,
> figure-source provenance, and an explicit controlled-data boundary.

The original public code remains available at [jiahanzhe-123/perioperative-glycemic-variability](https://github.com/jiahanzhe-123/perioperative-glycemic-variability). The added package is under [`analysis_extensions/`](analysis_extensions/).

## What is public here

- The original analysis code, synthetic workflow, tests, and aggregate QC records.
- Phase 1.6 and Phase 2A-2C source code.
- Aggregate, provenance-bearing extension results and candidate figures,
  including the public Figure 3 aggregate builder and final Phase 3B source map.

No patient/stay-level records, restricted raw or controlled inputs, real-data bootstrap replicates, manuscript files, local configuration, logs, or absolute-path manifests are included. The private reviewer data package remains local.

For the extension-specific reproduction boundary and result map, see [`analysis_extensions/README.md`](analysis_extensions/README.md) and [`analysis_extensions/provenance/release_scope.csv`](analysis_extensions/provenance/release_scope.csv).

# Perioperative glucose variability after open cardiac surgery

> **Original code snapshot: PUBLIC CODE RELEASE v1.0.1 (2026-08-04).** The
> original analysis code, synthetic workflow, frozen-result checks, and release
> documentation are preserved here; the JAHA v5 extension package is documented
> above and under `analysis_extensions/`.

Analysis code for a multi-database cohort study of early postoperative glucose
variability (GV) and mortality after open cardiac surgery.

**This repository contains code only — no patient-level data.** See
[DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

## Study purpose

To estimate the conditional association between glucose variability measured
in the first 24 hours after the index time anchor and subsequent mortality in
adults undergoing open cardiac surgery, using rigorously pre-specified
landmark analyses, and to examine the joint behavior of GV and the stress
hyperglycemia ratio (SHR) in a pre-specified secondary module.

## Databases and analysis hierarchy

| Database | Version | Role |
|---|---|---|
| **MIMIC-IV** | 3.1 | **Unique primary analysis** (Model B, 30-day) |
| **INSPIRE** | 1.4.2 | Secondary exact-timestamp cohort analysis; not formal external validation |
| **eICU-CRD** | 2.0 | Multicenter harmonized comparison — **not** a formal external validation; effects are never pooled across databases |

Each database has a different time anchor and outcome estimand
(index-date-anchored calendar windows in MIMIC; exact surgical completion
timestamp in INSPIRE; ICU admission offset in eICU). They are **not** the same
estimand and results are compared descriptively, never meta-analyzed.

Nothing in this repository establishes or claims causality, treatment
targets, or clinical thresholds.

## Repository structure

```text
├── config/                 # configuration templates + machine-readable model specs
├── sql/                    # database extraction: mimic/ eicu/ inspire/
├── src/common/             # shared path/config helpers (paths.R, paths.py)
├── analyses/
│   ├── 01_cohort_construction/
│   ├── 02_glucose_processing/
│   ├── 03_primary_mimic/   # primary Model B + MICE + supplement refits
│   ├── 04_time_dependent/  # 365-day time-dependence, PH diagnostics, figures
│   ├── 05_source_sensitivity/
│   ├── 06_inspire/         # INSPIRE outcome/time repairs + v5 analysis-of-record
│   ├── 07_eicu/            # eICU harmonized comparison
│   ├── 08_shr_component/   # SHR–GV joint module (MIMIC subset)
│   └── 09_quality_control/ # assembly, QC, synthetic workflow
├── scripts/                # entry points, validation, sensitive-info scan
├── tests/                  # pytest + testthat + frozen results YAML
├── data/                   # LOCAL ONLY (git-ignored); synthetic data committed
├── results/                # LOCAL ONLY (git-ignored); aggregate QC reports kept
├── docs/                   # workflow, provenance, reproducibility, release docs
└── archive/                # superseded scripts kept for audit (not used)
```

## Software environment

- **R** 4.5.0 — package versions pinned in `renv.lock`
  (survival, mice, rms, flexsurv, sandwich, lmtest, lme4, ggplot2, dplyr,
  tidyr, readr, Hmisc, patchwork, jsonlite, broom).
- **Python** ≥ 3.11 — see `requirements.txt` (pandas, numpy, scipy,
  matplotlib, pytest, pyyaml).
- **PostgreSQL** (≥ 14 recommended) for MIMIC-IV / eICU-CRD / INSPIRE
  extraction SQL.
- A conda alternative is provided in `environment.yml`.

## Installation

```bash
git clone <repo-url> && cd perioperative-glycemic-variability
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# R packages: either
Rscript -e 'install.packages("renv"); renv::restore()'   # uses renv.lock
# or install the packages listed in renv.lock manually.
```

## Configuration

```bash
cp config/config.example.yml config/paths.yml   # git-ignored
# edit config/paths.yml: point each key at YOUR authorized local data
```

All scripts resolve paths through `src/common/paths.R` / `paths.py`. No script
contains an absolute local path. `config/constants.yml` holds every analytic
constant (windows, glucose plausibility range, seeds); each constant is
computed once in its target cohort and frozen.

## Data requirements and versions

You must independently obtain authorized access to MIMIC-IV 3.1, eICU-CRD 2.0,
and INSPIRE 1.4.2 — see [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). The data
use agreements for all three datasets prohibit redistribution, which is why no
data can be included here. INSPIRE in particular is not a public database;
access is governed by its controllers and institutional approvals.

## Run order

What you can run depends on your access level
([docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) defines Levels 1–3):

```bash
make lint        # sensitive-info scan + config/schema checks (no data needed)
make test        # unit tests: landmark, dedup, GV, SHR, conversions (no data)
make synthetic   # end-to-end workflow on simulated data (no data)
make validate    # re-check frozen results (needs authorized data + paths.yml)

make mimic       # MIMIC pipeline (needs MIMIC-IV 3.1 in PostgreSQL)
make inspire     # INSPIRE pipeline (needs INSPIRE 1.4.2 + authorization)
make eicu        # eICU pipeline (needs eICU-CRD 2.0 in PostgreSQL)
make figures     # regenerate manuscript figures from analysis outputs
make all         # lint + test + synthetic + validate
```

Detailed step-by-step order: [docs/ANALYSIS_WORKFLOW.md](docs/ANALYSIS_WORKFLOW.md).

## Synthetic example

`make synthetic` generates a fully simulated cohort (seed 20260726, no real
patients), runs glucose cleaning → feature engineering → a Cox model and a
modified-Poisson model → a figure, and writes outputs under
`results/machine_readable/synthetic/`.

> Synthetic data are provided solely to test code execution and do not
> reproduce the study estimates.

## Generating tables and figures

- Manuscript tables/figures are produced by the scripts listed per exhibit in
  [docs/RESULT_PROVENANCE.md](docs/RESULT_PROVENANCE.md) and
  [docs/FINAL_OUTPUT_PROVENANCE.csv](docs/FINAL_OUTPUT_PROVENANCE.csv). The
  final figure builders live in
  `analyses/09_quality_control/figure_build/`.
- `make figures` regenerates the figure set from existing analysis outputs.

## Verifying frozen results

`make validate` compares your regenerated outputs against the frozen sentinel
values in `tests/expected/frozen_results.yml` (sample sizes and event counts
must match exactly; effect estimates within stated tolerances). A report is
written to `results/qc/frozen_validation_report.csv`.

## What can and cannot be reproduced

- **Without any clinical data (Level 1):** lint, unit tests, synthetic
  workflow, config/model-spec checks.
- **With MIMIC/eICU authorization (Level 2):** full MIMIC and eICU pipelines,
  their tables/figures, frozen-result checks for those modules.
- **Study-team controlled environment (Level 3):** the INSPIRE module and the
  complete manuscript output set. INSPIRE's governance does not permit
  redistribution, so this module cannot be independently re-run by the public.

## Citation

See [CITATION.cff](CITATION.cff). The archived code record is
[10.5281/zenodo.21791846](https://doi.org/10.5281/zenodo.21791846); the
version-specific Zenodo record for the original `v1.0.1` snapshot is
[10.5281/zenodo.21791847](https://doi.org/10.5281/zenodo.21791847).

## License

[MIT](LICENSE). Copyright is held by the Department of Cardiology, The Second
Qilu Hospital, Cheeloo College of Medicine, Shandong University, Jinan, China.
The license covers code only, not the datasets.

## Contact

- Hanzhe Jia — `jhz0223@outlook.com`
- Hao Zhang — `zhcuriosity@163.com`
- Xin Wang — `happy97101@126.com`
- Runde Cao — `13848749698@163.com`
- Lin Deng — `denglin25@mail.sysu.edu.cn`
- Yuanyuan Sun — `syy912@126.com`
- Xiaowei Zhang — `taianrenxiaowei@126.com`
- Xiang Ning — `iningxiang@qq.com`
- Jiangying Kuang — `710887707@qq.com`

## Known limitations

- INSPIRE results can only be re-generated inside the authorized environment.
- eICU-CRD lacks reliable procedure timestamps; its module uses ICU-admission
  anchoring and is a harmonized comparison, not an external validation.
- The 365-day MIMIC analysis exhibits a proportional-hazards violation for GV;
  interval-specific estimates are reported and must not be read as protection.
- Glucose sampling is informative; measurement-process sensitivity models are
  provided and interpreted as such.
- See [docs/OPEN_ISSUES.md](docs/OPEN_ISSUES.md) for unresolved items.
