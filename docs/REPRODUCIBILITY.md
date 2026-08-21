# Reproducibility Guide

Three reproducibility levels, depending on data access. Nobody — including
the public — can "one-click" reproduce every result, because the three
datasets are access-restricted.

## Level 1 — no clinical data (anyone)

```bash
git clone <repo-url> && cd perioperative-glycemic-variability
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# R 4.5 + packages from renv.lock (or environment.yml for conda)
make lint        # sensitive-info scan + config checks
make test        # unit tests (landmark, dedup, GV, SHR, conversions)
make synthetic   # simulated end-to-end pipeline → results/machine_readable/synthetic/
```

This verifies that the code is executable and that core data rules behave as
documented. It does not reproduce any study estimate.

## Level 2 — MIMIC-IV / eICU-CRD authorization

Prerequisites:

- MIMIC-IV **3.1** and eICU-CRD **2.0** loaded into PostgreSQL;
- `config/paths.yml` pointing at your databases and work directories.

```bash
make mimic       # SQL extraction → glucose processing → cohort → MICE → models
make eicu        # SQL extraction → harmonized modified-Poisson models
make figures     # manuscript figures
make validate    # frozen-result verification (MIMIC/eICU sentinels)
```

Expected runtime: hours for MIMIC extraction + MICE (m=50) on a workstation.
All randomness is seeded (seed 20260726). Exact package versions are pinned in
`renv.lock` / `requirements.txt`.

## Level 3 — study-team controlled environment (INSPIRE)

The INSPIRE module additionally requires INSPIRE **1.4.2** in PostgreSQL,
available only under the dataset controllers' governance and institutional
approvals:

```bash
make inspire     # import → outcome/time repairs → v5 analysis-of-record → QC
make validate    # adds INSPIRE sentinel checks
```

## What cannot be reproduced from this repository alone

- Any result requiring INSPIRE (Level 3) — we cannot redistribute the data
  or broker access.
- The internal outcome-reconciliation audit files, which contain
  patient-level linkage and remain in the controlled environment by design.

## Environment recovery

- `renv.lock` pins R package versions (R 4.5.0).
- `requirements.txt` (+ `environment.yml`) pins Python.
- SQL assumes PostgreSQL ≥ 14; window functions and CTEs are used throughout.
- If `flexsurv` is unavailable from your CRAN mirror at the pinned version,
  install the version in `renv.lock` from the CRAN archive.

## Troubleshooting

- *"config/paths.yml not found"* — copy `config/config.example.yml`.
- *MICE loggedEvents warnings* — expected for deterministic duplicates and
  linear complements; the exclusion rules are documented in the script and
  the supplement.
- *Random-intercept eICU model convergence warning* — documented; reported as
  a sensitivity, not the primary harmonized estimate.
