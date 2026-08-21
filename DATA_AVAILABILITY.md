# Data Availability

This repository contains **analysis code only — no patient-level data**.

## Databases and versions

| Database | Version | Role in the study | Access |
|---|---|---|---|
| MIMIC-IV | 3.1 | Unique primary analysis | PhysioNet credentialed access + signed DUA |
| eICU-CRD | 2.0 | Multicenter harmonized comparison (NOT a formal external validation) | PhysioNet credentialed access + signed DUA |
| INSPIRE | 1.4.2 | Secondary exact-timestamp cohort analysis; not formal external validation | Access governed by the INSPIRE data controllers and institutional governance; it is not a public download |

## How to obtain the data

- **MIMIC-IV**: complete the required training and apply on PhysioNet
  (https://physionet.org/content/mimiciv/). Place the extracted PostgreSQL
  database where your local `config/paths.yml` points.
- **eICU-CRD**: same process (https://physionet.org/content/eicu-crd/).
- **INSPIRE**: access is not handled by the authors of this repository. Contact
  the dataset controllers; we cannot redistribute any part of it.

## What you must provide locally

After obtaining authorized access, copy `config/config.example.yml` to
`config/paths.yml` (git-ignored) and edit the paths to your local databases
and working directories. All analysis code reads paths only from this file
(see `src/common/paths.R` / `src/common/paths.py`).

## What may be shared publicly

- Code, configuration templates, data dictionaries, and schema inventories.
- Aggregate statistics already published in the manuscript and supplement.
- Synthetic data in `data/synthetic/` (simulated; contains no real patients).

## What must NOT be shared

- Any patient-level or stay-level records, including intermediate analysis
  frames, imputed datasets, and bootstrap outputs derived from real data.
- Database credentials or local connection strings.
- Data dictionaries or materials whose redistribution is restricted by the
  source dataset's DUA.

This project does not support or condone circumventing database access
controls. If you cannot obtain access, you can still read all code and run the
synthetic workflow (`make synthetic`) and all unit tests (`make test`).
