# Contributing

Thank you for your interest. This repository accompanies a manuscript and its
primary purpose is auditability and reproducibility of frozen results.

## Ground rules

1. **Scientific results are frozen.** Pull requests that change cohort
   definitions, exposure windows, landmark definitions, outcome definitions,
   model specifications, or any reported number will not be merged while the
   manuscript is under review. Bug reports about code that fails to reproduce
   a frozen result are very welcome — please open an issue first.
2. **No patient-level data, ever.** Do not commit anything under `data/`
   (except `data/synthetic/`), any credentials, or local paths. `make lint`
   must pass.
3. **Keep modules independent.** MIMIC, INSPIRE, and eICU modules must remain
   runnable independently; there is intentionally no pooled cross-database
   analysis.

## Development setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# R 4.5 with packages listed in renv.lock
make lint && make test && make synthetic
```

## Pull requests

- Describe the motivation and link any issue.
- Confirm `make lint`, `make test`, and `make synthetic` pass.
- Confirm no frozen result changes (`make validate` in an authorized
  environment, if you have one — otherwise state that you did not run it).
- Do not include results generated from real patient data unless they are
  aggregate statistics already public in the manuscript.

## Code style

- Follow the existing style in each file; do not reformat unrelated code.
- R: no `setwd()`; paths only via `src/common/paths.R`.
- Python: paths only via `src/common/paths.py`.
