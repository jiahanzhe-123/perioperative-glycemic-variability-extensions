# Privacy Review

Date: 2026-08-22 (v1.1.0 public GitHub release). Reviewer: automated scan + manual audit.

## Automated controls

- `scripts/sensitive_scan.py` (run via `make lint` and CI) scans the working
  tree for: credential patterns, PostgreSQL password literals, `/Users/…` and
  `/home/…` paths, `/Volumes/…` paths, internal server addresses, private IPv4
  ranges, SSH private-key headers, and files > 2 MB. Result at freeze: **PASS**.
- `.gitignore` blocks `data/` (except `data/synthetic/`), `results/` (except
  aggregate `results/qc/`), `config/paths.yml`, `.env*`, credential
  directories, and binary analysis artifacts (`*.rds`, `*.RData`, `*.pkl`,
  `*.parquet`, database dumps).

## Manual audit findings

1. **Original local paths in the migration map** — `docs/FILE_MIGRATION_MAP.csv`
   originally recorded source paths under the author's home directory. These
   were sanitized to `~/…` form. The unsanitized map remains only in the local
   (non-public) workspace.
2. **Database credentials in legacy scripts** — internal scripts contained a
   PostgreSQL password for the local INSPIRE instance. During migration all
   credential literals were replaced with `[REDACTED-CREDENTIAL]` and
   connection details moved to `config/paths.yml` (git-ignored). The affected
   role was confirmed as `postgres` on the local Docker instance `inspire-pg`
   (database `inspire_v142`, host port 5434); the other local containers expose
   the same built-in role but were not implicated by the legacy credential
   finding. The password was rotated on 2026-08-04 using psql's interactive
   `\password` command. A TCP connection test succeeded with the new local
   credential. The new password is stored only in the ignored local
   `config/paths.yml`; it is not recorded in this repository or report.
3. **Patient-level data** — no patient-level or stay-level records are present
   in the repository. Synthetic data in `data/synthetic/` are simulated
   (seed 20260726) and not derived from real records.
4. **Aggregate results** — `results/qc/` and `tests/expected/frozen_results.yml`
   contain only aggregate counts and effect estimates already public in the
   manuscript; the reviewed Phase 3B bundle adds only aggregate extension rows
   and logical provenance labels; no small-cell patient-level information.
5. **Git history** — this repository is a fresh import (no legacy git history).
   The legacy workspaces were never published; a history-leak scan is therefore
   not applicable to the public repo. **[AUTHOR ACTION: when creating the
   GitHub remote, verify no legacy workspace is accidentally added as a
   parent/branch.]**

## Data-sharing statement accuracy

- MIMIC-IV and eICU-CRD require PhysioNet credentialing and DUAs.
- INSPIRE is not a public dataset; access is governed by its controllers.
- This repository does not facilitate bypassing any access control.

## Residual risk

- Reviewers re-running code in their own environment must keep their own
  outputs out of public channels; the repo's `.gitignore` does not protect a
  user's separate clone after they add real data. This is called out in
  `CONTRIBUTING.md`.

## Phase 3B extension boundary

- The public Phase 3B bundle contains only aggregate landscape rows, aggregate
  source-agreement summaries, the coefficient summary, sanitized logical
  provenance paths, and an R builder for the aggregate Figure 3.
- The local Figure 2 scatter/Bland–Altman source values and final manuscript
  files remain controlled and are not copied into the repository.
- A repository-wide scan was rerun after adding the Phase 3B files; no user
  home paths, credentials, restricted inputs, patient/stay identifiers, or
  real-data bootstrap replicates were added to tracked text/source files.
