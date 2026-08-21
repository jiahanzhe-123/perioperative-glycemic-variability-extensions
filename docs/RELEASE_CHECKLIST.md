# Release Checklist

Status: **RELEASE CANDIDATE v1.1.0 (2026-08-22)**. The original v1.0.1
snapshot remains archived; the Phase 3B extension bundle is prepared for
repository release and archival review.

## 1. Code audit
- [x] All scripts inventoried (`docs/CODE_INVENTORY.md`)
- [x] Migration map complete, 83/83 files (`docs/FILE_MIGRATION_MAP.csv`)
- [x] No hard-coded local paths in analysis code (paths via `config/paths.yml`)
- [x] Execution order documented (`docs/ANALYSIS_WORKFLOW.md`)

## 2. Sensitive information
- [x] Automated scan passes (`make lint` → `scripts/sensitive_scan.py`)
- [x] Credential literals redacted from migrated scripts
- [x] Local PostgreSQL password rotated on 2026-08-04 for role `postgres` on
      local `inspire-pg` (`inspire_v142`, host port 5434); the new credential is
      stored only in ignored `config/paths.yml` (see `docs/PRIVACY_REVIEW.md`)
- [x] Git history of the private review repo confirmed clean at first push
      (fresh import; `main` only; no legacy branches).

## 3. Data compliance
- [x] No patient-level data in the repository
- [x] `DATA_AVAILABILITY.md` states access requirements for all three datasets
- [x] Synthetic data clearly labelled as simulated
- [ ] **[AUTHOR]** Confirm no third-party code is included without attribution
      (all analysis scripts are study-team originals; database vendor example
      SQL was not copied)

## 4. Dependency / license
- [x] `renv.lock` (R 4.5.0) and `requirements.txt` / `environment.yml` present
- [x] MIT license choice and copyright holder confirmed in `LICENSE`
- [x] Author names and contact e-mails filled in `CITATION.cff`, `README.md`,
      `SECURITY.md`, and `CODE_OF_CONDUCT.md`
- [ ] **[AUTHOR]** Add ORCID identifiers to `CITATION.cff` if available

## 5. Frozen results
- [x] `tests/expected/frozen_results.yml` populated with true frozen values
- [x] `make validate` PASS in authorized environment (35/35 checks;
      `results/qc/frozen_validation_report.csv`)
- [x] Re-run `make validate` immediately before tagging (35/35 PASS on
      2026-08-04)

## 6. Documentation
- [x] README (purpose, hierarchy, structure, install, config, run order,
      synthetic, validation, limitations)
- [x] `docs/REPRODUCIBILITY.md` (Levels 1–3)
- [x] `docs/RESULT_PROVENANCE.md`
- [x] `docs/PRIVACY_REVIEW.md`
- [x] `docs/OPEN_ISSUES.md`
- [x] Exhibit-level provenance map complete (43/43 manuscript and
      supplementary exhibits; see `docs/FINAL_OUTPUT_PROVENANCE.csv`)

## 7. CI
- [x] GitHub Actions green on the DOI-metadata commit `6d9aad9`
      (lint; tests + synthetic; docs). Runs `30915718898`, `30915721243`, and
      `30915721274`.

## 8. Phase 3B extension checks
- [x] Public aggregate Phase 3B source tables copied with logical provenance paths
- [x] Patient-level Figure 2 source values excluded from the public release tier
- [x] Public aggregate Figure 3 builder added and limited to 18 aggregate rows
- [x] Final Figure 1/Figure 3 aggregate-safe files and panel-level source map added
- [x] Release manifest and file-hash manifest prepared in the controlled workspace

## 9. Release mechanics
- [x] Create private GitHub repo first
- [ ] **[AUTHOR]** Invite a co-author to review
- [x] Tags `v1.0.0` and `v1.0.1` created after the technical and privacy checks
      were complete
- [x] Archive on Zenodo and update `CITATION.cff`, README, and the manuscript
      code-availability statement with concept DOI `10.5281/zenodo.21791846`
      (version DOI `10.5281/zenodo.21791847`)
- [x] Flip repository to public
