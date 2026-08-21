#!/usr/bin/env Rscript

# Phase 2C evidence-lock QC. Read-only checks are performed against the
# registry, maps, manifest, upstream gates, and declared source files. This
# script does not fit models or calculate new inferential statistics.

phase2c_arg <- commandArgs(trailingOnly = FALSE)
phase2c_file_arg <- grep("^--file=", phase2c_arg, value = TRUE)
if (length(phase2c_file_arg) == 0L) stop("Phase 2C QC must be run with Rscript")
phase2c_script_path <- normalizePath(sub("^--file=", "", phase2c_file_arg[1]), mustWork = TRUE)
phase2c_workspace <- normalizePath(file.path(dirname(phase2c_script_path), "../.."), mustWork = TRUE)
source(file.path(phase2c_workspace, "scripts", "phase2b", "phase2b_common.R"))
suppressMessages(library(jsonlite))

cfg <- phase2b_load_config()
out <- function(...) phase2b_output_path(cfg, ...)
cfg_path <- function(key) phase2a_cfg_path(cfg, key)
read_csv <- function(path) phase2b_read(path)
sha256 <- function(path) phase2b_sha256(path)
failures <- character()
passes <- character()
check <- function(condition, pass_message, fail_message = pass_message) {
  if (isTRUE(condition)) passes <<- c(passes, pass_message) else failures <<- c(failures, fail_message)
  invisible(condition)
}
same_num <- function(a, b, tol = 1e-12) {
  aa <- suppressWarnings(as.numeric(a)[1]); bb <- suppressWarnings(as.numeric(b)[1])
  isTRUE(is.finite(aa) && is.finite(bb) && abs(aa - bb) <= tol * max(1, abs(aa), abs(bb)))
}
same_text <- function(a, b) identical(as.character(a)[1], as.character(b)[1])
read_lines_safe <- function(path) if (file.exists(path)) readLines(path, warn = FALSE, encoding = "UTF-8") else character()

registry_path <- out("evidence_lock", "master_evidence_registry.csv")
locked_landscape_path <- out("evidence_lock", "analytic_context_landscape_LOCKED.csv")
evidence_map_path <- out("evidence_lock", "manuscript_evidence_map.csv")
figure_map_path <- out("evidence_lock", "figure_source_map_LOCKED.csv")
manifest_path <- out("evidence_lock", "evidence_lock_manifest.json")
brief_path <- file.path(phase2c_workspace, "PHASE3_MANUSCRIPT_REWRITE_BRIEF.md")
report_path <- file.path(phase2c_workspace, "PHASE2C_EVIDENCE_LOCK_REPORT.md")

required_files <- c(registry_path, locked_landscape_path, evidence_map_path, figure_map_path, manifest_path, brief_path, report_path)
check(all(file.exists(required_files)), "All required Phase 2C lock files exist.", paste("Missing required Phase 2C files:", paste(required_files[!file.exists(required_files)], collapse = "; ")))

registry <- read_csv(registry_path)
landscape <- read_csv(out("machine_readable", "phase2b_analytic_context_landscape.csv"))
locked_landscape <- read_csv(locked_landscape_path)
evidence_map <- read_csv(evidence_map_path)
figure_map <- read_csv(figure_map_path)
manifest <- fromJSON(manifest_path, simplifyVector = FALSE)

registry_required <- c("evidence_id", "domain", "analysis", "population", "N", "events", "estimand", "effect_type", "estimate", "lower95", "upper95", "p_value", "secondary_statistic", "secondary_value", "source_status", "source_file", "source_row_or_key", "script_origin", "new_or_original", "phase_origin", "candidate_manuscript_location", "interpretive_role", "notes")
map_required <- c("evidence_id", "priority", "destination", "figure_or_table", "main_or_supplement", "reason")
figure_required <- c("figure", "panel", "source_data_file", "analysis_id", "N", "events", "statistic", "source_status", "script", "hash")
check(identical(names(registry), registry_required), "Registry columns match the Phase 2C contract.", "Registry columns do not match the Phase 2C contract.")
check(identical(names(evidence_map), map_required), "Manuscript evidence-map columns match the contract.", "Manuscript evidence-map columns do not match the contract.")
check(identical(names(figure_map), figure_required), "Figure source-map columns match the contract.", "Figure source-map columns do not match the contract.")
check(length(unique(registry$evidence_id)) == nrow(registry) && nrow(registry) > 0L, "Registry evidence IDs are unique and non-empty.", "Registry evidence IDs are duplicated or empty.")
check(length(unique(evidence_map$evidence_id)) == nrow(evidence_map), "Every evidence-map row has a unique evidence ID.", "Evidence-map evidence IDs are duplicated.")
check(all(evidence_map$evidence_id %in% registry$evidence_id), "Every evidence-map ID exists in the master registry.", "Evidence-map contains an ID absent from the master registry.")
check(nrow(landscape) == 18L && nrow(locked_landscape) == 18L, "The analytic-context landscape has 18 rows in both source and locked copy.", "The analytic-context landscape row count is not 18 in both copies.")
check(identical(as.character(landscape$context_id), as.character(locked_landscape$context_id)), "Locked landscape preserves context-row identity and order.", "Locked landscape context-row identity/order differs from Phase 2B source.")
check(identical(sha256(out("machine_readable", "phase2b_analytic_context_landscape.csv")), sha256(locked_landscape_path)), "Locked landscape SHA-256 matches the Phase 2B source exactly.", "Locked landscape SHA-256 does not match the Phase 2B source.")

expected_families <- c("MIMIC adjustment context", "MIMIC measurement-source context", "INSPIRE timing/adjustment context", "eICU adjustment context")
check(identical(unique(as.character(landscape$context_family)), expected_families), "Landscape family order is conceptual and matches the prespecified four-family order.", "Landscape family order is not the prespecified conceptual order.")
check(all(as.character(landscape$effect_type[landscape$database %in% c("MIMIC-IV", "INSPIRE")]) == "HR"), "MIMIC and INSPIRE landscape rows are labeled HR.", "A MIMIC or INSPIRE landscape row is not labeled HR.")
check(all(as.character(landscape$effect_type[landscape$database == "eICU-CRD"]) == "RR"), "eICU landscape rows are labeled RR.", "An eICU landscape row is not labeled RR.")
check(!any(grepl("365", as.character(landscape$context_id))), "Withdrawn/365-day results are absent from the main landscape.", "A 365-day result appears in the main landscape.")

nonempty <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
check(all(nonempty(registry$source_file)), "Every registry row has a source-file field.", "At least one registry row lacks a source-file field.")
check(all(nonempty(registry$source_row_or_key)), "Every registry row has a source-row/key field.", "At least one registry row lacks a source-row/key field.")
check(all(nonempty(registry$script_origin)), "Every registry row has a script-origin field.", "At least one registry row lacks a script-origin field.")
check(all(file.exists(as.character(registry$source_file))), "Every registry source file exists at evidence lock.", "At least one registry source file is missing.")
check(all(file.exists(as.character(registry$script_origin))), "Every registry script origin exists at evidence lock.", "At least one registry script origin is missing.")
check(all(as.character(registry$source_status) %in% c("LOCKED_RESULT_VERIFIED", "LOCKED_RESULT_REPRODUCED_PHASE2A", "LOCKED_CONTEXT_SUMMARY", "EXTENSION_RESULT_LOCKED", "NOT_ESTIMABLE", "NOT_RUN")), "Registry source-status values are controlled.", "Registry contains an uncontrolled source-status value.")

row_by_id <- function(id) {
  hit <- registry[registry$evidence_id == id, , drop = FALSE]
  if (nrow(hit) != 1L) return(NULL)
  hit[1, , drop = FALSE]
}
assert_registry_source <- function(id, path, key, N, events, estimate = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_) {
  row <- row_by_id(id)
  ok <- !is.null(row) && same_text(row$source_file, path) && grepl(key, as.character(row$source_row_or_key), fixed = TRUE) && same_num(row$N, N) && same_num(row$events, events)
  if (is.finite(estimate)) ok <- ok && same_num(row$estimate, estimate)
  if (is.finite(lower)) ok <- ok && same_num(row$lower95, lower)
  if (is.finite(upper)) ok <- ok && same_num(row$upper95, upper)
  if (is.finite(p)) ok <- ok && same_num(row$p_value, p)
  check(ok, paste("Registry/source concordance:", id), paste("Registry/source concordance failed:", id))
}

mice_path <- cfg_path("mimic_mice_results")
mice <- read_csv(mice_path)
mice_b30 <- mice[mice$model_id == "MICE_B_30d", , drop = FALSE]
assert_registry_source("MIMIC_MODEL_B_30D", mice_path, "MICE_B_30d", 10561, 296, mice_b30$HR_per10, mice_b30$lo_per10, mice_b30$hi_per10, mice_b30$P_per10)
source_path <- cfg_path("mimic_source_results")
source_results <- read_csv(source_path)
for (spec in list(c("MIMIC_SAMEPATIENT_POCT_HR", "SRC_samepatient_poct_30d"), c("MIMIC_SAMEPATIENT_LAB_HR", "SRC_samepatient_lab_30d"))) {
  sr <- source_results[source_results$model_id == spec[2], , drop = FALSE]
  assert_registry_source(spec[1], source_path, spec[2], 409, 49, sr$HR_per10, sr$lo, sr$hi, sr$P)
}
inspire_path <- cfg_path("inspire_primary_results")
inspire <- read_csv(inspire_path)
ir <- inspire[inspire$model_id == "ADMINV5_I2_30d", , drop = FALSE]
assert_registry_source("INSPIRE_I2_24H_30D", inspire_path, "ADMINV5_I2_30d", 1353, 27, ir$HR_per10, ir$lo, ir$hi, ir$P)
inspire48_path <- cfg_path("inspire_48h_results")
ir48 <- read_csv(inspire48_path)
ir48r <- ir48[ir48$model_id == "ADMINV5_I2_48h_landmark_30d", , drop = FALSE]
assert_registry_source("INSPIRE_CORRECTED_48H_30D", inspire48_path, "ADMINV5_I2_48h_landmark_30d", 1511, 31, ir48r$HR_per10, ir48r$lo, ir48r$hi, ir48r$P)
eicu_path <- cfg_path("eicu_locked_results")
eicu <- read_csv(eicu_path)
e3 <- eicu[eicu$model_id == "HARM_eICU-CRD_M3", , drop = FALSE]
e4 <- eicu[eicu$model_id == "HARM_eICU-CRD_M4", , drop = FALSE]
assert_registry_source("EICU_M3", eicu_path, "HARM_eICU-CRD_M3", 7115, 130, e3$RR_per10, e3$lo, e3$hi, e3$P)
assert_registry_source("EICU_M4", eicu_path, "HARM_eICU-CRD_M4", 7115, 130, e4$RR_per10, e4$lo, e4$hi, e4$P)

agreement_path <- out("machine_readable", "phase2a_same_patient_agreement.csv")
agreement <- read_csv(agreement_path)
agr <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "GV_SD", , drop = FALSE]
registry_agr <- row_by_id("P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE")
check(nrow(agr) == 1L && !is.null(registry_agr) && same_num(registry_agr$N, 452) && same_num(registry_agr$estimate, agr$mean_paired_difference) && same_num(registry_agr$lower95, agr$loa_lower95) && same_num(registry_agr$upper95, agr$loa_upper95), "Primary same-patient agreement registry values match Phase 2A output.", "Primary same-patient agreement registry values do not match Phase 2A output.")
boot_path <- out("machine_readable", "phase2a_same_patient_bootstrap_summary.csv")
boot <- read_csv(boot_path)[1, , drop = FALSE]
boot_registry <- row_by_id("P2A_BOOTSTRAP_OBSERVED_DELTA_BETA")
check(!is.null(boot_registry) && same_num(boot_registry$N, 409) && same_num(boot_registry$events, 49) && same_num(boot_registry$estimate, boot$observed_delta_beta) && same_num(boot_registry$lower95, boot$percentile_ci_delta_beta_lower) && same_num(boot_registry$upper95, boot$percentile_ci_delta_beta_upper), "Paired bootstrap registry values match Phase 2A output.", "Paired bootstrap registry values do not match Phase 2A output.")
log_path <- out("machine_readable", "phase2b_log_scale_agreement.csv")
log_scale <- read_csv(log_path)
log_primary <- log_scale[log_scale$cohort_id == "PRIMARY_FINAL_TARGET", , drop = FALSE]
log_registry <- row_by_id("P2B_LOG_PRIMARY_GEOMETRIC_MEAN_RATIO")
check(nrow(log_primary) == 1L && !is.null(log_registry) && same_num(log_registry$N, 447) && same_num(log_registry$estimate, log_primary$geometric_mean_ratio) && same_num(as.numeric(log_registry$secondary_value), NA_real_) == FALSE, "Primary log-scale ratio is registered with valid N=447.", "Primary log-scale ratio registry check failed.")
imbalance_path <- out("machine_readable", "phase2b_source_sampling_imbalance.csv")
imbalance <- read_csv(imbalance_path)
imb_primary <- imbalance[imbalance$cohort_id == "PRIMARY_FINAL_TARGET" & imbalance$variable == "delta_count", , drop = FALSE]
count_registry <- row_by_id("P2B_COUNT_PRIMARY_GV_CORRELATION")
check(nrow(imb_primary) == 1L && !is.null(count_registry) && same_num(count_registry$estimate, imb_primary$delta_gv_delta_count_pearson_r) && grepl("Spearman", as.character(count_registry$secondary_statistic), fixed = TRUE), "Primary count-imbalance correlation is registered with Pearson/Spearman provenance.", "Primary count-imbalance correlation registry check failed.")
span_registry <- row_by_id("P2B_SPAN_PRIMARY_NOT_ESTIMABLE")
check(!is.null(span_registry) && identical(as.character(span_registry$source_status), "NOT_ESTIMABLE"), "Source-specific span is explicitly locked as NOT_ESTIMABLE.", "Source-specific span is not explicitly locked as NOT_ESTIMABLE.")

figure_hash_check <- function(hash_field) {
  parts <- strsplit(as.character(hash_field), ";", fixed = TRUE)[[1]]
  if (!length(parts) || any(!nzchar(parts))) return(FALSE)
  vapply(parts, function(part) {
    split <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (length(split) != 2L || !file.exists(split[1])) return(FALSE)
    identical(sha256(split[1]), split[2])
  }, logical(1)) |> all()
}
check(all(vapply(figure_map$hash, figure_hash_check, logical(1))), "Every figure-source-map hash matches its declared source file(s).", "At least one figure-source-map hash does not match its declared source file(s).")
figure_source_files <- unlist(strsplit(as.character(figure_map$source_data_file[figure_map$source_status != "SPECIFICATION_ONLY"]), ";", fixed = TRUE))
check(all(file.exists(figure_source_files)), "Figure source-data files exist for all generated/locked panels.", "A figure source-data file is missing.")
candidate_a <- out("figures", "candidate", "phase2b_candidate_figure_A_source_dependence.png")
candidate_b <- out("figures", "candidate", "phase2b_candidate_figure_B_analytic_context_landscape.png")
check(file.exists(candidate_a) && file.exists(candidate_b), "Candidate Figure A/B files exist for evidence-lock provenance.", "Candidate Figure A or B is missing.")
figure4 <- figure_map[figure_map$figure == "Figure 4", , drop = FALSE]
check(nrow(figure4) == 1L && grepl("MOVE_TO_SUPPLEMENT", figure4$statistic, fixed = TRUE), "Original mean-glucose/absolute-risk figure is explicitly marked MOVE_TO_SUPPLEMENT.", "Figure 4 decision is not explicitly MOVE_TO_SUPPLEMENT.")

manifest_entries <- manifest$files
manifest_paths <- vapply(manifest_entries, function(x) as.character(x$path), character(1))
manifest_hashes_ok <- vapply(manifest_entries, function(x) file.exists(as.character(x$path)) && identical(sha256(as.character(x$path)), as.character(x$sha256)), logical(1))
check(all(manifest_hashes_ok), "All manifest-listed files exist and match their recorded SHA-256 hashes.", "At least one manifest-listed file is missing or has a changed SHA-256 hash.")
check(identical(as.character(manifest$decision), "GO_PHASE3_MANUSCRIPT_REWRITE"), "Manifest decision is GO_PHASE3_MANUSCRIPT_REWRITE.", "Manifest decision is not GO_PHASE3_MANUSCRIPT_REWRITE.")
check(identical(as.character(manifest$reproducibility$new_inferential_analysis_in_phase2c)[1], "FALSE"), "Manifest records no new Phase 2C inferential analysis.", "Manifest does not record the no-new-inference boundary.")
check(isTRUE(as.logical(manifest$public_repository$status_clean)), "Public authoritative repository is clean at evidence lock.", "Public authoritative repository is not clean at evidence lock.")

phase2c_build_script <- file.path(phase2c_workspace, "scripts", "phase2c", "01_build_evidence_lock.R")
script_text <- paste(read_lines_safe(phase2c_build_script), collapse = "\n")
forbidden <- c("coxph\\s*\\(", "glm\\s*\\(", "glmnet", "poisson\\s*\\(", "boot\\s*\\(", "cor\\s*\\(", "lm\\s*\\(")
check(!any(vapply(forbidden, grepl, logical(1), x = script_text, perl = TRUE)), "Phase 2C scripts contain no model-fitting or new inferential-function calls.", "A Phase 2C script appears to contain a prohibited inferential-function call.")
check(!any(grepl("^GO_PHASE[23]", read_lines_safe(report_path)[1:(length(read_lines_safe(report_path)) - 1L)])), "Evidence-lock report contains no premature gate token before its final line.", "Evidence-lock report contains a premature gate token.")
report_lines <- read_lines_safe(report_path)
check(length(report_lines) > 0L && identical(trimws(report_lines[length(report_lines)]), "GO_PHASE3_MANUSCRIPT_REWRITE"), "Evidence-lock report ends with the required GO status.", "Evidence-lock report does not end with the required GO status.")
brief_lines <- read_lines_safe(brief_path)
check(any(grepl("Measurement-Context Dependence of Routine-Care Glycemic Variability After Cardiac Surgery: A Multidatabase Cohort Study", brief_lines, fixed = TRUE)), "Rewrite brief contains the recommended title.", "Rewrite brief is missing the recommended title.")

decision <- if (length(failures) == 0L) "GO_PHASE3_MANUSCRIPT_REWRITE" else "HOLD_PHASE3_EVIDENCE_LOCK_FAIL"
qc_md <- c(
  "# Phase 2C Evidence-Lock QC",
  "",
  paste0("Decision: **", decision, "**"),
  "",
  "No new inferential analysis, model fitting, bootstrap, correlation, agreement, or specification search was performed by this QC script.",
  "",
  "## Passed checks",
  "",
  if (length(passes)) paste0("- ", passes) else "- None",
  "",
  "## Failed checks",
  "",
  if (length(failures)) paste0("- ", failures) else "- None",
  "",
  paste0("Registry rows: ", nrow(registry), "; locked landscape rows: ", nrow(locked_landscape), "; manifest files checked: ", length(manifest_entries), "."),
  paste0("Locked landscape SHA-256: ", sha256(locked_landscape_path)),
  "",
  decision
)
qc_md_path <- out("qc", "phase2c_evidence_lock_qc.md")
writeLines(qc_md, qc_md_path, useBytes = TRUE)
gate <- list(
  phase = "Phase 2C",
  decision = decision,
  hard_failure = length(failures) > 0L,
  modeling_performed = FALSE,
  new_inferential_analysis = FALSE,
  registry_rows = nrow(registry),
  landscape_rows = nrow(locked_landscape),
  checks_passed = length(passes),
  checks_failed = length(failures),
  failures = failures,
  locked_landscape_sha256 = sha256(locked_landscape_path),
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
gate_path <- out("qc", "phase2c_final_gate.json")
write_json(gate, gate_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
cat(decision, "\n")
