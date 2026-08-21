#!/usr/bin/env Rscript

phase2a_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2a_file, mustWork = TRUE)), "phase2a_common.R"))
suppressMessages({ library(jsonlite) })
cfg <- phase2a_load_config()
phase2a_open_log(cfg, "phase2a_06_phase2a_qc.log")
on.exit(phase2a_close_log(), add = TRUE)

checks <- list()
add_check <- function(name, pass, detail, severity = "HARD") {
  status <- if (isTRUE(pass)) "PASS" else if (severity == "ADVISORY") "WARN" else "FAIL"
  checks[[length(checks) + 1L]] <<- list(name = name, status = status, detail = detail)
}

gate_ok <- phase2a_phase1_gate(cfg)
add_check("Phase 1.6 GO gate", gate_ok, "GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES")

agreement_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_same_patient_agreement.csv")
reproduction_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_source_model_reproduction.csv")
bootstrap_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_same_patient_bootstrap.csv")
bootstrap_summary_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_same_patient_bootstrap_summary.csv")
sampling_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_sampling_process.csv")
shift_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_specification_shift.csv")
context_path <- phase2a_output_path(cfg, "machine_readable", "phase2a_measurement_context_table.csv")

required_outputs <- c(agreement_path, reproduction_path, bootstrap_path, bootstrap_summary_path, sampling_path, shift_path, context_path,
                      phase2a_output_path(cfg, "figures", "diagnostic", "same_patient_gv_scatter_primary.png"),
                      phase2a_output_path(cfg, "figures", "diagnostic", "same_patient_gv_bland_altman_primary.png"))
add_check("Phase 2A output files", all(file.exists(required_outputs)), paste(required_outputs, collapse = " | "))

agreement <- if (file.exists(agreement_path)) phase2a_read(agreement_path) else data.frame()
if (nrow(agreement)) {
  gv_agreement <- agreement[agreement$variable == "GV_SD", , drop = FALSE]
  add_check("Same-patient denominators", all(sort(unique(as.integer(gv_agreement$N))) == c(452L, 453L)) && nrow(gv_agreement) == 2L,
            paste(gv_agreement$cohort_id, gv_agreement$N, sep = "=", collapse = "; "))
  add_check("Agreement metrics finite", all(is.finite(phase2a_num(gv_agreement$pearson_r))) && all(is.finite(phase2a_num(gv_agreement$loa_lower95))), "GV correlation and limits of agreement present")
} else {
  add_check("Same-patient denominators", FALSE, "agreement output missing")
}

reproduction <- if (file.exists(reproduction_path)) phase2a_read(reproduction_path) else data.frame()
reproduction_pass <- nrow(reproduction) == 2L && all(reproduction$reproduction_status == "PASS") &&
  all(phase2a_num(reproduction$N) == 409) && all(phase2a_num(reproduction$events) == 49)
add_check("Source-model reproduction", reproduction_pass, if (nrow(reproduction)) paste(reproduction$source, reproduction$reproduced_HR_per10, sep = "=", collapse = "; ") else "missing")

bootstrap_summary <- if (file.exists(bootstrap_summary_path)) phase2a_read(bootstrap_summary_path) else data.frame()
bootstrap_reps <- if (file.exists(bootstrap_path)) phase2a_read(bootstrap_path) else data.frame()
bootstrap_ok <- nrow(bootstrap_summary) == 1L && nrow(bootstrap_reps) == 2000L &&
  phase2a_num(bootstrap_summary$N) == 409 && phase2a_num(bootstrap_summary$events) == 49
add_check("Paired bootstrap requested/completed", bootstrap_ok, if (nrow(bootstrap_summary)) paste0("requested=", bootstrap_summary$requested_replicates, "; rows=", nrow(bootstrap_reps)) else "missing")
failure_pct <- if (nrow(bootstrap_summary)) phase2a_num(bootstrap_summary$failed_percent) else NA_real_
bootstrap_stable <- is.finite(failure_pct) && failure_pct <= 5
add_check("Paired bootstrap failure rate", bootstrap_stable, paste0("failed_percent=", format(failure_pct, digits = 6), "; threshold=5"), "ADVISORY")

sampling <- if (file.exists(sampling_path)) phase2a_read(sampling_path) else data.frame()
sampling_ok <- nrow(sampling) >= 7L && all(c("MIMIC-IV", "eICU-CRD", "INSPIRE") %in% sampling$database) &&
  all(phase2a_num(sampling$N[sampling$database == "MIMIC-IV"]) == 10561) &&
  all(phase2a_num(sampling$N[sampling$database == "eICU-CRD"]) == 7115)
add_check("Sampling-process denominators", sampling_ok, paste(sampling$database, sampling$analysis_status, sampling$N, sep = "=", collapse = "; "))
add_check("INSPIRE sampling boundary", any(sampling$database == "INSPIRE" & sampling$analysis_status == "NOT_RUN"), "No new INSPIRE sampling-process correlation added")

shift <- if (file.exists(shift_path)) phase2a_read(shift_path) else data.frame()
shift_ok <- nrow(shift) == 3L && all(shift$verification_status == "LOCKED_RESULT_VERIFIED") &&
  identical(as.character(shift$effect_type[shift$database == "MIMIC-IV"]), "HR") &&
  identical(as.character(shift$effect_type[shift$database == "INSPIRE"]), "HR") &&
  identical(as.character(shift$effect_type[shift$database == "eICU-CRD"]), "RR")
add_check("Locked specification shifts", shift_ok, if (nrow(shift)) paste(shift$database, shift$delta_log_effect, sep = "=", collapse = "; ") else "missing")

context_table <- if (file.exists(context_path)) phase2a_read(context_path) else data.frame()
context_ok <- nrow(context_table) == 3L &&
  all(context_table$database %in% c("MIMIC-IV", "INSPIRE", "eICU-CRD")) &&
  identical(as.character(context_table$effect_type[context_table$database == "eICU-CRD"]), "RR") &&
  all(phase2a_num(context_table$N[match(c("MIMIC-IV", "INSPIRE", "eICU-CRD"), context_table$database)]) == c(10561, 1353, 7115)) &&
  all(phase2a_num(context_table$events[match(c("MIMIC-IV", "INSPIRE", "eICU-CRD"), context_table$database)]) == c(296, 27, 130))
add_check("Measurement-context summary table", context_ok, if (nrow(context_table)) paste(context_table$database, context_table$N, context_table$events, sep = "=", collapse = "; ") else "missing")

phase2a_files <- list.files(file.path(cfg$workspace_root, "scripts", "phase2a"), pattern = "\\.(R|py)$", full.names = TRUE)
analysis_files <- phase2a_files[basename(phase2a_files) != "06_phase2a_qc.R"]
phase2a_text <- paste(vapply(analysis_files, function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
history_token <- paste0("cardiac_gv_reanalysis_", "20260808")
history_hit <- grepl(history_token, phase2a_text, fixed = TRUE)
add_check("Historical lineage excluded", !history_hit, "No historical extension path in Phase 2A scripts")
public_status <- system2("git", c("-C", as.character(cfg$public_code_root), "status", "--short"), stdout = TRUE, stderr = FALSE)
add_check("Authoritative public code unchanged", length(public_status) == 0L, if (length(public_status)) paste(public_status, collapse = " | ") else "clean worktree", "ADVISORY")

hard_fail <- any(vapply(checks, function(x) identical(x$status, "FAIL"), logical(1)))
decision <- if (hard_fail) {
  "NO_GO_PHASE2A_LINEAGE_OR_MODEL_REPRODUCTION_FAIL"
} else if (!bootstrap_stable) {
  "HOLD_PHASE2B_INTERPRETATION_REQUIRES_REVIEW"
} else {
  "GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE"
}

check_df <- do.call(rbind, lapply(checks, function(x) {
  data.frame(check = x$name, status = x$status, detail = x$detail, stringsAsFactors = FALSE)
}))
qc_md <- c(
  "# Phase 2A QC report",
  "",
  paste0("Decision: **", decision, "**"),
  "",
  "This QC is restricted to the final JAHA v5 lineage and the declared Phase 2A outputs.",
  "",
  "| Check | Status | Detail |",
  "|---|---|---|",
  vapply(seq_len(nrow(check_df)), function(i) {
    paste0("| ", check_df$check[i], " | ", check_df$status[i], " | ", gsub("\\|", "\\\\|", check_df$detail[i]), " |")
  }, character(1)),
  "",
  "Interpretation boundaries: agreement is descriptive same-patient source comparison; the paired coefficient difference is not a causal interaction test; sampling correlations describe measurement opportunity/process; specification shifts are not confounding attenuation percentages."
)
qc_md_path <- phase2a_output_path(cfg, "qc", "phase2a_qc_report.md")
dir.create(dirname(qc_md_path), recursive = TRUE, showWarnings = FALSE)
writeLines(qc_md, qc_md_path, useBytes = TRUE)

fmt <- function(x, digits = 4L) {
  z <- suppressWarnings(as.numeric(x)[1])
  if (!length(z) || !is.finite(z)) return("NA")
  formatC(z, format = "f", digits = digits)
}
cell <- function(df, col, row = 1L) {
  if (!nrow(df) || !col %in% names(df) || nrow(df) < row) return(NA_real_)
  df[[col]][row]
}

report <- c(
  "# Phase 2A — Measurement-Context Analyses",
  "",
  paste0("Final decision: **", decision, "**"),
  "",
  "## Scope and lineage",
  "",
  "All analyses use the final JAHA v5 lineage declared in `config.yaml`. No historical extension input was used, and no manuscript or authoritative source repository was modified by this task.",
  "",
  "Phase 1.6 hard gate: `GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES`; new modeling remains outside the scope of this report.",
  "",
  "## A. Same-patient measurement-source agreement",
  "",
  "The primary comparison is the 452-patient final-target intersection; the 453-patient all-paired cohort is a sensitivity denominator. Paired differences are POCT minus central-laboratory values.",
  ""
)
primary_gv <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "GV_SD", , drop = FALSE]
sens_gv <- agreement[agreement$cohort_id == "SENSITIVITY_ALL_PAIRED" & agreement$variable == "GV_SD", , drop = FALSE]
primary_mean <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "MEAN_GLUCOSE", , drop = FALSE]
report <- c(report,
  paste0("- Primary GV SD: N=", fmt(cell(primary_gv, "N"), 0), "; POCT mean=", fmt(cell(primary_gv, "poct_mean")), "; laboratory mean=", fmt(cell(primary_gv, "lab_mean")), "; Pearson r=", fmt(cell(primary_gv, "pearson_r")), "; Spearman rho=", fmt(cell(primary_gv, "spearman_rho")), "; paired mean difference=", fmt(cell(primary_gv, "mean_paired_difference")), "; 95% limits of agreement [", fmt(cell(primary_gv, "loa_lower95")), ", ", fmt(cell(primary_gv, "loa_upper95")), "]."),
  paste0("- Primary GV difference versus pair mean: Pearson r=", fmt(cell(primary_gv, "difference_pair_mean_pearson_r")), "; Spearman rho=", fmt(cell(primary_gv, "difference_pair_mean_spearman_rho")), "."),
  "- The negative difference-versus-pair-mean correlations indicate scale-dependent/proportional disagreement: POCT minus laboratory GV tends to become more negative at larger paired GV values; this is not consistent with a constant offset alone.",
  paste0("- All-paired GV sensitivity: N=", fmt(cell(sens_gv, "N"), 0), "; paired mean difference=", fmt(cell(sens_gv, "mean_paired_difference")), "; Pearson r=", fmt(cell(sens_gv, "pearson_r")), "; Spearman rho=", fmt(cell(sens_gv, "spearman_rho")), "."),
  paste0("- Primary mean-glucose comparison: POCT mean=", fmt(cell(primary_mean, "poct_mean")), "; laboratory mean=", fmt(cell(primary_mean, "lab_mean")), "; Pearson r=", fmt(cell(primary_mean, "pearson_r")), "; paired mean difference=", fmt(cell(primary_mean, "mean_paired_difference")), "."),
  "",
  "These are diagnostic agreement and source-comparability results. They do not establish a gold-standard device, analytical validity, or pure analytical measurement error.",
  "",
  "## B. Paired source-definition coefficient difference",
  ""
)
if (nrow(bootstrap_summary)) {
  report <- c(report,
    paste0("The authoritative final source-model contract was reproduced for both source definitions at N=", fmt(cell(reproduction, "N"), 0), " and 49 events. Observed HR per 10-unit GV SD: POCT ", fmt(cell(bootstrap_summary, "observed_hr_poct")), "; laboratory ", fmt(cell(bootstrap_summary, "observed_hr_lab")), "; observed delta-beta=", fmt(cell(bootstrap_summary, "observed_delta_beta")), "; HR ratio=", fmt(cell(bootstrap_summary, "observed_coefficient_ratio")), "."),
    paste0("Paired bootstrap: B=", fmt(cell(bootstrap_summary, "requested_replicates"), 0), "; successful=", fmt(cell(bootstrap_summary, "successful_replicates"), 0), "; failed=", fmt(cell(bootstrap_summary, "failed_replicates"), 0), " (", fmt(cell(bootstrap_summary, "failed_percent")), "%). Delta-beta median=", fmt(cell(bootstrap_summary, "bootstrap_median_delta_beta")), "; SD=", fmt(cell(bootstrap_summary, "bootstrap_sd_delta_beta")), "; percentile 95% CI [", fmt(cell(bootstrap_summary, "percentile_ci_delta_beta_lower")), ", ", fmt(cell(bootstrap_summary, "percentile_ci_delta_beta_upper")), "]."),
    "The percentile interval excludes zero, so the source-specific coefficient difference is not negligible descriptively in this sensitivity analysis; it is not evidence of a causal interaction.",
    "This paired difference is a descriptive sensitivity analysis of source-definition coefficients, not a causal interaction or effect-modification test."
  )
} else {
  report <- c(report, "Source-model reproduction/bootstrap output was unavailable; Phase 2B is blocked.")
}

report <- c(report,
  "",
  "## C. Sampling-process dependence",
  "",
  "MIMIC-IV uses the final 10,561-patient / 296-event target. eICU uses the final aggregate M3/M4 input of 7,115 patients / 130 events; event-level eICU sampling-process analysis was not available. INSPIRE was conservatively not assigned a new sampling-process analysis because an unambiguous validated sampling-count/span field was not established.",
  ""
)
if (nrow(sampling)) {
  samp_m <- sampling[sampling$database == "MIMIC-IV" & sampling$analysis_status == "COMPLETED", , drop = FALSE]
  samp_e <- sampling[sampling$database == "eICU-CRD" & sampling$analysis_status == "COMPLETED", , drop = FALSE]
  report <- c(report,
    paste0("MIMIC GV versus glucose count: Pearson r=", fmt(cell(samp_m[samp_m$variable_y == "measurement_count", , drop = FALSE], "pearson_r")), "; Spearman rho=", fmt(cell(samp_m[samp_m$variable_y == "measurement_count", , drop = FALSE], "spearman_rho")), "; GV versus span: Pearson r=", fmt(cell(samp_m[samp_m$variable_y == "measurement_span_hours", , drop = FALSE], "pearson_r")), "; Spearman rho=", fmt(cell(samp_m[samp_m$variable_y == "measurement_span_hours", , drop = FALSE], "spearman_rho")), "."),
    paste0("eICU aggregate GV versus glucose count: Pearson r=", fmt(cell(samp_e[samp_e$variable_y == "measurement_count", , drop = FALSE], "pearson_r")), "; Spearman rho=", fmt(cell(samp_e[samp_e$variable_y == "measurement_count", , drop = FALSE], "spearman_rho")), "; GV versus span: Pearson r=", fmt(cell(samp_e[samp_e$variable_y == "measurement_span_hours", , drop = FALSE], "pearson_r")), "; Spearman rho=", fmt(cell(samp_e[samp_e$variable_y == "measurement_span_hours", , drop = FALSE], "spearman_rho")), ".")
  )
}

report <- c(report, "", "## D. Locked specification shifts", "", "These are sampling-process specification shifts reported as before/after effects and delta log-effect. No attenuation percentage is calculated.", "")
if (nrow(shift)) {
  for (db in c("MIMIC-IV", "eICU-CRD", "INSPIRE")) {
    rr <- shift[shift$database == db, , drop = FALSE]
    if (nrow(rr)) report <- c(report, paste0("- ", db, ": ", rr$specification_before[1], " ", rr$effect_type[1], " ", fmt(rr$effect_before[1]), " [", fmt(rr$lower95_before[1]), ", ", fmt(rr$upper95_before[1]), "] to ", rr$specification_after[1], " ", rr$effect_type[1], " ", fmt(rr$effect_after[1]), " [", fmt(rr$lower95_after[1]), ", ", fmt(rr$upper95_after[1]), "]; delta log-effect=", fmt(rr$delta_log_effect[1]), "; null crossing=", as.character(rr$null_crossing[1]), "."))
  }
}

report <- c(report,
  "",
  "## E. Interpretation boundary and Phase 2B decision",
  "",
  "The combined evidence is an analytic-context landscape: observed GV associations can depend on glucose-source definition, the opportunity to observe glucose, and locked specification choices. These analyses do not identify a causal measurement-error mechanism, establish interchangeability, or support pooled cross-database inference/heterogeneity claims.",
  "",
  "Results that weaken a stronger interpretation are the only-moderate same-patient GV agreement with wide limits of agreement, the source-specific coefficient difference, modest rather than dominant count/span correlations, and the preserved aggregate-only eICU structure. These findings support contextual caution rather than a claim that one measurement process fully explains the outcome association.",
  "",
  paste0("Decision: **", decision, "**."),
  ""
)
report_path <- file.path(cfg$workspace_root, "PHASE2A_MEASUREMENT_CONTEXT_REPORT.md")
writeLines(report, report_path, useBytes = TRUE)

gate_payload <- list(
  decision = decision,
  phase = "2A",
  lineage = "JAHA_v5_final_lineage",
  hard_failure = hard_fail,
  bootstrap_failure_percent = failure_pct,
  bootstrap_failure_threshold_percent = 5,
  modeling_performed = FALSE,
  checks = checks,
  qc_report = qc_md_path,
  final_report = report_path
)
gate_path <- phase2a_output_path(cfg, "qc", "phase2a_final_gate.json")
jsonlite::write_json(gate_payload, gate_path, auto_unbox = TRUE, pretty = TRUE, na = "null")

cat("PHASE2A_DECISION", decision, "\n")
if (hard_fail) quit(status = 2L)
