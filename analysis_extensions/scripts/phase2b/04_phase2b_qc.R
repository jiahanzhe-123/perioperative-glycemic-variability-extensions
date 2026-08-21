#!/usr/bin/env Rscript

phase2b_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2b_file, mustWork = TRUE)), "phase2b_common.R"))
cfg <- phase2b_load_config()
phase2a_open_log(cfg, "phase2b_04_phase2b_qc.log")
on.exit(phase2a_close_log(), add = TRUE)

checks <- list()
add_check <- function(name, pass, detail, severity = "HARD") {
  status <- if (isTRUE(pass)) "PASS" else if (severity == "ADVISORY") "WARN" else "FAIL"
  checks[[length(checks) + 1L]] <<- list(name = name, status = status, detail = detail)
}

gate_phase2a <- tryCatch({ phase2b_require_gates(cfg); TRUE }, error = function(e) FALSE)
add_check("Phase 2A GO gate", gate_phase2a, "GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE")

log_path <- phase2b_output_path(cfg, "machine_readable", "phase2b_log_scale_agreement.csv")
imbalance_path <- phase2b_output_path(cfg, "machine_readable", "phase2b_source_sampling_imbalance.csv")
imbalance_summary_path <- phase2b_output_path(cfg, "machine_readable", "phase2b_source_sampling_imbalance_summary.csv")
landscape_path <- phase2b_output_path(cfg, "machine_readable", "phase2b_analytic_context_landscape.csv")
provenance_path <- phase2b_output_path(cfg, "machine_readable", "phase2b_result_provenance.csv")
fig_a_path <- phase2b_output_path(cfg, "figures", "candidate", "phase2b_candidate_figure_A_source_dependence.png")
fig_b_path <- phase2b_output_path(cfg, "figures", "candidate", "phase2b_candidate_figure_B_analytic_context_landscape.png")
required_outputs <- c(log_path, imbalance_path, imbalance_summary_path, landscape_path, provenance_path, fig_a_path, fig_b_path)
add_check("Phase 2B output files", all(file.exists(required_outputs)), paste(required_outputs, collapse = " | "))

log_scale <- if (file.exists(log_path)) phase2b_read(log_path) else data.frame()
log_ok <- nrow(log_scale) == 2L && all(sort(phase2b_num(log_scale$N) + phase2b_num(log_scale$positive_pairs_excluded)) == c(452, 453)) &&
  all(is.finite(phase2b_num(log_scale$mean_log_ratio))) &&
  all(is.finite(phase2b_num(log_scale$ratio_loa_lower95))) &&
  all(phase2b_num(log_scale$positive_pairs_excluded) >= 0)
add_check("Log-scale agreement closure", log_ok, if (nrow(log_scale)) paste(log_scale$cohort_id, log_scale$N, "valid; excluded", log_scale$positive_pairs_excluded, sep = " ", collapse = "; ") else "missing")

imbalance <- if (file.exists(imbalance_path)) phase2b_read(imbalance_path) else data.frame()
imbalance_summary <- if (file.exists(imbalance_summary_path)) phase2b_read(imbalance_summary_path) else data.frame()
count_ok <- nrow(imbalance) == 4L && all(imbalance$variable[imbalance$analysis_status == "COMPLETED"] == "delta_count") &&
  all(phase2b_num(imbalance$N[imbalance$variable == "delta_count"]) %in% c(452, 453))
span_boundary_ok <- nrow(imbalance) == 4L && all(imbalance$variable[imbalance$analysis_status == "NOT_ESTIMABLE"] == "delta_span") &&
  nrow(imbalance_summary) == 1L && identical(as.character(imbalance_summary$overall_status), "NOT_ESTIMABLE")
add_check("Source-count imbalance component", count_ok, "Validated n_src_poct/n_src_lab component present for primary and sensitivity cohorts")
add_check("Source-span NOT_ESTIMABLE boundary", span_boundary_ok, "No source-specific span field in validated same-patient file; no undocumented reconstruction")

landscape <- if (file.exists(landscape_path)) phase2b_read(landscape_path) else data.frame()
provenance <- if (file.exists(provenance_path)) phase2b_read(provenance_path) else data.frame()
required_columns <- c("context_id", "database", "context_family", "cohort_definition", "time_anchor", "exposure_window",
                      "source_definition", "gv_definition", "adjustment_model", "outcome", "effect_type", "N", "events",
                      "estimate", "lower95", "upper95", "p_value", "estimate_origin", "new_or_locked", "interpretive_note", "fit_status")
landscape_schema_ok <- nrow(landscape) == 18L && all(required_columns %in% names(landscape)) &&
  all(!duplicated(landscape$context_id)) && all(landscape$new_or_locked == "LOCKED") &&
  all(landscape$fit_status %in% c("LOCKED_RESULT_VERIFIED", "LOCKED_RESULT_REPRODUCED_PHASE2A"))
add_check("Landscape schema and locked status", landscape_schema_ok, paste0("rows=", nrow(landscape), "; required columns present; no duplicate context IDs"))

effect_type_ok <- nrow(landscape) > 0L &&
  all(landscape$effect_type[landscape$database %in% c("MIMIC-IV", "INSPIRE")] == "HR") &&
  all(landscape$effect_type[landscape$database == "eICU-CRD"] == "RR")
add_check("Effect-type boundaries", effect_type_ok, "MIMIC/INSPIRE=HR; eICU=RR")

denom <- function(id, n, events) {
  if (!nrow(landscape)) return(FALSE)
  rr <- landscape[landscape$context_id == id, , drop = FALSE]
  nrow(rr) == 1L && phase2b_num(rr$N) == n && phase2b_num(rr$events) == events
}
denom_ok <- denom("MIMIC_ADJUSTMENT_B_30D", 10561, 296) &&
  denom("MIMIC_SOURCE_SAME_PATIENT_POCT_30D", 409, 49) &&
  denom("MIMIC_SOURCE_SAME_PATIENT_LAB_30D", 409, 49) &&
  denom("EICU_M3_24H_AGGREGATE", 7115, 130) &&
  denom("EICU_M4_24H_AGGREGATE", 7115, 130) &&
  denom("INSPIRE_I2_24H_30D", 1353, 27)
add_check("Required final-lineage denominators", denom_ok, "MIMIC 10561/296; same-patient 409/49; eICU M3/M4 7115/130; INSPIRE I2 1353/27")

finite_ok <- nrow(landscape) == 18L && all(is.finite(phase2b_num(landscape$estimate))) &&
  all(is.finite(phase2b_num(landscape$lower95))) && all(is.finite(phase2b_num(landscape$upper95))) &&
  all(is.finite(phase2b_num(landscape$p_value)))
add_check("Landscape numeric fields finite", finite_ok, "All estimates, confidence limits, and p-values are finite")

provenance_schema_ok <- nrow(provenance) == nrow(landscape) &&
  all(provenance$context_id == landscape$context_id) &&
  all(!is.na(provenance$origin_sha256)) && all(file.exists(provenance$origin_path))
hash_ok <- if (provenance_schema_ok) all(vapply(seq_len(nrow(provenance)), function(i) {
  identical(as.character(provenance$origin_sha256[i]), phase2b_sha256(provenance$origin_path[i]))
}, logical(1))) else FALSE
add_check("Result provenance mapping", provenance_schema_ok && hash_ok, "Every landscape row maps to an existing hashed locked result")

main_text <- if (nrow(landscape)) paste(c(landscape$context_id, landscape$context_family, landscape$gv_definition, landscape$adjustment_model), collapse = " | ") else ""
forbidden_metric_hit <- grepl("(^|[^A-Za-z])(CV|ARV|MAD|IQR)([^A-Za-z]|$)", main_text, perl = TRUE)
add_check("SD-GV-only main landscape", !forbidden_metric_hit, "Main landscape contains SD-GV rows only")

withdrawn_hit <- if (nrow(landscape)) grepl("365d|withdrawn|historical|\\bold\\b", paste(landscape$context_id, landscape$outcome, collapse = " | "), ignore.case = TRUE, perl = TRUE) else FALSE
add_check("Final valid context restriction", !withdrawn_hit, "No withdrawn 365-day or historical context row in landscape")

analysis_files <- list.files(file.path(cfg$workspace_root, "scripts", "phase2b"), pattern = "\\.(R|py)$", full.names = TRUE)
analysis_files <- analysis_files[basename(analysis_files) != "04_phase2b_qc.R"]
analysis_text <- paste(vapply(analysis_files, function(p) paste(readLines(p, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
history_hit <- grepl(paste0("cardiac_gv_reanalysis_", "20260808"), analysis_text, fixed = TRUE)
add_check("Historical lineage excluded", !history_hit, "No historical extension input token in Phase 2B analysis scripts")
public_status <- system2("git", c("-C", as.character(cfg$public_code_root), "status", "--short"), stdout = TRUE, stderr = FALSE)
add_check("Authoritative public code unchanged", length(public_status) == 0L, if (length(public_status)) paste(public_status, collapse = " | ") else "clean worktree", "ADVISORY")

hard_fail <- any(vapply(checks, function(x) identical(x$status, "FAIL"), logical(1)))
decision <- if (hard_fail) "HOLD_PHASE3_STORY_REQUIRES_REVISION" else "GO_PHASE3_MANUSCRIPT_REWRITE"

check_df <- do.call(rbind, lapply(checks, function(x) {
  data.frame(check = x$name, status = x$status, detail = x$detail, stringsAsFactors = FALSE)
}))
qc_md <- c(
  "# Phase 2B landscape QC",
  "",
  paste0("Decision: **", decision, "**"),
  "",
  "The Phase 2B package is restricted to the final JAHA v5 lineage and is organized as an analytic-context landscape, not an exhaustive specification curve.",
  "",
  "| Check | Status | Detail |",
  "|---|---|---|",
  vapply(seq_len(nrow(check_df)), function(i) paste0("| ", check_df$check[i], " | ", check_df$status[i], " | ", gsub("\\|", "\\\\|", check_df$detail[i]), " |"), character(1)),
  "",
  "Interpretation boundary: source-defined GV disagreement is descriptive; source-count correlations are not causal adjustment; source-specific coefficients are not causal interactions; eICU rows are aggregate-feature based; no estimates are pooled and no cross-database heterogeneity statistic is calculated."
)
qc_path <- phase2b_output_path(cfg, "qc", "phase2b_landscape_qc.md")
dir.create(dirname(qc_path), recursive = TRUE, showWarnings = FALSE)
writeLines(qc_md, qc_path, useBytes = TRUE)

row_by_id <- function(id) {
  if (!nrow(landscape)) return(landscape)
  landscape[landscape$context_id == id, , drop = FALSE]
}
val <- function(id, col) {
  rr <- row_by_id(id)
  if (!nrow(rr) || !col %in% names(rr)) return(NA_real_)
  phase2b_num(rr[[col]])[1]
}
log_move <- function(id_after, id_before) {
  a <- val(id_after, "estimate")
  b <- val(id_before, "estimate")
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) NA_real_ else log(a) - log(b)
}

primary_agreement <- phase2b_read(phase2b_output_path(cfg, "machine_readable", "phase2a_same_patient_agreement.csv"))
primary_agreement <- primary_agreement[primary_agreement$cohort_id == "PRIMARY_FINAL_TARGET" & primary_agreement$variable == "GV_SD", , drop = FALSE]
primary_log <- log_scale[log_scale$cohort_id == "PRIMARY_FINAL_TARGET", , drop = FALSE]
count_primary <- imbalance[imbalance$cohort_id == "PRIMARY_FINAL_TARGET" & imbalance$variable == "delta_count", , drop = FALSE]
count_sens <- imbalance[imbalance$cohort_id == "SENSITIVITY_ALL_PAIRED" & imbalance$variable == "delta_count", , drop = FALSE]
source_delta <- log_move("MIMIC_SOURCE_SAME_PATIENT_POCT_30D", "MIMIC_SOURCE_SAME_PATIENT_LAB_30D")
mimic_bc_delta <- log_move("MIMIC_ADJUSTMENT_C_30D", "MIMIC_ADJUSTMENT_B_30D")
eicu_m34_delta <- log_move("EICU_M4_24H_AGGREGATE", "EICU_M3_24H_AGGREGATE")
inspire_i23_delta <- log_move("INSPIRE_I3_24H_30D", "INSPIRE_I2_24H_30D")
movement <- c(source_delta, mimic_bc_delta, eicu_m34_delta, inspire_i23_delta)
movement_names <- c("same-patient source definition", "MIMIC Model B to C", "eICU M3 to M4", "INSPIRE I2 to I3")
largest_idx <- if (all(!is.finite(movement))) NA_integer_ else which.max(abs(movement))

orig_r <- phase2b_scalar(primary_agreement$difference_pair_mean_pearson_r)
log_r <- phase2b_scalar(primary_log$log_ratio_pair_mean_pearson_r)
orig_rho <- phase2b_scalar(primary_agreement$difference_pair_mean_spearman_rho)
log_rho <- phase2b_scalar(primary_log$log_ratio_pair_mean_spearman_rho)
scale_conclusion <- if (is.finite(orig_r) && is.finite(log_r) && sign(orig_r) == sign(log_r)) {
  "The log-scale sensitivity retains the direction of the original-scale scale-dependence signal and does not overturn the Phase 2A disagreement conclusion."
} else {
  "The log-scale sensitivity changes the direction of the original-scale signal; the two scales must therefore be reported side by side without replacing the original-scale result."
}

report <- c(
  "# Phase 2B — Analytic-Context Landscape and Source-Dependence Closure",
  "",
  "## Material Passport",
  "",
  "- Phase: 2B analytic-context landscape",
  "- Lineage: JAHA v5 final lineage only",
  "- Inputs: Phase 1.6 validated context, locked final result files, final same-patient source file",
  "- Status: reproducible extension package; no manuscript or authoritative code modification",
  "- Decision: Phase 3 manuscript rewrite is permitted only under the bounded framing below",
  "",
  paste0("## Final decision: ", decision),
  "",
  "## 1. Multiplicative-scale source agreement",
  "",
  paste0("Primary log-scale comparison used N=", if (nrow(primary_log)) primary_log$N[1] else "NA", " positive paired observations. Mean log ratio=", phase2b_fmt(if (nrow(primary_log)) primary_log$mean_log_ratio[1] else NA), "; median log ratio=", phase2b_fmt(if (nrow(primary_log)) primary_log$median_log_ratio[1] else NA), "; geometric mean ratio=", phase2b_fmt(if (nrow(primary_log)) primary_log$geometric_mean_ratio[1] else NA), "."),
  paste0("95% limits on the log scale [", phase2b_fmt(if (nrow(primary_log)) primary_log$log_loa_lower95[1] else NA), ", ", phase2b_fmt(if (nrow(primary_log)) primary_log$log_loa_upper95[1] else NA), "]; exponentiated ratio limits [", phase2b_fmt(if (nrow(primary_log)) primary_log$ratio_loa_lower95[1] else NA), ", ", phase2b_fmt(if (nrow(primary_log)) primary_log$ratio_loa_upper95[1] else NA), "]."),
  paste0("Correlation of log ratio with log pair mean: Pearson r=", phase2b_fmt(log_r), "; Spearman rho=", phase2b_fmt(log_rho), ". Original-scale difference/pair-mean correlations were Pearson r=", phase2b_fmt(orig_r), " and Spearman rho=", phase2b_fmt(orig_rho), ". ", scale_conclusion),
  "",
  "The log-scale analysis is a sensitivity representation, not a replacement for the original-scale Bland–Altman analysis. It does not turn either routine-care source into a reference standard."
)

report <- c(report,
  "",
  "## 2. Sampling-count/span imbalance versus source GV disagreement",
  "",
  paste0("For the primary 452-patient cohort, POCT minus laboratory source-count difference had mean=", phase2b_fmt(if (nrow(count_primary)) count_primary$mean[1] else NA), "; SD=", phase2b_fmt(if (nrow(count_primary)) count_primary$sd[1] else NA), "; median=", phase2b_fmt(if (nrow(count_primary)) count_primary$median[1] else NA), "; IQR=", phase2b_fmt(if (nrow(count_primary)) count_primary$iqr[1] else NA), "; Q5/Q95=[", phase2b_fmt(if (nrow(count_primary)) count_primary$q5[1] else NA), ", ", phase2b_fmt(if (nrow(count_primary)) count_primary$q95[1] else NA), "]."),
  paste0("Correlation of delta GV with delta count: primary Pearson r=", phase2b_fmt(if (nrow(count_primary)) count_primary$delta_gv_delta_count_pearson_r[1] else NA), "; Spearman rho=", phase2b_fmt(if (nrow(count_primary)) count_primary$delta_gv_delta_count_spearman_rho[1] else NA), "; sensitivity Pearson r=", phase2b_fmt(if (nrow(count_sens)) count_sens$delta_gv_delta_count_pearson_r[1] else NA), "; Spearman rho=", phase2b_fmt(if (nrow(count_sens)) count_sens$delta_gv_delta_count_spearman_rho[1] else NA), "."),
  "This near-zero source-difference/count-difference relationship is distinct from the positive-but-modest patient-level GV-versus-count correlations reported in Phase 2A.",
  "Source-specific span is `NOT_ESTIMABLE`: the validated final same-patient source file contains source counts but no source-specific span fields. No undocumented span was reconstructed. Therefore the available count component does not establish that sampling imbalance explains source-GV disagreement, and the full count/span question remains incomplete rather than solved.",
  "",
  "These are descriptive source-process comparisons, not causal adjustment models and not evidence that count or span is a confounder."
)

largest_text <- if (is.na(largest_idx)) {
  "No coefficient movement could be calculated because required landscape rows were unavailable."
} else {
  paste0("The largest absolute movement among the prespecified comparisons is the ", movement_names[largest_idx], ": delta log-effect=", phase2b_fmt(movement[largest_idx]), ".")
}
report <- c(report,
  "",
  "## 3. Analytic-context landscape",
  "",
  "The master table contains 18 ordered, SD-GV-only rows grouped by context family: MIMIC adjustment (A/B/C), MIMIC measurement source (priority, POCT-only, laboratory-only, blood-gas-only, common-source, and the two same-patient rows), INSPIRE timing/adjustment (I1/I2/I3 and corrected 48-hour landmark), and eICU adjustment (M1-M4). Rows are presented in conceptual order, never sorted by effect size.",
  "",
  paste0("Same-patient outcome-model estimates were POCT HR=", phase2b_fmt(val("MIMIC_SOURCE_SAME_PATIENT_POCT_30D", "estimate")), " and laboratory HR=", phase2b_fmt(val("MIMIC_SOURCE_SAME_PATIENT_LAB_30D", "estimate")), "; their paired log-effect difference was ", phase2b_fmt(source_delta), " (ratio descriptor ", phase2b_fmt(exp(source_delta)), ")."),
  paste0("For comparison, MIMIC B→C delta log-effect=", phase2b_fmt(mimic_bc_delta), "; eICU M3→M4=", phase2b_fmt(eicu_m34_delta), "; INSPIRE I2→I3=", phase2b_fmt(inspire_i23_delta), ". ", largest_text),
  "No common summary effect, pooled estimate, or formal cross-database heterogeneity statistic was calculated. eICU estimates remain RR and MIMIC/INSPIRE estimates remain HR."
)

report <- c(report,
  "",
  "## 4. Does the evidence support measurement-context dependence as the central framing?",
  "",
  "Yes, with a bounded formulation. The evidence supports describing routine-care SD-GV as measurement-context dependent: source-defined GV shows only moderate within-patient agreement, disagreement is scale-dependent, source definition changes the estimated mortality coefficient in an identical-patient comparison, and simple count relationships are positive but modest.",
  "",
  "The evidence does not support the stronger claim that sampling intensity explains the GV association. No single observed process fully explains the heterogeneity, and source-specific differences are not causal interactions."
)

report <- c(report,
  "",
  "## 5. Main manuscript versus supplement",
  "",
  "Main manuscript candidate content:",
  "- A compact same-patient source-dependence panel: original-scale scatter/Bland–Altman, with the log-scale result identified as sensitivity rather than replacement.",
  "- The identical-patient POCT versus laboratory mortality coefficients and paired delta-beta result, explicitly labelled descriptive source-definition sensitivity.",
  "- A compact four-panel analytic-context landscape or selected context-family summary, retaining HR/RR labels and conceptual ordering.",
  "",
  "Supplement candidate content:",
  "- The full 18-row machine-readable landscape and provenance map.",
  "- Full log-scale agreement statistics, source-count imbalance distributions/correlations, the NOT_ESTIMABLE span boundary, all source-specific rows, and QC/hash evidence.",
  "- The complete candidate figures and exact locked-result selectors."
)

report <- c(report,
  "",
  "## 6. Findings that materially weaken the framing",
  "",
  "The moderate Pearson/Spearman source agreement and broad original-scale limits of agreement weaken any claim of interchangeable GV sources. The log-scale sensitivity does not erase the disagreement pattern. The validated file lacks source-specific span, so the sampling-imbalance closure is incomplete. The count relationships are modest, eICU remains aggregate-only, and the source-specific coefficient difference cannot be interpreted causally. These are reasons for a cautious measurement-context framing, not reasons to call POCT wrong or laboratory true.",
  "",
  "## 7. Final candidate title and central claim",
  "",
  "Candidate title: **Routine-Care Glycemic Variability After Cardiac Surgery: Dependence on Measurement Source, Sampling Opportunity, and Analytic Context**",
  "",
  "Central claim: **In routine-care cardiac-surgery data, SD-based glycemic variability is a measurement-context-dependent biomarker whose source, observation opportunity, and analytic specification can change its apparent mortality association, while no single observed process fully explains the differences.**",
  ""
)

report <- c(report, decision)
report_path <- file.path(cfg$workspace_root, "PHASE2B_ANALYTIC_CONTEXT_REPORT.md")
writeLines(report, report_path, useBytes = TRUE)

gate_payload <- list(
  decision = decision,
  phase = "2B",
  lineage = "JAHA_v5_final_lineage",
  hard_failure = hard_fail,
  source_span_status = if (nrow(imbalance_summary)) as.character(imbalance_summary$overall_status[1]) else "MISSING",
  landscape_rows = nrow(landscape),
  modeling_performed = FALSE,
  checks = checks,
  qc_report = qc_path,
  final_report = report_path
)
gate_path <- phase2b_output_path(cfg, "qc", "phase2b_final_gate.json")
jsonlite::write_json(gate_payload, gate_path, auto_unbox = TRUE, pretty = TRUE, na = "null")

cat("PHASE2B_DECISION", decision, "\n")
if (hard_fail) quit(status = 2L)
