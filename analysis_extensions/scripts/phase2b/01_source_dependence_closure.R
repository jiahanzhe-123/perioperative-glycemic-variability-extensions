#!/usr/bin/env Rscript

phase2b_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2b_file, mustWork = TRUE)), "phase2b_common.R"))
cfg <- phase2b_load_config()
phase2a_open_log(cfg, "phase2b_01_source_dependence_closure.log")
on.exit(phase2a_close_log(), add = TRUE)

phase2b_require_gates(cfg)
data <- phase2a_source_data(cfg)
sp <- data$sp
target_ids <- unique(phase2a_target_ids(cfg))
primary <- sp[sp$stay_id %in% target_ids, , drop = FALSE]
candidate <- sp
if (nrow(primary) != 452L || nrow(candidate) != 453L) stop("PHASE2B_SAME_PATIENT_DENOMINATOR_FAIL")
phase2a_require_columns(sp, c("stay_id", "gv_poct", "gv_lab", "n_src_poct", "n_src_lab"), "Final same-patient source file")

log_row <- function(dd, cohort_id, cohort_definition) {
  poct <- phase2b_num(dd$gv_poct)
  lab <- phase2b_num(dd$gv_lab)
  keep <- is.finite(poct) & is.finite(lab) & poct > 0 & lab > 0
  if (sum(keep) < 2L) stop("LOG_SCALE_AGREEMENT_NOT_ESTIMABLE")
  lp <- log(poct[keep])
  ll <- log(lab[keep])
  log_ratio <- lp - ll
  log_pair_mean <- (lp + ll) / 2
  mean_log_ratio <- mean(log_ratio)
  sd_log_ratio <- sd(log_ratio)
  data.frame(
    cohort_id = cohort_id,
    cohort_definition = cohort_definition,
    N = sum(keep),
    positive_pairs_excluded = sum(!keep),
    mean_log_ratio = mean_log_ratio,
    median_log_ratio = median(log_ratio),
    geometric_mean_ratio = exp(mean_log_ratio),
    sd_log_ratio = sd_log_ratio,
    log_loa_lower95 = mean_log_ratio - 1.96 * sd_log_ratio,
    log_loa_upper95 = mean_log_ratio + 1.96 * sd_log_ratio,
    ratio_loa_lower95 = exp(mean_log_ratio - 1.96 * sd_log_ratio),
    ratio_loa_upper95 = exp(mean_log_ratio + 1.96 * sd_log_ratio),
    log_ratio_pair_mean_pearson_r = cor(log_ratio, log_pair_mean, method = "pearson"),
    log_ratio_pair_mean_spearman_rho = cor(log_ratio, log_pair_mean, method = "spearman"),
    log_scale_definition = "log(POCT GV / central-laboratory GV); pair mean = mean(log POCT GV, log laboratory GV)",
    analysis_status = "COMPLETED",
    interpretation_boundary = "multiplicative-scale descriptive agreement; original-scale Bland-Altman remains primary",
    lineage_version = "JAHA_v5_final_lineage",
    source_file = phase2a_cfg_path(cfg, "mimic_samepatient_source"),
    stringsAsFactors = FALSE
  )
}

log_scale <- do.call(rbind, list(
  log_row(primary, "PRIMARY_FINAL_TARGET", "452 paired patients in final MIMIC target"),
  log_row(candidate, "SENSITIVITY_ALL_PAIRED", "453 final source-paired candidates")
))
phase2b_write(log_scale, phase2b_output_path(cfg, "machine_readable", "phase2b_log_scale_agreement.csv"))

summary_row <- function(dd, cohort_id, cohort_definition) {
  gv <- phase2b_num(dd$gv_poct) - phase2b_num(dd$gv_lab)
  dc <- phase2b_num(dd$n_src_poct) - phase2b_num(dd$n_src_lab)
  keep <- is.finite(gv) & is.finite(dc)
  if (sum(keep) < 3L) stop("SOURCE_COUNT_IMBALANCE_NOT_ESTIMABLE")
  q <- quantile(dc[keep], c(.05, .25, .50, .75, .95), names = FALSE, type = 7)
  data.frame(
    cohort_id = cohort_id,
    cohort_definition = cohort_definition,
    N = sum(keep),
    variable = "delta_count",
    mean = mean(dc[keep]),
    sd = sd(dc[keep]),
    median = q[3],
    iqr = q[4] - q[2],
    q5 = q[1],
    q95 = q[5],
    delta_gv_delta_count_pearson_r = cor(gv[keep], dc[keep], method = "pearson"),
    delta_gv_delta_count_spearman_rho = cor(gv[keep], dc[keep], method = "spearman"),
    analysis_status = "COMPLETED",
    unavailable_reason = NA_character_,
    interpretation_boundary = "descriptive source-count imbalance; not a causal adjustment model",
    source_file = phase2a_cfg_path(cfg, "mimic_samepatient_source"),
    lineage_version = "JAHA_v5_final_lineage",
    stringsAsFactors = FALSE
  )
}

span_row <- function(dd, cohort_id, cohort_definition) {
  data.frame(
    cohort_id = cohort_id,
    cohort_definition = cohort_definition,
    N = nrow(dd),
    variable = "delta_span",
    mean = NA_real_, sd = NA_real_, median = NA_real_, iqr = NA_real_, q5 = NA_real_, q95 = NA_real_,
    delta_gv_delta_count_pearson_r = NA_real_,
    delta_gv_delta_count_spearman_rho = NA_real_,
    analysis_status = "NOT_ESTIMABLE",
    unavailable_reason = "Validated final same-patient source file provides n_src_poct and n_src_lab but no source-specific span fields; no undocumented span reconstruction performed.",
    interpretation_boundary = "source-span imbalance not estimable from declared final source file",
    source_file = phase2a_cfg_path(cfg, "mimic_samepatient_source"),
    lineage_version = "JAHA_v5_final_lineage",
    stringsAsFactors = FALSE
  )
}

imbalance <- do.call(rbind, list(
  summary_row(primary, "PRIMARY_FINAL_TARGET", "452 paired patients in final MIMIC target"),
  span_row(primary, "PRIMARY_FINAL_TARGET", "452 paired patients in final MIMIC target"),
  summary_row(candidate, "SENSITIVITY_ALL_PAIRED", "453 final source-paired candidates"),
  span_row(candidate, "SENSITIVITY_ALL_PAIRED", "453 final source-paired candidates")
))
phase2b_write(imbalance, phase2b_output_path(cfg, "machine_readable", "phase2b_source_sampling_imbalance.csv"))

imbalance_summary <- data.frame(
  analysis = "sampling imbalance versus source GV disagreement",
  primary_N = nrow(primary),
  sensitivity_N = nrow(candidate),
  delta_count_status = "COMPLETED",
  delta_span_status = "NOT_ESTIMABLE",
  overall_status = "NOT_ESTIMABLE",
  reason = "Source-specific count is available and analyzed; source-specific span is absent from the validated same-patient file and was not reconstructed.",
  interpretation_boundary = "descriptive source-process comparison; no causal adjustment or measurement-error model",
  source_file = phase2a_cfg_path(cfg, "mimic_samepatient_source"),
  lineage_version = "JAHA_v5_final_lineage",
  stringsAsFactors = FALSE
)
phase2b_write(imbalance_summary, phase2b_output_path(cfg, "machine_readable", "phase2b_source_sampling_imbalance_summary.csv"))

cat("PHASE2B_SOURCE_CLOSURE_DONE\n")
cat("LOG_SCALE_PRIMARY_GEOMETRIC_RATIO", log_scale$geometric_mean_ratio[1], "\n")
cat("SOURCE_SPAN_STATUS NOT_ESTIMABLE\n")
