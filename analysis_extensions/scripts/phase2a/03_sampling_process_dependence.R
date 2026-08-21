#!/usr/bin/env Rscript

phase2a_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2a_file, mustWork = TRUE)), "phase2a_common.R"))
suppressMessages({ library(jsonlite) })
cfg <- phase2a_load_config()
phase2a_open_log(cfg, "phase2a_03_sampling_process_dependence.log")
on.exit(phase2a_close_log(), add = TRUE)

if (!phase2a_phase1_gate(cfg)) stop("PHASE1_6_GATE_FAIL")
context <- phase2a_read(phase2a_output_path(cfg, "machine_readable", "final_measurement_context.csv"))
phase2a_require_columns(context, c("database", "analysis_context", "gv_sd", "mean_glucose", "measurement_count", "measurement_span_hours"), "Phase 1.6 context")

correlation_row <- function(dd, database, context_id, x_name, y_name, measurement_structure) {
  x <- phase2a_num(dd[[x_name]])
  y <- phase2a_num(dd[[y_name]])
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) stop("Insufficient complete rows for ", database, " ", x_name, " vs ", y_name)
  data.frame(
    database = database,
    analysis_context = context_id,
    variable_x = "GV_SD",
    variable_y = y_name,
    N = sum(keep),
    pearson_r = cor(x[keep], y[keep], method = "pearson"),
    spearman_rho = cor(x[keep], y[keep], method = "spearman"),
    measurement_structure = measurement_structure,
    analysis_status = "COMPLETED",
    interpretation_boundary = "descriptive correlation; not confounding or causation",
    stringsAsFactors = FALSE
  )
}

rows <- list()
mimic <- context[context$database == "MIMIC-IV" & context$analysis_context == "MIMIC_primary_24h_day1_landmark", , drop = FALSE]
eicu <- context[context$database == "eICU-CRD" & context$analysis_context == "eICU_harmonized_M4", , drop = FALSE]
rows[[length(rows) + 1L]] <- correlation_row(mimic, "MIMIC-IV", "MIMIC_primary_24h_day1_landmark", "gv_sd", "measurement_count", "final source-priority series; patient-level count")
rows[[length(rows) + 1L]] <- correlation_row(mimic, "MIMIC-IV", "MIMIC_primary_24h_day1_landmark", "gv_sd", "measurement_span_hours", "final source-priority series; patient-level span")
rows[[length(rows) + 1L]] <- correlation_row(mimic, "MIMIC-IV", "MIMIC_primary_24h_day1_landmark", "gv_sd", "mean_glucose", "final source-priority series; optional mean-glucose benchmark")
rows[[length(rows) + 1L]] <- correlation_row(eicu, "eICU-CRD", "eICU_harmonized_M4", "gv_sd", "measurement_count", "preserved aggregate all-source feature; aggregate count")
rows[[length(rows) + 1L]] <- correlation_row(eicu, "eICU-CRD", "eICU_harmonized_M4", "gv_sd", "measurement_span_hours", "preserved aggregate all-source feature; aggregate span")
rows[[length(rows) + 1L]] <- correlation_row(eicu, "eICU-CRD", "eICU_harmonized_M4", "gv_sd", "mean_glucose", "preserved aggregate all-source feature; optional mean-glucose benchmark")
rows[[length(rows) + 1L]] <- data.frame(
  database = "INSPIRE", analysis_context = "INSPIRE_final_v5", variable_x = "GV_SD", variable_y = NA_character_,
  N = 1353L, pearson_r = NA_real_, spearman_rho = NA_real_,
  measurement_structure = "final validated count/span fields retained for context table",
  analysis_status = "NOT_RUN",
  interpretation_boundary = "No new INSPIRE sampling-process analysis added in Phase 2A; no unsupported extrapolation",
  stringsAsFactors = FALSE
)
sampling <- do.call(rbind, rows)
sampling$lineage_version <- "JAHA_v5_final_lineage"
phase2a_write(sampling, phase2a_output_path(cfg, "machine_readable", "phase2a_sampling_process.csv"))
cat("PHASE2A_SAMPLING_DONE\n")
