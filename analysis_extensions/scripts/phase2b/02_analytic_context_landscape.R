#!/usr/bin/env Rscript

phase2b_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2b_file, mustWork = TRUE)), "phase2b_common.R"))
cfg <- phase2b_load_config()
phase2a_open_log(cfg, "phase2b_02_analytic_context_landscape.log")
on.exit(phase2a_close_log(), add = TRUE)

phase2b_require_gates(cfg)
mice <- phase2b_read(phase2a_cfg_path(cfg, "mimic_mice_results"))
source_results <- phase2b_read(phase2a_cfg_path(cfg, "mimic_source_results"))
inspire_30 <- phase2b_read(phase2a_cfg_path(cfg, "inspire_primary_results"))
inspire_48 <- phase2b_read(phase2a_cfg_path(cfg, "inspire_48h_results"))
eicu <- phase2b_read(phase2a_cfg_path(cfg, "eicu_locked_results"))

mice_path <- phase2a_cfg_path(cfg, "mimic_mice_results")
source_path <- phase2a_cfg_path(cfg, "mimic_source_results")
inspire_30_path <- phase2a_cfg_path(cfg, "inspire_primary_results")
inspire_48_path <- phase2a_cfg_path(cfg, "inspire_48h_results")
eicu_path <- phase2a_cfg_path(cfg, "eicu_locked_results")

rows <- list()
provenance <- list()

add_locked <- function(context_id, database, context_family, cohort_definition, time_anchor, exposure_window,
                       source_definition, gv_definition, adjustment_model, outcome, effect_type,
                       result, result_id, origin_path, estimate_col, lower_col, upper_col, p_col,
                       interpretive_note, figure_panel, figure_order, same_patient_pair = FALSE,
                       fit_status = "LOCKED_RESULT_VERIFIED", verified_by = NA_character_,
                       verified_by_selector = NA_character_) {
  if (nrow(result) != 1L) stop("LANDSCAPE_RESULT_ROW_MISSING: ", result_id)
  estimate <- phase2b_scalar(result[[estimate_col]])
  lower <- phase2b_scalar(result[[lower_col]])
  upper <- phase2b_scalar(result[[upper_col]])
  p_value <- phase2b_scalar(result[[p_col]])
  N <- phase2b_scalar(result$N)
  events <- phase2b_scalar(result$events)
  if (any(!is.finite(c(estimate, lower, upper, p_value, N, events)))) stop("LANDSCAPE_RESULT_NONFINITE: ", result_id)
  origin <- paste0(origin_path, "#", result_id)
  rows[[length(rows) + 1L]] <<- data.frame(
    context_id = context_id,
    database = database,
    context_family = context_family,
    cohort_definition = cohort_definition,
    time_anchor = time_anchor,
    exposure_window = exposure_window,
    source_definition = source_definition,
    gv_definition = gv_definition,
    adjustment_model = adjustment_model,
    outcome = outcome,
    effect_type = effect_type,
    N = as.integer(N),
    events = as.integer(events),
    estimate = estimate,
    lower95 = lower,
    upper95 = upper,
    p_value = p_value,
    estimate_origin = origin,
    new_or_locked = "LOCKED",
    interpretive_note = interpretive_note,
    fit_status = fit_status,
    figure_panel = figure_panel,
    figure_order = figure_order,
    same_patient_pair = same_patient_pair,
    lineage_version = "JAHA_v5_final_lineage",
    stringsAsFactors = FALSE
  )
  verified_path <- if (is.na(verified_by)) NA_character_ else verified_by
  provenance[[length(provenance) + 1L]] <<- data.frame(
    context_id = context_id,
    estimate_origin = origin,
    origin_type = "LOCKED_RESULT",
    origin_path = origin_path,
    origin_sha256 = phase2b_sha256(origin_path),
    row_selector = paste0("model_id=", result_id),
    verified_by = verified_path,
    verified_by_sha256 = if (is.na(verified_path)) NA_character_ else phase2b_sha256(verified_path),
    verified_by_selector = verified_by_selector,
    verification_status = fit_status,
    lineage_version = "JAHA_v5_final_lineage",
    stringsAsFactors = FALSE
  )
}

gv_definition_mimic <- "sample SD of retained glucose values; per 10 mg/dL GV"
gv_definition_inspire <- "sample SD of retained glucose values; per 10 mg/dL GV"
gv_definition_eicu <- "aggregate glucose_sd_24h; per 10 mg/dL GV"

for (i in seq_along(c("A", "B", "C"))) {
  id <- paste0("MIMIC_ADJUSTMENT_", c("A", "B", "C")[i], "_30D")
  model_id <- paste0("MICE_", c("A", "B", "C")[i], "_30d")
  rr <- mice[mice$model_id == model_id, , drop = FALSE]
  add_locked(
    id, "MIMIC-IV", "MIMIC adjustment context", "final MIMIC primary target",
    "surgery t0 to postoperative day-1 landmark", "final 24-hour source-priority glucose window",
    "central laboratory > blood gas > POCT > ICU-charted source-priority series",
    gv_definition_mimic, paste0("Model ", c("A", "B", "C")[i]),
    "30-day all-cause mortality after day-1 landmark", "HR", rr, model_id, mice_path,
    "HR_per10", "lo_per10", "hi_per10", "P_per10",
    "Locked MIMIC adjustment context; Model C is retained as a sampling-process specification shift from Model B.",
    "MIMIC_adjustment", i
  )
}

add_locked(
  "MIMIC_SOURCE_PRIORITY_30D", "MIMIC-IV", "MIMIC measurement-source context", "final MIMIC primary target",
  "surgery t0 to postoperative day-1 landmark", "final 24-hour source-priority glucose window",
  "central laboratory > blood gas > POCT > ICU-charted source-priority series",
  gv_definition_mimic, "Model B (primary priority series)",
  "30-day all-cause mortality after day-1 landmark", "HR", mice[mice$model_id == "MICE_B_30d", , drop = FALSE],
  "MICE_B_30d", mice_path, "HR_per10", "lo_per10", "hi_per10", "P_per10",
  "Primary priority-series estimate; included in the source-context family to anchor source comparisons.",
  "MIMIC_source", 1
)

source_specs <- list(
  list(id = "MIMIC_SOURCE_POCT_ONLY_30D", model_id = "SRC_poct_30d", label = "POCT-only series", order = 2, same = FALSE,
       note = "Source-specific final priority-derived series; not a same-patient paired estimate."),
  list(id = "MIMIC_SOURCE_CENTRAL_LAB_ONLY_30D", model_id = "SRC_centrallab_30d", label = "central-laboratory-only series", order = 3, same = FALSE,
       note = "Source-specific final priority-derived series; not a same-patient paired estimate."),
  list(id = "MIMIC_SOURCE_BLOOD_GAS_ONLY_30D", model_id = "SRC_bloodgas_30d", label = "blood-gas-only series", order = 4, same = FALSE,
       note = "Source-specific final priority-derived series; not a same-patient paired estimate."),
  list(id = "MIMIC_SOURCE_COMMON_30D", model_id = "SRC_common_30d", label = "common-source (POCT + central lab) series", order = 5, same = FALSE,
       note = "Common-source series; cohort differs from the same-patient comparison cohort."),
  list(id = "MIMIC_SOURCE_SAME_PATIENT_POCT_30D", model_id = "SRC_samepatient_poct_30d", label = "same-patient POCT-derived GV", order = 6, same = TRUE,
       note = "Identical 409 patients and 49 events as the same-patient laboratory row; source-definition coefficient difference is descriptive, not causal."),
  list(id = "MIMIC_SOURCE_SAME_PATIENT_LAB_30D", model_id = "SRC_samepatient_lab_30d", label = "same-patient central-laboratory-derived GV", order = 7, same = TRUE,
       note = "Identical 409 patients and 49 events as the same-patient POCT row; source-definition coefficient difference is descriptive, not causal.")
)
for (spec in source_specs) {
  rr <- source_results[source_results$model_id == spec$model_id, , drop = FALSE]
  verified_by <- if (spec$same) phase2b_output_path(cfg, "machine_readable", "phase2a_source_model_reproduction.csv") else NA_character_
  verified_selector <- if (spec$same) paste0("source=", if (grepl("poct", spec$model_id)) "POCT" else "central_laboratory", "; reproduction_status=PASS") else NA_character_
  add_locked(
    spec$id, "MIMIC-IV", "MIMIC measurement-source context", as.character(rr$cohort[1]),
    "surgery t0 to postoperative day-1 landmark", "final 24-hour source-specific glucose window",
    spec$label, gv_definition_mimic, "Model B (same as primary)",
    "30-day all-cause mortality after day-1 landmark", "HR", rr, spec$model_id, source_path,
    "HR_per10", "lo", "hi", "P", spec$note, "MIMIC_source", spec$order,
    same_patient_pair = spec$same,
    fit_status = if (spec$same) "LOCKED_RESULT_REPRODUCED_PHASE2A" else "LOCKED_RESULT_VERIFIED",
    verified_by = verified_by, verified_by_selector = verified_selector
  )
}

inspire_specs <- list(
  list(id = "INSPIRE_I1_24H_30D", model_id = "ADMINV5_I1_30d", label = "I1", anchor = "operation end (opend), 24-hour landmark", window = "0-24-hour post-opend glucose window", order = 1, note = "Final v5 24-hour timing/adjustment context."),
  list(id = "INSPIRE_I2_24H_30D", model_id = "ADMINV5_I2_30d", label = "I2", anchor = "operation end (opend), 24-hour landmark", window = "0-24-hour post-opend glucose window", order = 2, note = "Final v5 designated 24-hour model."),
  list(id = "INSPIRE_I3_24H_30D", model_id = "ADMINV5_I3_30d", label = "I3", anchor = "operation end (opend), 24-hour landmark", window = "0-24-hour post-opend glucose window", order = 3, note = "Final v5 I3 result verified; retained as a locked specification context."),
  list(id = "INSPIRE_I2_48H_CORRECTED_30D", model_id = "ADMINV5_I2_48h_landmark_30d", label = "I2 corrected 48h", anchor = "operation end (opend), 48-hour landmark", window = "0-48-hour post-opend glucose window", order = 4, note = "Corrected final v5 48-hour landmark; risk begins at the 48-hour landmark; no withdrawn 365-day result included.")
)
for (spec in inspire_specs) {
  if (grepl("48h", spec$model_id)) {
    rr <- inspire_48[inspire_48$model_id == spec$model_id, , drop = FALSE]
    origin <- inspire_48_path
  } else {
    rr <- inspire_30[inspire_30$model_id == spec$model_id, , drop = FALSE]
    origin <- inspire_30_path
  }
  add_locked(
    spec$id, "INSPIRE", "INSPIRE timing/adjustment context", "final INSPIRE validated frame",
    spec$anchor, spec$window, "final INSPIRE retained glucose features", gv_definition_inspire,
    spec$label, "30-day all-cause mortality after INSPIRE landmark", "HR", rr, spec$model_id, origin,
    "HR_per10", "lo", "hi", "P", spec$note, "INSPIRE_timing_adjustment", spec$order
  )
}

eicu_specs <- lapply(c("M1", "M2", "M3", "M4"), function(m) list(model_id = paste0("HARM_eICU-CRD_", m), label = m, order = match(m, c("M1", "M2", "M3", "M4"))))
for (spec in eicu_specs) {
  rr <- eicu[eicu$model_id == spec$model_id, , drop = FALSE]
  cohort_note <- if (spec$label == "M1") {
    "M1 uses the broader locked aggregate frame (N/events differ from final M2-M4 frame)."
  } else {
    "M2-M4 use the final aggregate frame N=7,115/events=130; M4 adds count and span to the locked aggregate context."
  }
  add_locked(
    paste0("EICU_", spec$label, "_24H_AGGREGATE"), "eICU-CRD", "eICU adjustment context",
    paste0("locked aggregate frame for ", spec$label), "preserved aggregate 24-hour landmark frame",
    "preserved aggregate 24-hour glucose feature window", "preserved aggregate all-source glucose features",
    gv_definition_eicu, paste0("modified Poisson ", spec$label), "post-landmark hospital mortality", "RR",
    rr, spec$model_id, eicu_path, "RR_per10", "lo", "hi", "P",
    paste0(cohort_note, " Event-level glucose reconstruction is unavailable; RR is not interchangeable with MIMIC/INSPIRE HR."),
    "eICU_adjustment", spec$order
  )
}

landscape <- do.call(rbind, rows)
landscape <- landscape[order(match(landscape$figure_panel, c("MIMIC_adjustment", "MIMIC_source", "INSPIRE_timing_adjustment", "eICU_adjustment")), landscape$figure_order), , drop = FALSE]
if (nrow(landscape) != 18L) stop("PHASE2B_LANDSCAPE_ROW_COUNT_FAIL")
if (any(!is.finite(landscape$estimate)) || any(!is.finite(landscape$lower95)) || any(!is.finite(landscape$upper95)) || any(!is.finite(landscape$p_value))) stop("PHASE2B_LANDSCAPE_NONFINITE")
phase2b_write(landscape, phase2b_output_path(cfg, "machine_readable", "phase2b_analytic_context_landscape.csv"))

provenance_table <- do.call(rbind, provenance)
provenance_table <- provenance_table[match(landscape$context_id, provenance_table$context_id), , drop = FALSE]
phase2b_write(provenance_table, phase2b_output_path(cfg, "machine_readable", "phase2b_result_provenance.csv"))

cat("PHASE2B_LANDSCAPE_DONE\n")
cat("LANDSCAPE_ROWS", nrow(landscape), "\n")
