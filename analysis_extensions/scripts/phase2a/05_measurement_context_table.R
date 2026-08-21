#!/usr/bin/env Rscript

phase2a_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2a_file, mustWork = TRUE)), "phase2a_common.R"))
suppressMessages({ library(jsonlite) })
cfg <- phase2a_load_config()
phase2a_open_log(cfg, "phase2a_05_measurement_context_table.log")
on.exit(phase2a_close_log(), add = TRUE)

if (!phase2a_phase1_gate(cfg)) stop("PHASE1_6_GATE_FAIL")
context <- phase2a_read(phase2a_output_path(cfg, "machine_readable", "final_measurement_context.csv"))
spec <- phase2a_read(phase2a_output_path(cfg, "machine_readable", "phase2a_specification_shift.csv"))
mice <- phase2a_read(phase2a_cfg_path(cfg, "mimic_mice_results"))
eicu_locked <- phase2a_read(phase2a_cfg_path(cfg, "eicu_locked_results"))
inspire_locked <- phase2a_read(phase2a_cfg_path(cfg, "inspire_primary_results"))

summary_measurement <- function(dd, count_col = "measurement_count", span_col = "measurement_span_hours") {
  count <- phase2a_num(dd[[count_col]])
  span <- phase2a_num(dd[[span_col]])
  qs <- function(x) {
    x <- x[is.finite(x)]
    q <- quantile(x, c(.25, .50, .75), names = FALSE, type = 7)
    c(q1 = q[1], median = q[2], q3 = q[3])
  }
  qc <- qs(count); qs_span <- qs(span)
  list(
    count = sprintf("median %.2f (Q1 %.2f, Q3 %.2f)", qc["median"], qc["q1"], qc["q3"]),
    span = sprintf("median %.2f h (Q1 %.2f, Q3 %.2f h)", qs_span["median"], qs_span["q1"], qs_span["q3"])
  )
}

make_context <- function(database, clinical_setting, time_anchor, exposure_window, landmark, source_structure, sampling_structure, gv_definition, dd, outcome, effect_type, core_model, effect, lower, upper, boundary, secondary_note, event_col = "event_30d") {
  events <- sum(phase2a_num(dd[[event_col]]) == 1, na.rm = TRUE)
  ms <- summary_measurement(dd)
  data.frame(
    database = database,
    clinical_setting = clinical_setting,
    time_anchor = time_anchor,
    exposure_window = exposure_window,
    landmark = landmark,
    glucose_source_structure = source_structure,
    sampling_structure = sampling_structure,
    GV_definition = gv_definition,
    measurement_count_summary = ms$count,
    measurement_span_summary = ms$span,
    outcome_definition = outcome,
    effect_type = effect_type,
    core_adjustment_model = core_model,
    locked_effect = effect,
    locked_lower95 = lower,
    locked_upper95 = upper,
    N = nrow(dd),
    events = events,
    secondary_context_note = secondary_note,
    interpretive_boundary = boundary,
    lineage_version = "JAHA_v5_final_lineage",
    stringsAsFactors = FALSE
  )
}

mimic <- context[context$database == "MIMIC-IV" & context$analysis_context == "MIMIC_primary_24h_day1_landmark", , drop = FALSE]
eicu <- context[context$database == "eICU-CRD" & context$analysis_context == "eICU_harmonized_M4", , drop = FALSE]
inspire <- context[context$database == "INSPIRE" & context$analysis_context == "INSPIRE_24h_landmark", , drop = FALSE]
if (nrow(mimic) != 10561L || nrow(eicu) != 7115L || nrow(inspire) != 1353L) stop("MEASUREMENT_CONTEXT_DENOMINATOR_FAIL")
mb <- mice[mice$model_id == "MICE_B_30d", , drop = FALSE]
em4 <- eicu_locked[eicu_locked$model_id == "HARM_eICU-CRD_M4", , drop = FALSE]
ii2 <- inspire_locked[inspire_locked$model_id == "ADMINV5_I2_30d", , drop = FALSE]
if (any(vapply(list(mb, em4, ii2), nrow, integer(1)) != 1L)) stop("MEASUREMENT_CONTEXT_LOCKED_RESULT_FAIL")

rows <- list(
  make_context(
    "MIMIC-IV", "adult cardiac-surgery routine care", "surgery t0 to postoperative day-1 landmark", "final 24-hour source-priority glucose window", "day-1 landmark",
    "central laboratory > blood gas > POCT > ICU-charted source-priority series", "patient-level routine-care count and time span", "sample SD of retained glucose values",
    mimic, "30-day all-cause mortality after day-1 landmark", "HR", "final Model B", mb$HR_per10, mb$lo_per10, mb$hi_per10,
    "GV is a source- and sampling-dependent routine-care representation; do not treat as a gold-standard biological exposure.", "Model C is retained as a locked specification shift; 365-day events are available in the final MIMIC target."),
  make_context(
    "INSPIRE", "cardiac-surgery routine care", "operation end (opend)", "0-24-hour post-opend glucose window", "24-hour opend landmark",
    "final INSPIRE retained glucose features", "patient-level count and span retained; linked mortality coverage attested through day 30", "sample SD of retained glucose values",
    inspire, "30-day all-cause mortality after 24-hour landmark", "HR", "final I2: GV + flexible mean glucose + age", ii2$HR_per10, ii2$lo, ii2$hi,
    "Day-30 coverage is source-provenance attested; no patient-level coverage-end field is present; 365-day endpoint is not reintroduced.", "The final 48-hour landmark context is separately validated at N=1511/events=31; no pooling with the 24-hour frame."),
  make_context(
    "eICU-CRD", "critical-care cardiac-surgery routine care", "preserved aggregate 24-hour landmark frame", "preserved aggregate 24-hour glucose feature window", "final aggregate in_landmark mask",
    "preserved aggregate all-source glucose features", "aggregate count and span only; no event-level glucose file preserved", "aggregate glucose_sd_24h (sample-SD field in locked aggregate input)",
    eicu, "post-landmark hospital mortality", "RR", "final M4: GV + age + sex + procedure + diabetes + creatinine + APACHE + flexible mean glucose + log count + span", em4$RR_per10, em4$lo, em4$hi,
    "eICU effect is an RR, not an HR; aggregate feature context cannot support event-level reconstruction or direct interchangeability with MIMIC/INSPIRE.", "M3 is the locked before-specification context; M4 adds count and span.", event_col = "event_hospital_post_landmark")
)
table <- do.call(rbind, rows)
phase2a_write(table, phase2a_output_path(cfg, "machine_readable", "phase2a_measurement_context_table.csv"))
cat("PHASE2A_CONTEXT_TABLE_DONE\n")
