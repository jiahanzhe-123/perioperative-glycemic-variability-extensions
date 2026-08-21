#!/usr/bin/env Rscript

# Phase 2C evidence-lock assembly.
# This script only reads locked Phase 1.6 / Phase 2A / Phase 2B artifacts and
# writes provenance-bearing registries and planning documents. It contains no
# model-fitting, resampling, correlation, agreement, or specification search.

phase2c_arg <- commandArgs(trailingOnly = FALSE)
phase2c_file_arg <- grep("^--file=", phase2c_arg, value = TRUE)
if (length(phase2c_file_arg) == 0L) stop("Phase 2C scripts must be run with Rscript")
phase2c_script_path <- normalizePath(sub("^--file=", "", phase2c_file_arg[1]), mustWork = TRUE)
phase2c_workspace <- normalizePath(file.path(dirname(phase2c_script_path), "../.."), mustWork = TRUE)

source(file.path(phase2c_workspace, "scripts", "phase2b", "phase2b_common.R"))
suppressMessages(library(jsonlite))

cfg <- phase2b_load_config()
phase2b_require_gates(cfg)

out <- function(...) phase2b_output_path(cfg, ...)
cfg_path <- function(key) phase2a_cfg_path(cfg, key)
script_path <- function(...) normalizePath(file.path(phase2c_workspace, "scripts", ...), mustWork = FALSE)
read_csv <- function(path) phase2b_read(path)
num <- function(x) suppressWarnings(as.numeric(as.character(x[1])))
text_num <- function(x) {
  z <- num(x)
  if (!is.finite(z)) return(NA_character_)
  format(z, digits = 17, scientific = TRUE, trim = TRUE)
}
text_value <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_character_)
  paste(vapply(x, text_num, character(1)), collapse = ";")
}
pick_one <- function(data, column, value, label = column) {
  hit <- data[as.character(data[[column]]) == as.character(value), , drop = FALSE]
  if (nrow(hit) != 1L) stop("Expected exactly one row for ", label, "=", value, "; found ", nrow(hit))
  hit[1, , drop = FALSE]
}
sha256 <- function(path) phase2b_sha256(path)

phase2a_script_agreement <- script_path("phase2a", "01_same_patient_agreement.R")
phase2a_script_bootstrap <- script_path("phase2a", "02_same_patient_paired_bootstrap.R")
phase2a_script_sampling <- script_path("phase2a", "03_sampling_process_dependence.R")
phase2a_script_shift <- script_path("phase2a", "04_model_specification_shift.R")
phase2b_script_landscape <- script_path("phase2b", "02_analytic_context_landscape.R")
phase2b_script_figures <- script_path("phase2b", "03_candidate_figures.R")

agreement_path <- out("machine_readable", "phase2a_same_patient_agreement.csv")
bootstrap_path <- out("machine_readable", "phase2a_same_patient_bootstrap_summary.csv")
sampling_path <- out("machine_readable", "phase2a_sampling_process.csv")
shift_path <- out("machine_readable", "phase2a_specification_shift.csv")
source_model_path <- out("machine_readable", "phase2a_source_model_reproduction.csv")
log_path <- out("machine_readable", "phase2b_log_scale_agreement.csv")
imbalance_path <- out("machine_readable", "phase2b_source_sampling_imbalance.csv")
landscape_path <- out("machine_readable", "phase2b_analytic_context_landscape.csv")
landscape_prov_path <- out("machine_readable", "phase2b_result_provenance.csv")
summary_path <- out("qc", "final_reconstruction_summary.json")

agreement <- read_csv(agreement_path)
bootstrap <- read_csv(bootstrap_path)
sampling <- read_csv(sampling_path)
shift <- read_csv(shift_path)
source_model <- read_csv(source_model_path)
log_scale <- read_csv(log_path)
imbalance <- read_csv(imbalance_path)
landscape <- read_csv(landscape_path)
mice <- read_csv(cfg_path("mimic_mice_results"))
primary <- read_csv(cfg_path("mimic_primary_results"))
source_results <- read_csv(cfg_path("mimic_source_results"))
inspire <- read_csv(cfg_path("inspire_primary_results"))
inspire48 <- read_csv(cfg_path("inspire_48h_results"))
eicu <- read_csv(cfg_path("eicu_locked_results"))
absolute_risk <- read_csv(file.path(cfg_path("analysis_record_root"), "results", "12_absolute_risk_results.csv"))
ph_diag <- read_csv(file.path(cfg_path("analysis_record_root"), "results", "07_ph_diagnostics.csv"))
reconstruction <- fromJSON(summary_path)

registry_cols <- c(
  "evidence_id", "domain", "analysis", "population", "N", "events", "estimand",
  "effect_type", "estimate", "lower95", "upper95", "p_value", "secondary_statistic",
  "secondary_value", "source_status", "source_file", "source_row_or_key", "script_origin",
  "new_or_original", "phase_origin", "candidate_manuscript_location", "interpretive_role", "notes"
)

registry_row <- function(evidence_id, domain, analysis, population, N = NA_real_, events = NA_real_,
                         estimand = NA_character_, effect_type = NA_character_, estimate = NA_real_,
                         lower95 = NA_real_, upper95 = NA_real_, p_value = NA_real_,
                         secondary_statistic = NA_character_, secondary_value = NA_character_,
                         source_status = NA_character_, source_file = NA_character_,
                         source_row_or_key = NA_character_, script_origin = NA_character_,
                         new_or_original = NA_character_, phase_origin = NA_character_,
                         candidate_manuscript_location = NA_character_, interpretive_role = NA_character_,
                         notes = NA_character_) {
  data.frame(
    evidence_id = as.character(evidence_id), domain = as.character(domain),
    analysis = as.character(analysis), population = as.character(population),
    N = as.numeric(N), events = as.numeric(events), estimand = as.character(estimand),
    effect_type = as.character(effect_type), estimate = as.numeric(estimate),
    lower95 = as.numeric(lower95), upper95 = as.numeric(upper95), p_value = as.numeric(p_value),
    secondary_statistic = as.character(secondary_statistic), secondary_value = as.character(secondary_value),
    source_status = as.character(source_status), source_file = as.character(source_file),
    source_row_or_key = as.character(source_row_or_key), script_origin = as.character(script_origin),
    new_or_original = as.character(new_or_original), phase_origin = as.character(phase_origin),
    candidate_manuscript_location = as.character(candidate_manuscript_location),
    interpretive_role = as.character(interpretive_role), notes = as.character(notes),
    stringsAsFactors = FALSE
  )[registry_cols]
}

registry <- registry_row("__placeholder__", "assembly", "placeholder", "none")
registry <- registry[0, , drop = FALSE]
add <- function(...) registry <<- rbind(registry, registry_row(...))

model_notes <- paste(
  "Per 10 mg/dL SD-GV coefficient; source row is the locked final-lineage result.",
  "Observational association; no causal interpretation.",
  sep = " "
)

add_mice <- function(row, id, horizon, model_label, location, role, notes = model_notes) {
  horizon_label <- if (horizon == "30d") "30-day" else "365-day"
  add(
    id, "MIMIC adjustment context", paste0("MIMIC Model ", model_label, " ", horizon_label),
    "final MIMIC primary target", row$N, row$events,
    paste0("MIMIC day-1 landmark mortality by ", horizon_label), "HR per 10 mg/dL GV",
    row$HR_per10, row$lo_per10, row$hi_per10, row$P_per10,
    "n_imputations_fit", text_num(row$n_imputations_fit), "LOCKED_RESULT_VERIFIED",
    cfg_path("mimic_mice_results"), paste0("model_id=", row$model_id), cfg_path("public_primary_model_script"),
    "original_final_lineage", "Phase 1.6", location, role, notes
  )
}

mice_a30 <- pick_one(mice, "model_id", "MICE_A_30d")
mice_b30 <- pick_one(mice, "model_id", "MICE_B_30d")
mice_c30 <- pick_one(mice, "model_id", "MICE_C_30d")
mice_a365 <- pick_one(mice, "model_id", "MICE_A_365d")
mice_b365 <- pick_one(mice, "model_id", "MICE_B_365d")
mice_c365 <- pick_one(mice, "model_id", "MICE_C_365d")
add_mice(mice_a30, "MIMIC_MODEL_A_30D", "30d", "A", "Results; MIMIC adjustment context", "Main adjustment context")
add_mice(mice_b30, "MIMIC_MODEL_B_30D", "30d", "B", "Results; primary MIMIC estimate", "Primary locked MIMIC estimate")
add_mice(mice_c30, "MIMIC_MODEL_C_30D", "30d", "C", "Results; MIMIC adjustment context", "Measurement-process specification shift")
add_mice(mice_a365, "MIMIC_MODEL_A_365D", "365d", "A", "Supplement; MIMIC 365-day audit", "Supplementary time-horizon context")
add_mice(mice_b365, "MIMIC_MODEL_B_365D", "365d", "B", "Supplement; MIMIC 365-day audit", "Descriptive average under documented PH violation",
         "365-day Model B is not a primary effect: the final manuscript provenance documents non-proportional hazards; retain only as a prespecified descriptive average.")
add_mice(mice_c365, "MIMIC_MODEL_C_365D", "365d", "C", "Supplement; MIMIC 365-day audit", "Supplementary time-horizon context",
         "365-day Model C is retained only as a descriptive specification sensitivity under the documented PH limitation.")

add(
  "MIMIC_PRIMARY_COHORT_30D", "MIMIC cohort context", "Final MIMIC primary cohort", "final MIMIC primary target",
  reconstruction$mimic$target_n, reconstruction$mimic$events_30d, "day-1 landmark cohort size and 30-day events",
  "cohort N/events", reconstruction$mimic$target_n, NA_real_, NA_real_, NA_real_, "events", text_num(reconstruction$mimic$events_30d),
  "LOCKED_CONTEXT_SUMMARY", summary_path, "mimic.target_n; mimic.events_30d", script_path("03_final_lineage_qc.py"),
  "original_final_lineage", "Phase 1.6", "Methods; Results", "Cohort denominator", "N and 30-day event count are fixed lineage facts."
)
add(
  "MIMIC_PRIMARY_COHORT_365D", "MIMIC cohort context", "Final MIMIC primary cohort", "final MIMIC primary target",
  reconstruction$mimic$target_n, reconstruction$mimic$events_365d, "day-1 landmark cohort size and 365-day events",
  "cohort N/events", reconstruction$mimic$target_n, NA_real_, NA_real_, NA_real_, "events", text_num(reconstruction$mimic$events_365d),
  "LOCKED_CONTEXT_SUMMARY", summary_path, "mimic.target_n; mimic.events_365d", script_path("03_final_lineage_qc.py"),
  "original_final_lineage", "Phase 1.6", "Supplement; MIMIC 365-day audit", "Cohort denominator", "365-day event count is retained for the time-horizon audit only."
)

add_context_summary <- function(id, analysis, estimate, secondary_statistic, secondary_value, location, role, notes, source_key) {
  add(
    id, "MIMIC measurement context", analysis, "final MIMIC primary target", reconstruction$mimic$target_n,
    reconstruction$mimic$events_30d, "descriptive final-lineage measurement context", "descriptive context statistic",
    estimate, NA_real_, NA_real_, NA_real_, secondary_statistic, secondary_value,
    "LOCKED_CONTEXT_SUMMARY", summary_path, source_key, script_path("03_final_lineage_qc.py"),
    "original_final_lineage", "Phase 1.6", location, role, notes
  )
}
add_context_summary("MIMIC_MEAN_GV", "Final source-priority GV mean", reconstruction$mimic$gv_mean,
                   "units", "mg/dL; sample SD-derived GV", "Methods; Results", "Primary GV scale context",
                   "Mean GV is a descriptive feature summary, not an external reference value.", "mimic.gv_mean")
add_context_summary("MIMIC_MEAN_GLUCOSE", "Final source-priority mean glucose", reconstruction$mimic$mean_glucose_mean,
                   "units", "mg/dL", "Methods; Results", "Concurrent glucose context",
                   "Mean glucose is retained to define the conditioning context for Model B.", "mimic.mean_glucose_mean")
add_context_summary("MIMIC_MEASUREMENT_COUNT_MEDIAN", "Routine-care measurement count", reconstruction$mimic$count_median,
                   "Q1;Q3", paste(text_num(reconstruction$mimic$count_q1), text_num(reconstruction$mimic$count_q3), sep = ";"),
                   "Methods; Results", "Sampling opportunity context",
                   "Median and IQR of retained source-priority measurements; descriptive only.", "mimic.count_median;mimic.count_q1;mimic.count_q3")
add_context_summary("MIMIC_MEASUREMENT_SPAN_MEDIAN", "Routine-care measurement span", reconstruction$mimic$span_median,
                   "Q1;Q3", paste(text_num(reconstruction$mimic$span_q1), text_num(reconstruction$mimic$span_q3), sep = ";"),
                   "Methods; Results", "Sampling opportunity context",
                   "Median and IQR of retained source-priority span in hours; descriptive only.", "mimic.span_median;mimic.span_q1;mimic.span_q3")

add(
  "ORIG_MIMIC_ABSOLUTE_RISK_30D", "MIMIC outcome presentation", "Model-standardized GV quartile risk contrast", "complete-case standardization frame",
  absolute_risk$N[absolute_risk$model_id == "ABSRISK_30d"], absolute_risk$events[absolute_risk$model_id == "ABSRISK_30d"],
  "30-day model-standardized Q75 versus Q25 risk difference", "risk difference",
  absolute_risk$rd[absolute_risk$model_id == "ABSRISK_30d"], absolute_risk$rd_lo[absolute_risk$model_id == "ABSRISK_30d"],
  absolute_risk$rd_hi[absolute_risk$model_id == "ABSRISK_30d"], NA_real_, "risk_q25;risk_q75",
  text_value(c(absolute_risk$risk_q25[absolute_risk$model_id == "ABSRISK_30d"], absolute_risk$risk_q75[absolute_risk$model_id == "ABSRISK_30d"])),
  "LOCKED_RESULT_VERIFIED", file.path(cfg_path("analysis_record_root"), "results", "12_absolute_risk_results.csv"),
  "model_id=ABSRISK_30d", cfg_path("public_primary_model_script"), "original_final_lineage", "Phase 1.6",
  "Supplement; original Figure 3 evidence", "Valid secondary absolute-risk presentation",
  "Retain as a supplementary descriptive contrast; it is not the revised manuscript's central estimand."
)

ph365_gv <- ph_diag[ph_diag$model_id == "ModelB_365d" & ph_diag$term == "gv10", , drop = FALSE]
if (nrow(ph365_gv) != 1L) stop("Expected one ModelB_365d gv10 PH row")
ph365_global <- ph_diag[ph_diag$model_id == "ModelB_365d" & ph_diag$term == "GLOBAL", , drop = FALSE]
if (nrow(ph365_global) != 1L) stop("Expected one ModelB_365d GLOBAL PH row")
add(
  "ORIG_MIMIC_365D_PH_DIAGNOSTIC", "MIMIC time-context audit", "365-day Model B proportional-hazards diagnostic", "complete-case diagnostic frame",
  9751, 633, "diagnostic evidence for 365-day average HR", "PH diagnostic p-value", ph365_gv$p, NA_real_, NA_real_, NA_real_,
  "global_p", text_num(ph365_global$p), "LOCKED_RESULT_VERIFIED", file.path(cfg_path("analysis_record_root"), "results", "07_ph_diagnostics.csv"),
  "model_id=ModelB_365d; term=gv10", cfg_path("public_primary_model_script"), "original_final_lineage", "Phase 1.6",
  "Supplement; MIMIC 365-day audit", "Boundary condition", "The GV and global PH diagnostics support demoting the 365-day average HR from the revised main story."
)

add_source_model <- function(model_id, evidence_id, source_label, location, role, notes) {
  row <- pick_one(source_results, "model_id", model_id)
  same_pair <- grepl("samepatient", model_id)
  add(
    evidence_id, "MIMIC measurement-source context", paste0("Final source-defined GV: ", source_label),
    row$cohort, row$N, row$events, "30-day mortality after day-1 landmark", "HR per 10 mg/dL GV",
    row$HR_per10, row$lo, row$hi, row$P, "gv_sd_within", text_num(row$gv_sd_within),
    if (same_pair) "LOCKED_RESULT_REPRODUCED_PHASE2A" else "LOCKED_RESULT_VERIFIED", cfg_path("mimic_source_results"),
    paste0("model_id=", model_id), cfg_path("public_source_model_script"),
    if (same_pair) "original_final_lineage_verified_phase2a" else "original_final_lineage", if (same_pair) "Phase 1.6 + Phase 2A" else "Phase 1.6",
    location, role, notes
  )
}
add_source_model("SRC_poct_30d", "MIMIC_SOURCE_POCT_ONLY_30D", "POCT-only", "Results; Figure 2", "Source-specific context",
                "Source-specific cohort; not a same-patient paired comparison.")
add_source_model("SRC_centrallab_30d", "MIMIC_SOURCE_CENTRAL_LAB_ONLY_30D", "central-laboratory-only", "Results; Figure 2", "Source-specific context",
                "Source-specific cohort; not a same-patient paired comparison.")
add_source_model("SRC_bloodgas_30d", "MIMIC_SOURCE_BLOOD_GAS_ONLY_30D", "blood-gas-only", "Results; Figure 2", "Source-specific context",
                "Source-specific cohort; not a same-patient paired comparison.")
add_source_model("SRC_common_30d", "MIMIC_SOURCE_COMMON_30D", "common-source", "Results; Figure 2", "Source-specific context",
                "Common-source cohort; not a same-patient paired comparison.")
add_source_model("SRC_samepatient_poct_30d", "MIMIC_SAMEPATIENT_POCT_HR", "same-patient POCT", "Results; Figure 2", "Identical-patient source comparison",
                 "Same N=409 and 49 events as the laboratory row; source is not a reference standard.")
add_source_model("SRC_samepatient_lab_30d", "MIMIC_SAMEPATIENT_LAB_HR", "same-patient laboratory", "Results; Figure 2", "Identical-patient source comparison",
                 "Same N=409 and 49 events as the POCT row; source is not a reference standard.")

add_inspire <- function(model_id, evidence_id, location, role, notes) {
  row <- pick_one(inspire, "model_id", model_id)
  add(
    evidence_id, "INSPIRE timing/adjustment context", paste0("INSPIRE ", row$model_id), row$frame, row$N, row$events,
    "30-day all-cause mortality after operation-end landmark", "HR per 10 mg/dL GV", row$HR_per10, row$lo, row$hi, row$P,
    "ph_global", text_num(row$ph_global), "LOCKED_RESULT_VERIFIED", cfg_path("inspire_primary_results"), paste0("model_id=", model_id),
    cfg_path("controlled_inspire_runner"), "original_final_lineage", "Phase 1.6", location, role, notes
  )
}
add_inspire("ADMINV5_I1_30d", "INSPIRE_I1_24H_30D", "Results; Figure 3", "Timing/adjustment context", "Final v5 24-hour operation-end landmark; event-limited.")
add_inspire("ADMINV5_I2_30d", "INSPIRE_I2_24H_30D", "Results; Figure 3", "Primary INSPIRE context", "Final v5 I2; day-30 coverage is source-provenance attested and no patient-level coverage-end field is present.")
add_inspire("ADMINV5_I3_30d", "INSPIRE_I3_24H_30D", "Results; Figure 3", "Measurement-process specification shift", "Final v5 I3; count/span are process-context terms, not causal confounders.")
row48 <- pick_one(inspire48, "model_id", "ADMINV5_I2_48h_landmark_30d")
add(
  "INSPIRE_CORRECTED_48H_30D", "INSPIRE timing context", "Corrected INSPIRE 48-hour landmark sensitivity", row48$frame, row48$N, row48$events,
  "30-day all-cause mortality after corrected 48-hour operation-end landmark", "HR per 10 mg/dL GV", row48$HR_per10, row48$lo, row48$hi, row48$P,
  "ph_global", text_num(row48$ph_global), "LOCKED_RESULT_VERIFIED", cfg_path("inspire_48h_results"), "model_id=ADMINV5_I2_48h_landmark_30d",
  cfg_path("controlled_inspire_runner"), "original_final_lineage", "Phase 1.6", "Supplement; Figure 3 timing sensitivity", "Corrected 48-hour landmark sensitivity; no exposure-risk overlap claim is made."
)

add_eicu <- function(model_id, evidence_id, location, role, notes) {
  row <- pick_one(eicu, "model_id", model_id)
  add(
    evidence_id, "eICU adjustment context", paste0("eICU ", row$model), "final locked aggregate eICU frame", row$N, row$events,
    "post-landmark hospital mortality", "RR per 10 mg/dL GV", row$RR_per10, row$lo, row$hi, row$P,
    "variance_estimator", row$variance, "LOCKED_RESULT_VERIFIED", cfg_path("eicu_locked_results"), paste0("model_id=", model_id),
    cfg_path("public_eicu_model_script"), "original_final_lineage", "Phase 1.6", location, role, notes
  )
}
add_eicu("HARM_eICU-CRD_M1", "EICU_M1", "Supplement; eICU context", "Adjustment context", "M1 uses the broader locked aggregate frame; effect type is RR.")
add_eicu("HARM_eICU-CRD_M2", "EICU_M2", "Supplement; eICU context", "Adjustment context", "M2 uses the final aggregate frame; effect type is RR.")
add_eicu("HARM_eICU-CRD_M3", "EICU_M3", "Results; Figure 3", "Primary eICU landscape context", "Final aggregate frame N=7,115/events=130; RR is not interchangeable with MIMIC/INSPIRE HR.")
add_eicu("HARM_eICU-CRD_M4", "EICU_M4", "Results; Figure 3", "Measurement-process specification shift", "M4 adds count and span to the preserved aggregate context; this is not a causal adjustment model.")

add_agreement <- function(row, cohort_tag, variable_label, id_prefix, location, role) {
  N <- row$N
  pop <- row$cohort_definition
  base <- paste0("same-patient source agreement: ", variable_label, " (", cohort_tag, ")")
  status <- "EXTENSION_RESULT_LOCKED"
  add(paste0(id_prefix, "_POCT_MEAN"), "same-patient source dependence", base, pop, N, NA_real_, "paired source-defined descriptive comparison", "POCT mean",
      row$poct_mean, NA_real_, NA_real_, NA_real_, "laboratory_mean", text_num(row$lab_mean), status, agreement_path,
      paste0("cohort_id=", row$cohort_id, "; variable=", row$variable), phase2a_script_agreement, "new_extension", "Phase 2A", location, role,
      "Means are source-defined routine-care summaries; neither source is a reference standard.")
  add(paste0(id_prefix, "_LAB_MEAN"), "same-patient source dependence", base, pop, N, NA_real_, "paired source-defined descriptive comparison", "laboratory mean",
      row$lab_mean, NA_real_, NA_real_, NA_real_, "POCT_mean", text_num(row$poct_mean), status, agreement_path,
      paste0("cohort_id=", row$cohort_id, "; variable=", row$variable), phase2a_script_agreement, "new_extension", "Phase 2A", location, role,
      "Means are source-defined routine-care summaries; neither source is a reference standard.")
  add(paste0(id_prefix, "_PEARSON"), "same-patient source dependence", base, pop, N, NA_real_, "paired source-defined descriptive comparison", "Pearson correlation",
      row$pearson_r, NA_real_, NA_real_, NA_real_, "Spearman_rho", text_num(row$spearman_rho), status, agreement_path,
      paste0("cohort_id=", row$cohort_id, "; variable=", row$variable), phase2a_script_agreement, "new_extension", "Phase 2A", location, role,
      "Correlation describes agreement/association only; it does not establish interchangeability.")
  add(paste0(id_prefix, "_SPEARMAN"), "same-patient source dependence", base, pop, N, NA_real_, NA_real_, "Spearman correlation",
      row$spearman_rho, NA_real_, NA_real_, NA_real_, "Pearson_r", text_num(row$pearson_r), status, agreement_path,
      paste0("cohort_id=", row$cohort_id, "; variable=", row$variable), phase2a_script_agreement, "new_extension", "Phase 2A", location, role,
      "Rank correlation describes agreement/association only; it does not establish interchangeability.")
  add(paste0(id_prefix, "_MEAN_DIFFERENCE"), "same-patient source dependence", base, pop, N, NA_real_, "POCT minus laboratory source-defined difference", "mean paired difference",
      row$mean_paired_difference, row$loa_lower95, row$loa_upper95, NA_real_, "sd_paired_difference", text_num(row$sd_paired_difference), status, agreement_path,
      paste0("cohort_id=", row$cohort_id, "; variable=", row$variable), phase2a_script_agreement, "new_extension", "Phase 2A", location, role,
      "lower95/upper95 are 95% Bland-Altman limits of agreement, not confidence limits for the mean.")
  add(paste0(id_prefix, "_PAIR_MEAN_CORRELATION"), "same-patient source dependence", base, pop, N, NA_real_, "paired source-defined difference versus pair mean", "Pearson correlation",
      row$difference_pair_mean_pearson_r, NA_real_, NA_real_, NA_real_, "Spearman_rho", text_num(row$difference_pair_mean_spearman_rho), status, agreement_path,
      paste0("cohort_id=", row$cohort_id, "; variable=", row$variable), phase2a_script_agreement, "new_extension", "Phase 2A", location, role,
      "Scale-dependent disagreement diagnostic; not a causal model.")
}

ag_gv_primary <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "GV_SD", , drop = FALSE]
ag_gv_sens <- agreement[agreement$cohort_id == "SENSITIVITY_ALL_PAIRED" & agreement$variable == "GV_SD", , drop = FALSE]
ag_mean_primary <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "MEAN_GLUCOSE", , drop = FALSE]
ag_mean_sens <- agreement[agreement$cohort_id == "SENSITIVITY_ALL_PAIRED" & agreement$variable == "MEAN_GLUCOSE", , drop = FALSE]
if (nrow(ag_gv_primary) != 1L || nrow(ag_gv_sens) != 1L || nrow(ag_mean_primary) != 1L || nrow(ag_mean_sens) != 1L) stop("Agreement rows are not uniquely identified")
add_agreement(ag_gv_primary, "primary N=452", "GV SD", "P2A_AGREEMENT_PRIMARY_GV", "Results; Figure 2", "Primary source-agreement evidence")
add_agreement(ag_gv_sens, "sensitivity N=453", "GV SD", "P2A_AGREEMENT_SENSITIVITY_GV", "Supplement; source-agreement sensitivity", "All source-paired candidate sensitivity")
add_agreement(ag_mean_primary, "primary N=452", "mean glucose", "P2A_AGREEMENT_PRIMARY_MEAN_GLUCOSE", "Supplement; source-agreement context", "Mean-glucose source context")
add_agreement(ag_mean_sens, "sensitivity N=453", "mean glucose", "P2A_AGREEMENT_SENSITIVITY_MEAN_GLUCOSE", "Supplement; source-agreement sensitivity", "Mean-glucose source context")

boot <- bootstrap[1, , drop = FALSE]
add(
  "P2A_BOOTSTRAP_OBSERVED_DELTA_BETA", "same-patient source dependence", "Paired source-definition coefficient-difference sensitivity analysis", "identical same-patient outcome-model cohort",
  boot$N, boot$events, "difference between paired source-specific log-HR coefficients", "delta-beta",
  boot$observed_delta_beta, boot$percentile_ci_delta_beta_lower, boot$percentile_ci_delta_beta_upper, NA_real_,
  "successful_replicates;failed_replicates", paste(text_num(boot$successful_replicates), text_num(boot$failed_replicates), sep = ";"),
  "EXTENSION_RESULT_LOCKED", bootstrap_path, "analysis=paired source-definition coefficient-difference sensitivity analysis", phase2a_script_bootstrap,
  "new_extension", "Phase 2A", "Results; Figure 2", "Paired source-definition coefficient-difference sensitivity analysis",
  "Paired bootstrap B=2,000 with 2,000 successful and 0 failed replicates; not an interaction test and not causal."
)
add(
  "P2A_BOOTSTRAP_COEFFICIENT_RATIO", "same-patient source dependence", "Paired source-definition coefficient-difference sensitivity analysis", "identical same-patient outcome-model cohort",
  boot$N, boot$events, "exponentiated paired log-HR coefficient difference", "paired coefficient ratio",
  boot$observed_coefficient_ratio, NA_real_, NA_real_, NA_real_, "POCT_HR;laboratory_HR", paste(text_num(boot$observed_hr_poct), text_num(boot$observed_hr_lab), sep = ";"),
  "EXTENSION_RESULT_LOCKED", bootstrap_path, "analysis=paired source-definition coefficient-difference sensitivity analysis", phase2a_script_bootstrap,
  "new_extension", "Phase 2A", "Results; Figure 2", "Paired source-definition coefficient-difference sensitivity analysis",
  "Exponentiated delta-beta is a descriptive paired coefficient ratio; not a causal interaction contrast."
)
add(
  "P2A_BOOTSTRAP_MEDIAN_DELTA_BETA", "same-patient source dependence", "Paired source-definition coefficient-difference sensitivity analysis", "identical same-patient outcome-model cohort",
  boot$N, boot$events, "bootstrap distribution of paired log-HR coefficient difference", "bootstrap median delta-beta",
  boot$bootstrap_median_delta_beta, NA_real_, NA_real_, NA_real_, "bootstrap_sd", text_num(boot$bootstrap_sd_delta_beta),
  "EXTENSION_RESULT_LOCKED", bootstrap_path, "analysis=paired source-definition coefficient-difference sensitivity analysis", phase2a_script_bootstrap,
  "new_extension", "Phase 2A", "Supplement; paired coefficient sensitivity", "Distributional bootstrap summary",
  "Bootstrap distribution summary; no new inferential interpretation."
)

for (i in seq_len(nrow(sampling))) {
  row <- sampling[i, , drop = FALSE]
  id <- paste0("P2A_SAMPLING_", gsub("[^A-Za-z0-9]+", "_", row$database), "_", gsub("[^A-Za-z0-9]+", "_", row$variable_y))
  completed <- identical(as.character(row$analysis_status), "COMPLETED")
  add(
    id, "sampling opportunity context", paste0("Descriptive GV versus ", row$variable_y, " relationship"), row$measurement_structure,
    row$N, NA_real_, "descriptive correlation within declared final context", "Pearson correlation",
    if (completed) row$pearson_r else NA_real_, NA_real_, NA_real_, NA_real_, "Spearman_rho", if (completed) text_num(row$spearman_rho) else NA_character_,
    if (completed) "EXTENSION_RESULT_LOCKED" else "NOT_RUN", sampling_path, paste0("database=", row$database, "; variable_y=", row$variable_y), phase2a_script_sampling,
    "new_extension", "Phase 2A", if (completed) "Supplement; sampling-process context" else "Supplement; declared non-estimable context",
    "Descriptive relationship only; not a confounder adjustment, causal model, or proof that sampling intensity explains GV association.",
    if (completed) "Completed Phase 2A correlation." else as.character(row$interpretation_boundary)
  )
}

for (i in seq_len(nrow(shift))) {
  row <- shift[i, , drop = FALSE]
  stem <- paste0("P2A_SHIFT_", gsub("[^A-Za-z0-9]+", "_", row$database), "_", gsub("[^A-Za-z0-9]+", "_", row$comparison))
  base_note <- "Locked specification movement; not a causal attenuation or confounding estimate."
  add(
    paste0(stem, "_BEFORE"), "analytic-context landscape", paste0(row$database, " ", row$comparison, " before"), "declared final context", row$N_before, row$events_before,
    paste0(row$database, " declared outcome context"), row$effect_type, row$effect_before, row$lower95_before, row$upper95_before, NA_real_,
    "effect_after", text_num(row$effect_after), "LOCKED_RESULT_VERIFIED", row$source_before, paste0("specification=", row$specification_before), phase2a_script_shift,
    "new_extension", "Phase 2A", "Supplement; analytic-context movement", "Before-context coefficient",
    paste(base_note, "After-context CI and movement are in the paired registry row.")
  )
  add(
    paste0(stem, "_AFTER"), "analytic-context landscape", paste0(row$database, " ", row$comparison, " after"), "declared final context", row$N_after, row$events_after,
    paste0(row$database, " declared outcome context"), row$effect_type, row$effect_after, row$lower95_after, row$upper95_after, NA_real_,
    "effect_before", text_num(row$effect_before), "LOCKED_RESULT_VERIFIED", row$source_after, paste0("specification=", row$specification_after), phase2a_script_shift,
    "new_extension", "Phase 2A", "Supplement; analytic-context movement", "After-context coefficient",
    paste(base_note, "Before-context CI and movement are in the paired registry row.")
  )
  add(
    paste0(stem, "_DELTA_LOG"), "analytic-context landscape", paste0(row$database, " ", row$comparison, " coefficient movement"), "declared final context", row$N_after, row$events_after,
    paste0(row$database, " declared outcome context"), paste0(row$effect_type, " log movement"), row$delta_log_effect, NA_real_, NA_real_, NA_real_,
    "null_crossing", as.character(row$null_crossing), "EXTENSION_RESULT_LOCKED", shift_path, paste0("database=", row$database, "; comparison=", row$comparison), phase2a_script_shift,
    "new_extension", "Phase 2A", "Supplement; analytic-context movement", "Coefficient movement diagnostic",
    paste(base_note, "The null-crossing flag is a descriptive interval property.")
  )
}

add_log_stat <- function(row, suffix, analysis, effect_type, estimate, location, role, note, secondary_statistic = NA_character_, secondary_value = NA_character_) {
  add(
    paste0("P2B_LOG_", ifelse(row$cohort_id == "PRIMARY_FINAL_TARGET", "PRIMARY", "SENSITIVITY"), "_", suffix),
    "same-patient source dependence", analysis, row$cohort_definition, row$N, NA_real_, "multiplicative-scale source-defined agreement", effect_type,
    estimate, NA_real_, NA_real_, NA_real_, secondary_statistic, secondary_value, "EXTENSION_RESULT_LOCKED", log_path,
    paste0("cohort_id=", row$cohort_id), script_path("phase2b", "01_source_dependence_closure.R"), "new_extension", "Phase 2B", location, role, note
  )
}
for (i in seq_len(nrow(log_scale))) {
  row <- log_scale[i, , drop = FALSE]
  cohort_short <- if (row$cohort_id == "PRIMARY_FINAL_TARGET") "primary" else "sensitivity"
  location <- if (cohort_short == "primary") "Results; Figure 2" else "Supplement; multiplicative-scale sensitivity"
  role <- if (cohort_short == "primary") "Multiplicative-scale sensitivity" else "Sensitivity agreement"
  note <- "Log-scale sensitivity does not replace the original-scale Bland-Altman analysis; neither source is a reference standard."
  add_log_stat(row, "VALID_N", "Multiplicative-scale agreement valid-pair count", "valid positive pairs", row$N, location, role, note,
               "positive_pairs_excluded", text_num(row$positive_pairs_excluded))
  add_log_stat(row, "MEAN_LOG_RATIO", "Multiplicative-scale agreement", "mean log ratio", row$mean_log_ratio, location, role, note)
  add_log_stat(row, "MEDIAN_LOG_RATIO", "Multiplicative-scale agreement", "median log ratio", row$median_log_ratio, location, role, note)
  add_log_stat(row, "GEOMETRIC_MEAN_RATIO", "Multiplicative-scale agreement", "geometric mean ratio", row$geometric_mean_ratio, location, role, note)
  add_log_stat(row, "LOG_LOA_LOWER", "Multiplicative-scale agreement", "lower 95% limit on log scale", row$log_loa_lower95, location, role, note)
  add_log_stat(row, "LOG_LOA_UPPER", "Multiplicative-scale agreement", "upper 95% limit on log scale", row$log_loa_upper95, location, role, note)
  add_log_stat(row, "RATIO_LOA_LOWER", "Multiplicative-scale agreement", "lower exponentiated 95% limit", row$ratio_loa_lower95, location, role, note)
  add_log_stat(row, "RATIO_LOA_UPPER", "Multiplicative-scale agreement", "upper exponentiated 95% limit", row$ratio_loa_upper95, location, role, note)
  add_log_stat(row, "PAIR_MEAN_PEARSON", "Multiplicative-scale agreement", "log-ratio versus log-pair-mean Pearson correlation", row$log_ratio_pair_mean_pearson_r, location, role, note,
               "Spearman_rho", text_num(row$log_ratio_pair_mean_spearman_rho))
}

for (i in seq_len(nrow(imbalance))) {
  row <- imbalance[i, , drop = FALSE]
  cohort_short <- if (row$cohort_id == "PRIMARY_FINAL_TARGET") "PRIMARY" else "SENSITIVITY"
  if (row$variable == "delta_count") {
    add(paste0("P2B_COUNT_", cohort_short, "_MEAN"), "sampling opportunity context", "Source-specific measurement-count imbalance", row$cohort_definition, row$N, NA_real_,
        "descriptive source-count difference", "mean delta_count", row$mean, NA_real_, NA_real_, NA_real_, "sd_delta_count", text_num(row$sd),
        "EXTENSION_RESULT_LOCKED", imbalance_path, paste0("cohort_id=", row$cohort_id, "; variable=delta_count"), script_path("phase2b", "01_source_dependence_closure.R"),
        "new_extension", "Phase 2B", "Supplement; sampling-process context", "Measurement-count imbalance distribution", "Count imbalance is descriptive; it is not a confounder or causal adjustment variable.")
    add(paste0("P2B_COUNT_", cohort_short, "_DISTRIBUTION"), "sampling opportunity context", "Source-specific measurement-count imbalance", row$cohort_definition, row$N, NA_real_,
        "descriptive source-count difference", "median delta_count", row$median, row$q5, row$q95, NA_real_, "IQR", text_num(row$iqr),
        "EXTENSION_RESULT_LOCKED", imbalance_path, paste0("cohort_id=", row$cohort_id, "; variable=delta_count"), script_path("phase2b", "01_source_dependence_closure.R"),
        "new_extension", "Phase 2B", "Supplement; sampling-process context", "Measurement-count imbalance distribution", "lower95/upper95 fields are q5/q95, not confidence limits.")
    add(paste0("P2B_COUNT_", cohort_short, "_GV_CORRELATION"), "sampling opportunity context", "Source-defined GV disagreement versus measurement-count imbalance", row$cohort_definition, row$N, NA_real_,
        "descriptive source-defined relationship", "Pearson correlation", row$delta_gv_delta_count_pearson_r, NA_real_, NA_real_, NA_real_, "Spearman_rho", text_num(row$delta_gv_delta_count_spearman_rho),
        "EXTENSION_RESULT_LOCKED", imbalance_path, paste0("cohort_id=", row$cohort_id, "; variable=delta_count"), script_path("phase2b", "01_source_dependence_closure.R"),
        "new_extension", "Phase 2B", "Supplement; sampling-process context", "Descriptive correlation only", "Count imbalance alone did not account for observed source-defined GV disagreement; this does not rule out unmeasured sampling mechanisms.")
  } else {
    add(paste0("P2B_SPAN_", cohort_short, "_NOT_ESTIMABLE"), "sampling opportunity context", "Source-specific measurement-span imbalance", row$cohort_definition, row$N, NA_real_,
        "declared source-span comparison", "NOT_ESTIMABLE", NA_real_, NA_real_, NA_real_, NA_real_, "unavailable_reason", as.character(row$unavailable_reason),
        "NOT_ESTIMABLE", imbalance_path, paste0("cohort_id=", row$cohort_id, "; variable=delta_span"), script_path("phase2b", "01_source_dependence_closure.R"),
        "new_extension", "Phase 2B", "Supplement; declared non-estimable context", "Boundary condition", "Validated final same-patient source file has source counts but no source-specific span fields; no undocumented reconstruction performed.")
  }
}

split_origin <- function(origin) {
  pos <- regexpr("#", as.character(origin), fixed = TRUE)[1]
  if (pos < 0L) return(c(path = as.character(origin), key = NA_character_))
  c(path = substr(as.character(origin), 1L, pos - 1L), key = substr(as.character(origin), pos + 1L, nchar(as.character(origin))))
}
for (i in seq_len(nrow(landscape))) {
  row <- landscape[i, , drop = FALSE]
  origin <- split_origin(row$estimate_origin)
  add(
    paste0("LANDSCAPE_", row$context_id), paste0(row$database, " analytic-context landscape"), row$context_id, row$cohort_definition,
    row$N, row$events, row$outcome, row$effect_type, row$estimate, row$lower95, row$upper95, row$p_value,
    "same_patient_pair", as.character(row$same_patient_pair), ifelse(row$fit_status == "LOCKED_RESULT_VERIFIED", "LOCKED_RESULT_VERIFIED", row$fit_status),
    origin[["path"]], origin[["key"]], phase2b_script_landscape, "locked_landscape_row", "Phase 2B",
    paste0("Candidate Figure 3; ", row$figure_panel), "Context-family landscape result",
    paste(row$interpretive_note, "Effect type is preserved as HR or RR; no pooling or cross-database heterogeneity test.")
  )
}

registry <- registry[order(registry$evidence_id), , drop = FALSE]
registry_path <- out("evidence_lock", "master_evidence_registry.csv")
phase2a_write(registry, registry_path)

locked_landscape_path <- out("evidence_lock", "analytic_context_landscape_LOCKED.csv")
if (!file.copy(landscape_path, locked_landscape_path, overwrite = TRUE)) stop("Could not copy analytic-context landscape lock")

main_ids <- c(
  "MIMIC_PRIMARY_COHORT_30D", "MIMIC_MEAN_GV", "MIMIC_MEAN_GLUCOSE", "MIMIC_MEASUREMENT_COUNT_MEDIAN", "MIMIC_MEASUREMENT_SPAN_MEDIAN",
  "MIMIC_MODEL_A_30D", "MIMIC_MODEL_B_30D", "MIMIC_MODEL_C_30D", "MIMIC_SAMEPATIENT_POCT_HR", "MIMIC_SAMEPATIENT_LAB_HR",
  "P2A_BOOTSTRAP_OBSERVED_DELTA_BETA", "P2A_BOOTSTRAP_COEFFICIENT_RATIO", "P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE",
  "P2A_AGREEMENT_PRIMARY_GV_PEARSON", "P2A_AGREEMENT_PRIMARY_GV_SPEARMAN", "P2B_LOG_PRIMARY_GEOMETRIC_MEAN_RATIO",
  "P2B_LOG_PRIMARY_RATIO_LOA_LOWER", "P2B_LOG_PRIMARY_RATIO_LOA_UPPER", "INSPIRE_I2_24H_30D", "INSPIRE_I3_24H_30D",
  "INSPIRE_CORRECTED_48H_30D", "EICU_M3", "EICU_M4"
)
landscape_ids <- registry$evidence_id[grepl("^LANDSCAPE_", registry$evidence_id)]
main_ids <- unique(c(main_ids, landscape_ids, "MIMIC_SOURCE_POCT_ONLY_30D", "MIMIC_SOURCE_CENTRAL_LAB_ONLY_30D", "MIMIC_SOURCE_BLOOD_GAS_ONLY_30D", "MIMIC_SOURCE_COMMON_30D"))

map_rows <- lapply(registry$evidence_id, function(id) {
  is_main <- id %in% main_ids
  data.frame(
    evidence_id = id,
    priority = ifelse(is_main, 1L, 2L),
    destination = ifelse(is_main, ifelse(grepl("^LANDSCAPE_", id), "Results; Candidate Figure 3", "Results; revised central evidence"), "Supplementary Results / evidence tables"),
    figure_or_table = ifelse(is_main && grepl("^LANDSCAPE_", id), "Candidate Figure 3", ifelse(is_main, "Candidate Figure 2 or main Results table", "Supplementary evidence table / audit")),
    main_or_supplement = ifelse(is_main, "main", "supplement"),
    reason = ifelse(is_main,
                    "Directly supports the revised measurement-context framing or provides the locked context-family display.",
                    "Retain for transparency, lineage reconciliation, boundary conditions, or sensitivity context; not central to the revised claim."),
    stringsAsFactors = FALSE
  )
})
evidence_map <- do.call(rbind, map_rows)
evidence_map_path <- out("evidence_lock", "manuscript_evidence_map.csv")
phase2a_write(evidence_map, evidence_map_path)

hash_string <- function(paths) {
  paths <- as.character(paths)
  paste(vapply(paths, function(p) paste0(p, "=", sha256(p)), character(1)), collapse = ";")
}
candidate_a <- out("figures", "candidate", "phase2b_candidate_figure_A_source_dependence.png")
candidate_b <- out("figures", "candidate", "phase2b_candidate_figure_B_analytic_context_landscape.png")
samepatient_input <- cfg_path("mimic_samepatient_source")
manuscript_root <- cfg_path("manuscript_provenance_root")
figure_map <- data.frame(
  figure = c("Figure 1", "Figure 2", "Figure 2", "Figure 2", "Figure 2", "Figure 3", "Figure 3", "Figure 3", "Figure 3", "Figure 4"),
  panel = c("framework specification (not generated in Phase 2C)", "A: same-patient source scatter", "B: original-scale Bland-Altman", "C: multiplicative-scale sensitivity", "D: same-patient mortality coefficients", "A: MIMIC adjustment context", "B: MIMIC source context", "C: INSPIRE timing/adjustment context", "D: eICU adjustment context", "decision"),
  source_data_file = c(
    paste(manuscript_root, "JAHA_Data_Set_S1_Figure_Source_Map.csv", sep = "/"),
    paste(samepatient_input, agreement_path, sep = ";"),
    paste(samepatient_input, agreement_path, sep = ";"),
    paste(samepatient_input, log_path, sep = ";"),
    paste(cfg_path("mimic_source_results"), source_model_path, bootstrap_path, sep = ";"),
    landscape_path, landscape_path, landscape_path, landscape_path,
    paste(file.path(manuscript_root, "JAHA_Data_Set_S1_Figure_Source_Map.csv"), file.path(cfg_path("analysis_record_root"), "results", "12_absolute_risk_results.csv"), sep = ";")
  ),
  analysis_id = c("FRAMEWORK_SPECIFICATION", "P2A_AGREEMENT_PRIMARY_GV", "P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE", "P2B_LOG_PRIMARY_GEOMETRIC_MEAN_RATIO", "MIMIC_SAMEPATIENT_POCT_HR;MIMIC_SAMEPATIENT_LAB_HR;P2A_BOOTSTRAP_COEFFICIENT_RATIO", "MIMIC_ADJUSTMENT_CONTEXT", "MIMIC_SOURCE_CONTEXT", "INSPIRE_CONTEXT", "EICU_CONTEXT", "ORIG_MIMIC_ABSOLUTE_RISK_30D"),
  N = c(NA, 452, 452, 447, 409, NA, NA, NA, NA, 9751),
  events = c(NA, NA, NA, NA, 49, NA, NA, NA, NA, 240),
  statistic = c(
    "Database roles, anchors, exposure windows, source and outcome differences specified for Phase 3 figure construction.",
    "POCT/laboratory GV scatter; Pearson and Spearman agreement.",
    "Mean difference and original-scale 95% limits of agreement.",
    "Mean/median log ratio, geometric ratio and exponentiated limits; positive pairs only.",
    "POCT HR, laboratory HR, paired delta-beta ratio and bootstrap interval; identical patients/events marked.",
    "Model A/B/C in conceptual order; HR; N/events shown per row.",
    "Priority series, source-specific series and same-patient POCT/laboratory pair; HR; pair visually distinct.",
    "I1/I2/I3 and corrected 48-hour landmark; HR; no withdrawn 365-day result.",
    "M1/M2/M3/M4 context family; RR; no pooling with HR.",
    "Original Figure 3 risk contrast; decision is MOVE_TO_SUPPLEMENT."
  ),
  source_status = c("SPECIFICATION_ONLY", "PHASE2A_LOCKED", "PHASE2A_LOCKED", "PHASE2B_LOCKED", "PHASE2A_LOCKED", "PHASE2B_LOCKED", "PHASE2B_LOCKED", "PHASE2B_LOCKED", "PHASE2B_LOCKED", "ORIGINAL_SUBMISSION_PROVENANCE"),
  script = c(NA, phase2b_script_figures, phase2b_script_figures, phase2b_script_figures, phase2b_script_figures, phase2b_script_figures, phase2b_script_figures, phase2b_script_figures, phase2b_script_figures, file.path(manuscript_root, "JAHA_Data_Set_S1_Figure_Source_Map.csv")),
  hash = c(
    hash_string(c(file.path(manuscript_root, "JAHA_Data_Set_S1_Figure_Source_Map.csv"), out("machine_readable", "final_measurement_context.csv"))),
    hash_string(c(samepatient_input, agreement_path)),
    hash_string(c(samepatient_input, agreement_path)),
    hash_string(c(samepatient_input, log_path)),
    hash_string(c(cfg_path("mimic_source_results"), source_model_path, bootstrap_path)),
    hash_string(landscape_path), hash_string(landscape_path), hash_string(landscape_path), hash_string(landscape_path),
    hash_string(c(file.path(manuscript_root, "JAHA_Data_Set_S1_Figure_Source_Map.csv"), file.path(cfg_path("analysis_record_root"), "results", "12_absolute_risk_results.csv")))
  ),
  stringsAsFactors = FALSE
)
figure_map_path <- out("evidence_lock", "figure_source_map_LOCKED.csv")
phase2a_write(figure_map, figure_map_path)

brief_path <- file.path(phase2c_workspace, "PHASE3_MANUSCRIPT_REWRITE_BRIEF.md")
brief_lines <- c(
  "# Phase 3 Manuscript Rewrite Brief",
  "",
  "Status: evidence locked in Phase 2C; this file is a rewrite specification, not a manuscript edit.",
  "",
  "## Recommended title",
  "",
  "**Measurement-Context Dependence of Routine-Care Glycemic Variability After Cardiac Surgery: A Multidatabase Cohort Study**",
  "",
  "The alternative title, *Routine-Care Glycemic Variability After Cardiac Surgery: Dependence on Measurement Source and Analytic Context*, is accurate but less direct about the multidatabase design and is not preferred. The phrase ‘sampling opportunity’ is not placed in the title because source-specific span was not estimable and count imbalance alone did not account for the observed disagreement.",
  "",
  "## Central claim",
  "",
  "In routine-care cardiac-surgery data, SD-based glycemic variability was measurement-context dependent: measurement source materially altered both derived GV values and the estimated mortality association in identical patients, whereas observed measurement-count differences alone did not account for the source-defined disagreement and source-specific span was not estimable.",
  "",
  "This is an observational, noncausal claim. It does not designate POCT as wrong, laboratory measurements as true, measurement count/span as confounders, source coefficients as causal interactions, or database differences as biological heterogeneity.",
  "",
  "## Manuscript structure",
  "",
  "### Introduction",
  "",
  "1. Establish the clinical interest in postoperative glycemic variability while clarifying that routine-care GV is a derived biomarker whose value depends on the observed glucose series, source composition, timing anchor, exposure window, and outcome-risk origin.",
  "2. Summarize why existing GV literature can be heterogeneous when routine-care glucose is sampled and summarized differently.",
  "3. State the measurement-source and analytic-context problem: a mortality coefficient can change when the observed measurement context changes, even when the patient cohort is held constant.",
  "4. Define the objective as a bounded measurement-context assessment across final MIMIC, INSPIRE, and eICU lineage contexts, not an exhaustive multiverse.",
  "",
  "### Results",
  "",
  "1. Start with the final context table: MIMIC N=10,561/296 at 30 days, INSPIRE N=1,353/27, corrected 48-hour N=1,511/31, and eICU M3/M4 N=7,115/130; state each anchor, exposure window, source structure, and HR versus RR estimand.",
  "2. Make same-patient source dependence the first empirical section: N=452 primary paired patients, N=453 sensitivity candidates; original-scale agreement, scale dependence, and the identical N=409/49 mortality comparison.",
  "3. Report the multiplicative-scale analysis as a sensitivity to the original-scale Bland-Altman result. It does not replace the original scale and its limits remain broad.",
  "4. Report measurement-count imbalance as a limited descriptive closure: count relationships are modest and source-specific span is NOT_ESTIMABLE. Do not claim that sampling intensity explains the association.",
  "5. Present the analytic-context landscape in conceptual family order: MIMIC adjustment, MIMIC source, INSPIRE timing/adjustment, and eICU adjustment. Preserve HR/RR labels and do not pool.",
  "",
  "### Discussion",
  "",
  "1. Lead with the revised interpretation: source-defined routine-care GV behaves as a measurement-context-dependent biomarker.",
  "2. Explain why identical-patient source coefficients are stronger evidence for measurement-context dependence than cross-database coefficient comparison, while remaining observational and noncausal.",
  "3. Bound the sampling interpretation: count imbalance alone was weakly related to GV disagreement; source-specific span was not estimable; no single observed process fully explains heterogeneity.",
  "4. State database limitations, eICU aggregate-only boundaries, INSPIRE day-30 coverage attestation, event-limited INSPIRE estimates, and the non-equivalence of RR and HR.",
  "5. Explain implications for interpreting and reporting routine-care GV without treating it as a gold-standard biological exposure.",
  "6. State clinical implications cautiously: context should accompany GV reporting and interpretation, while no causal treatment recommendation follows from these observational associations.",
  "7. State limitations, including routine-care sampling, source composition, non-estimable source-specific span, aggregate-only eICU inputs, event-limited INSPIRE estimates, and lack of a common cross-database estimand.",
  "8. Close with future protocolized sampling and CGM validation as the appropriate route for testing transportability of GV as a measurement construct.",
  "9. Relegate 365-day MIMIC average-HR material to the supplement because proportional hazards were violated; do not reintroduce withdrawn INSPIRE 365-day results.",
  "",
  "## Main manuscript versus supplement",
  "",
  "- Main: final context/estimand table; Figure 2 candidate same-patient source dependence; the identical-patient mortality coefficients and paired delta-beta sensitivity; Candidate Figure 3 displaying all 18 context-family rows in conceptual order; concise mean-glucose, count, span, source composition, exposure-window, anchor, and risk-origin descriptions.",
  "- Supplement: sensitivity N=453 agreement; full log-scale values; all sampling-count correlations and NOT_ESTIMABLE span rows; specification movements; the full machine-readable/tabular 18-row landscape; source-restricted cohort details; alternative GV metrics; MIMIC 365-day PH/time-varying material; MICE diagnostics and model sequences; the original mean-glucose/absolute-risk contrast (MOVE_TO_SUPPLEMENT); full INSPIRE boundary audit.",
  "- Figure 1 is a Phase 3 framework schematic only at this lock: MIMIC date-anchored perioperative context, INSPIRE exact operation-end timestamp, eICU ICU-admission/aggregate context, and source/sampling/outcome differences. It is not generated or locked as a bitmap in Phase 2C.",
  "- Candidate Figure 2 and Candidate Figure 3 remain labeled CANDIDATE until manuscript production and rendered-artifact QA.",
  "",
  "## Reporting minimum",
  "",
  "For every displayed GV result, state the GV definition, mean glucose context, measurement count and span, source composition, exposure window, time anchor, landmark/risk origin, outcome definition, effect type, and whether the sampling structure reflects protocolized sampling or routine care. These are minimum reporting elements suggested by these findings, not claims of formal external standards.",
  "",
  "## Language blacklist",
  "",
  "Avoid: ‘GV has no prognostic value’; ‘GV is invalid’; ‘measurement error caused previous findings’; ‘POCT is inaccurate’; ‘laboratory GV is the true GV’; ‘sampling frequency confounded the association’; ‘sampling intensity explains the GV association’; ‘count/span are confounders’; ‘independent effect’ when the estimand is incremental association; ‘causal effect’; ‘causal interaction’; ‘mediation’; ‘device effect’; ‘all heterogeneity is explained by measurement process’; ‘biological heterogeneity across databases’; ‘external validation’ when the estimand is not transport-equivalent; and any pooled/common summary effect across databases.",
  "",
  "## Evidence anchors",
  "",
  "Primary evidence IDs: MIMIC_MODEL_B_30D; MIMIC_SAMEPATIENT_POCT_HR; MIMIC_SAMEPATIENT_LAB_HR; P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE; P2A_BOOTSTRAP_OBSERVED_DELTA_BETA; P2B_LOG_PRIMARY_GEOMETRIC_MEAN_RATIO; P2B_COUNT_PRIMARY_GV_CORRELATION; P2B_SPAN_PRIMARY_NOT_ESTIMABLE; and the LANDSCAPE_* rows in the locked registry."
)
writeLines(brief_lines, brief_path, useBytes = TRUE)

report_path <- file.path(phase2c_workspace, "PHASE2C_EVIDENCE_LOCK_REPORT.md")
report_num <- function(x, digits = 4L) {
  z <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(z)) return("NA")
  formatC(z, format = "f", digits = digits, big.mark = ifelse(digits == 0L, ",", ""))
}
report_row <- function(data, column, value) {
  hit <- data[as.character(data[[column]]) == as.character(value), , drop = FALSE]
  if (nrow(hit) != 1L) stop("Report source row is not unique: ", column, "=", value)
  hit[1, , drop = FALSE]
}
report_agreement <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "GV_SD", , drop = FALSE]
if (nrow(report_agreement) != 1L) stop("Report source agreement row is not unique")
report_log <- report_row(log_scale, "cohort_id", "PRIMARY_FINAL_TARGET")
report_imbalance <- imbalance[imbalance$cohort_id == "PRIMARY_FINAL_TARGET" & imbalance$variable == "delta_count", , drop = FALSE]
report_shift_mimic <- report_row(shift, "database", "MIMIC-IV")
report_shift_eicu <- report_row(shift, "database", "eICU-CRD")
report_shift_inspire <- report_row(shift, "database", "INSPIRE")
report_source_poct <- report_row(source_results, "model_id", "SRC_samepatient_poct_30d")
report_source_lab <- report_row(source_results, "model_id", "SRC_samepatient_lab_30d")
report_i2 <- report_row(inspire, "model_id", "ADMINV5_I2_30d")
report_i2_48 <- report_row(inspire48, "model_id", "ADMINV5_I2_48h_landmark_30d")
report_eicu_m3 <- report_row(eicu, "model_id", "HARM_eICU-CRD_M3")
report_eicu_m4 <- report_row(eicu, "model_id", "HARM_eICU-CRD_M4")
report_lines <- c(
  "# PHASE2C Evidence Lock Report",
  "",
  "Phase 2C assembled a provenance-bearing evidence package from the final Phase 1.6 lineage and locked Phase 2A/2B outputs. No new inferential analysis was performed in this phase.",
  "",
  "## Material Passport",
  "",
  "- Workspace: `JAHA_v5_analysis_extensions`.",
  "- Upstream gates: `GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES`, `GO_PHASE2B_ANALYTIC_CONTEXT_LANDSCAPE`, and the Phase 2B `GO_PHASE3_MANUSCRIPT_REWRITE` gate were required before assembly.",
  "- Primary registry: `outputs/evidence_lock/master_evidence_registry.csv`.",
  "- Locked landscape: `outputs/evidence_lock/analytic_context_landscape_LOCKED.csv`.",
  "- Figure source lock: `outputs/evidence_lock/figure_source_map_LOCKED.csv`.",
  "- Reproducibility manifest: `outputs/evidence_lock/evidence_lock_manifest.json`.",
  "- The submitted manuscript and controlled inputs were read-only provenance sources; they were not modified.",
  "",
  "## 1. Can every proposed main numerical claim be traced to a locked result or an extension result?",
  "",
  "Yes, subject to the Phase 2C QC gate. The master registry assigns a unique evidence ID, population, estimand/effect type, numeric fields, source file, source row/key, script origin, lineage status, and manuscript location to every proposed main numerical claim. The primary traceability chain is: MIMIC Model A/B/C and final context summary; same-patient Phase 2A agreement and source-model reproduction; paired bootstrap; Phase 2B log-scale and count-imbalance outputs; and the 18-row Phase 2B analytic-context landscape. The original mean-glucose/absolute-risk result is also registered but explicitly marked supplementary.",
  "",
  "## 2. Are the proposed main claims reproducible from frozen scripts and inputs?",
  "",
  "Yes for the locked evidence package. The manifest records SHA-256 hashes for the Phase 2A/2B machine outputs, candidate figures, scripts, final-lineage inputs, submitted provenance files, registries, maps, and rewrite brief. Reproduction here means re-reading and verifying the locked artifacts; Phase 2C does not rerun models. Candidate figure bitmap regeneration remains a Phase 3 production task because the submitted figure builders are not present in the public-backup repository.",
  "",
  "## 3. What are the core 5–10 results that deserve main-manuscript emphasis?",
  "",
  paste0("1. Final context identity and estimands: MIMIC primary N=", report_num(reconstruction$mimic$target_n, 0), "/", report_num(reconstruction$mimic$events_30d, 0), "; INSPIRE I2 N=", report_num(report_i2$N, 0), "/", report_num(report_i2$events, 0), " and corrected 48-hour N=", report_num(report_i2_48$N, 0), "/", report_num(report_i2_48$events, 0), "; eICU M3/M4 N=", report_num(report_eicu_m3$N, 0), "/", report_num(report_eicu_m3$events, 0), ". Source IDs: `MIMIC_PRIMARY_COHORT_30D`, `INSPIRE_I2_24H_30D`, `INSPIRE_CORRECTED_48H_30D`, `EICU_M3`, `EICU_M4`.", sep = ""),
  paste0("2. Same-patient GV agreement on the original scale: POCT mean ", report_num(report_agreement$poct_mean), ", laboratory mean ", report_num(report_agreement$lab_mean), ", Pearson r ", report_num(report_agreement$pearson_r), ", Spearman rho ", report_num(report_agreement$spearman_rho), "; mean difference ", report_num(report_agreement$mean_paired_difference), " mg/dL with limits ", report_num(report_agreement$loa_lower95), " to ", report_num(report_agreement$loa_upper95), " mg/dL. Evidence IDs: `P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE`, `P2A_AGREEMENT_PRIMARY_GV_PEARSON`, `P2A_AGREEMENT_PRIMARY_GV_SPEARMAN`.", sep = ""),
  paste0("3. Scale dependence of disagreement: primary positive-pair effective N=", report_num(report_log$N, 0), ", mean log ratio ", report_num(report_log$mean_log_ratio), ", median ", report_num(report_log$median_log_ratio), ", geometric mean ratio ", report_num(report_log$geometric_mean_ratio), "; exponentiated limits ", report_num(report_log$ratio_loa_lower95), " to ", report_num(report_log$ratio_loa_upper95), "; log-ratio versus log-pair-mean Pearson/Spearman ", report_num(report_log$log_ratio_pair_mean_pearson_r), "/", report_num(report_log$log_ratio_pair_mean_spearman_rho), ". Evidence IDs: `P2B_LOG_PRIMARY_*`; this sensitivity does not replace the original scale.", sep = ""),
  paste0("4. Identical-patient mortality coefficients: POCT HR ", report_num(report_source_poct$HR_per10), " (N=409/49) versus laboratory HR ", report_num(report_source_lab$HR_per10), " (N=409/49). Evidence IDs: `MIMIC_SAMEPATIENT_POCT_HR` and `MIMIC_SAMEPATIENT_LAB_HR`; the identical N/events comparison is marked in the candidate figure source lock.", sep = ""),
  paste0("5. Paired coefficient-difference sensitivity: delta-beta ", report_num(boot$observed_delta_beta), " (95% percentile interval ", report_num(boot$percentile_ci_delta_beta_lower), " to ", report_num(boot$percentile_ci_delta_beta_upper), "), exponentiated ratio ", report_num(boot$observed_coefficient_ratio), ", bootstrap median delta-beta ", report_num(boot$bootstrap_median_delta_beta), "; B=", report_num(boot$requested_replicates, 0), ", successful=", report_num(boot$successful_replicates, 0), ", failed=", report_num(boot$failed_replicates, 0), ". Evidence IDs: `P2A_BOOTSTRAP_OBSERVED_DELTA_BETA` and `P2A_BOOTSTRAP_COEFFICIENT_RATIO`; source-definition descriptive evidence, not a causal interaction.", sep = ""),
  paste0("6. MIMIC adjustment context: `MIMIC_MODEL_A_30D`, `MIMIC_MODEL_B_30D`, and `MIMIC_MODEL_C_30D` in conceptual order; the final primary Model B is HR ", report_num(mice_b30$HR_per10), " (", report_num(mice_b30$lo_per10), "–", report_num(mice_b30$hi_per10), "), P=", report_num(mice_b30$P_per10), ", with Model C as the locked context shift.", sep = ""),
  paste0("7. Source-specific MIMIC context: priority series, POCT-only, central-laboratory-only, blood-gas-only, common-source, and the same-patient pair in the `MIMIC_SOURCE_*` rows. The count-imbalance closure shows delta-count mean ", report_num(report_imbalance$mean), " (SD ", report_num(report_imbalance$sd), "), median ", report_num(report_imbalance$median, 0), ", IQR ", report_num(report_imbalance$iqr, 0), ", q5–q95 ", report_num(report_imbalance$q5, 0), "–", report_num(report_imbalance$q95, 0), "; source-defined GV disagreement Pearson/Spearman ", report_num(report_imbalance$delta_gv_delta_count_pearson_r), "/", report_num(report_imbalance$delta_gv_delta_count_spearman_rho), ". Evidence ID: `P2B_COUNT_PRIMARY_GV_CORRELATION`; source-specific span is `P2B_SPAN_PRIMARY_NOT_ESTIMABLE`.", sep = ""),
  paste0("8. Analytic-context landscape and context shifts: all `LANDSCAPE_*` rows are separated into MIMIC adjustment, MIMIC source, INSPIRE timing/adjustment, and eICU adjustment families with HR/RR preserved. Locked shifts are MIMIC B→C delta log effect ", report_num(report_shift_mimic$delta_log_effect), ", eICU M3→M4 ", report_num(report_shift_eicu$delta_log_effect), ", and INSPIRE I2→I3 ", report_num(report_shift_inspire$delta_log_effect), "; no pooling or heterogeneity test.", sep = ""),
  "",
  "## 4. Which analyses belong in the supplement?",
  "",
  "The supplement should contain the N=453 paired sensitivity; full original- and log-scale agreement details; mean-glucose paired context; all count correlations and the source-span NOT_ESTIMABLE rows; the three specification-shift diagnostics; the full machine-readable/tabular 18-row landscape that underlies Candidate Figure 3; source-restricted cohort details; alternative GV metrics; MIMIC 365-day PH and time-varying material; MICE/model diagnostics; the original mean-glucose/absolute-risk contrast; and the full INSPIRE coverage/timing boundary audit. The evidence map records this assignment row by row.",
  "",
  "## 5. Which prior main-manuscript results should be demoted?",
  "",
  "The original mean-glucose conditioning/absolute-risk Figure 3 result (`ORIG_MIMIC_ABSOLUTE_RISK_30D`) should move to the supplement: it is valid secondary evidence but is not the revised story's central estimand. The MIMIC 365-day average-HR sequence should move to the supplement and carry the documented PH violation (`ORIG_MIMIC_365D_PH_DIAGNOSTIC`; `MIMIC_MODEL_B_365D`). The old cross-database non-pooled display should be refocused into the conceptual landscape rather than treated as a common validation effect. Withdrawn INSPIRE 365-day and historical coverage-gate results remain excluded.",
  "",
  "## 6. Should source dependence be the centerpiece?",
  "",
  "Yes. The same-patient comparison provides the cleanest evidence that changing the observed measurement source can change derived GV and its mortality coefficient while holding the patient/event cohort fixed. The interpretation remains measurement-context dependence, not source accuracy and not a causal source-by-exposure interaction.",
  "",
  "## 7. Is ‘sampling opportunity’ too strong for the title?",
  "",
  "Yes. It is too strong as a title-level causal-sounding explanation because source-specific span was not estimable and measurement-count imbalance alone did not account for source-defined GV disagreement. Sampling opportunity should remain a bounded context dimension in the Results and Supplement, not the title's headline mechanism.",
  "",
  "## Findings that could weaken the framing",
  "",
  "The broad original-scale limits of agreement and broad multiplicative ratio limits (`P2A_AGREEMENT_PRIMARY_GV_MEAN_DIFFERENCE`; `P2B_LOG_PRIMARY_RATIO_LOA_LOWER`; `P2B_LOG_PRIMARY_RATIO_LOA_UPPER`) limit claims of precise interchangeability. Source-specific span is not estimable (`P2B_SPAN_PRIMARY_NOT_ESTIMABLE`), the count closure is descriptive rather than explanatory (`P2B_COUNT_PRIMARY_GV_CORRELATION`), INSPIRE estimates are event-limited, and eICU preserves aggregate rather than event-level glucose inputs. These findings weaken generalization and mechanism claims, but they do not materially overturn the narrower measurement-context framing.",
  "",
  "## 8. What exact title is recommended?",
  "",
  "**Measurement-Context Dependence of Routine-Care Glycemic Variability After Cardiac Surgery: A Multidatabase Cohort Study**",
  "",
  "## 9. What exact one-sentence central claim is recommended?",
  "",
  "In routine-care cardiac-surgery data, SD-based glycemic variability was measurement-context dependent: measurement source materially altered both derived GV values and the estimated mortality association in identical patients, whereas observed measurement-count differences alone did not account for the source-defined disagreement and source-specific span was not estimable.",
  "",
  "## 10. What contradictions exist between the old manuscript and the revised evidence package?",
  "",
  "- The old primary null/mean-glucose framing is not numerically invalid for final MIMIC Model B, but it is no longer the central interpretation; it is subordinate to the source-dependence evidence.",
  "- The old absolute-risk Figure 3 remains a valid complete-case secondary contrast, but its manuscript role changes from main figure to supplement.",
  "- The old 365-day average-HR material remains an audit result only because final provenance documents PH violation; it must not be presented as an uncomplicated long-horizon primary effect.",
  "- The submitted Figure 5 provenance included a harmonized MIMIC comparison frame of N=8,117/128, whereas the final Phase 2B landscape uses the final MIMIC primary/source contexts and the locked eICU M3/M4 N=7,115/130 frame. These are different non-pooled analytic contexts, not a single cross-database estimand and not a biological contradiction.",
  "- The INSPIRE final I1/I2/I3 and corrected 48-hour results are retained with the final v5 administrative-censoring and coverage-attestation boundaries; withdrawn 365-day results and historical wide-mask/discharge-censoring outputs are not part of the revised evidence.",
  "- The candidate figures are provenance-locked but remain candidate production artifacts; missing public-backup figure builders are a Phase 3 production task, not a license to introduce manually typed numbers.",
  "",
  "## Lock interpretation",
  "",
  "The total evidence supports ‘measurement-context dependence’ as the revised manuscript's central framing. It does not support the stronger claim that sampling intensity explains the GV association, a reference-standard hierarchy between POCT and laboratory measurements, causal source interactions, database pooling, or biological heterogeneity claims.",
  "",
  "GO_PHASE3_MANUSCRIPT_REWRITE"
)
writeLines(report_lines, report_path, useBytes = TRUE)

generated_paths <- c(
  registry_path, locked_landscape_path, evidence_map_path, figure_map_path, brief_path, report_path,
  landscape_prov_path, log_path, imbalance_path, out("machine_readable", "phase2b_source_sampling_imbalance_summary.csv"),
  agreement_path, bootstrap_path, sampling_path, shift_path, source_model_path
)
phase2a_outputs <- c(
  agreement_path, out("machine_readable", "phase2a_same_patient_bootstrap.csv"), bootstrap_path, sampling_path,
  shift_path, source_model_path, out("machine_readable", "phase2a_measurement_context_table.csv"), out("machine_readable", "phase2_estimand_inventory.csv")
)
phase2b_outputs <- c(landscape_path, log_path, imbalance_path, out("machine_readable", "phase2b_source_sampling_imbalance_summary.csv"), landscape_prov_path)
candidate_figures <- c(candidate_a, candidate_b)
phase2a_scripts <- list.files(file.path(phase2c_workspace, "scripts", "phase2a"), pattern = "\\.(R|r)$", full.names = TRUE)
phase2b_scripts <- list.files(file.path(phase2c_workspace, "scripts", "phase2b"), pattern = "\\.(R|r)$", full.names = TRUE)
phase2c_scripts <- list.files(file.path(phase2c_workspace, "scripts", "phase2c"), pattern = "\\.(R|r)$", full.names = TRUE)
final_lineage_inputs <- c(
  cfg_path("mimic_analysis_base"), cfg_path("mimic_features_priority"), cfg_path("mimic_series_priority"), cfg_path("mimic_samepatient_source"),
  cfg_path("mimic_primary_results"), cfg_path("mimic_mice_results"), cfg_path("mimic_source_results"),
  cfg_path("eicu_aggregate_input"), cfg_path("eicu_locked_results"), cfg_path("inspire_base"), cfg_path("inspire_outcome_30d"),
  cfg_path("inspire_glucose_features_anchor"), cfg_path("inspire_cohort_48h"), cfg_path("inspire_primary_results"), cfg_path("inspire_48h_results"),
  cfg_path("inspire_time_rule_qc"), cfg_path("inspire_analysis_manifest"), cfg_path("public_primary_model_script"),
  cfg_path("public_source_builder_script"), cfg_path("public_source_model_script"), cfg_path("public_eicu_model_script"),
  cfg_path("controlled_inspire_runner"), file.path(cfg_path("analysis_record_root"), "results", "12_absolute_risk_results.csv"),
  file.path(cfg_path("analysis_record_root"), "results", "07_ph_diagnostics.csv"), file.path(cfg_path("analysis_record_root"), "results", "PH365_INTERVAL_MICE_POOLED.csv")
)
submitted_inputs <- list.files(cfg_path("manuscript_provenance_root"), pattern = "^(JAHA_Data_Set_S[123].*\\.csv|JAHA_Revised_Manuscript\\.(docx|pdf)|JAHA_Supplemental_Material\\.(docx|pdf))$", full.names = TRUE)
manifest_paths <- unique(c(
  generated_paths, phase2a_outputs, phase2b_outputs, candidate_figures, phase2a_scripts, phase2b_scripts, phase2c_scripts,
  script_path("00_environment_check.py"), script_path("01_build_final_measurement_context.py"), script_path("02_build_phase2_estimand_inventory.py"),
  script_path("03_final_lineage_qc.py"), file.path(phase2c_workspace, "config.yaml"), final_lineage_inputs, submitted_inputs
))
manifest_paths <- manifest_paths[file.exists(manifest_paths)]
file_entries <- lapply(manifest_paths, function(path) {
  kind <- if (path %in% generated_paths) "phase2c_generated_or_locked" else if (path %in% phase2a_outputs) "phase2a_machine_output" else if (path %in% phase2b_outputs) "phase2b_machine_output" else if (path %in% candidate_figures) "candidate_figure" else if (path %in% submitted_inputs) "submitted_manuscript_provenance" else if (path %in% final_lineage_inputs) "final_lineage_locked_input_or_script" else "workspace_script_or_config"
  list(path = path, sha256 = sha256(path), kind = kind, status = "present_at_evidence_lock")
})
public_root <- cfg_path("public_code_root")
public_commit <- tryCatch(system2("git", c("-C", public_root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_)
public_status <- tryCatch(system2("git", c("-C", public_root, "status", "--short"), stdout = TRUE, stderr = FALSE), error = function(e) character())
package_version_safe <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(packageVersion(pkg))
}
manifest <- list(
  phase = "Phase 2C — Evidence Lock Before Manuscript Rewrite",
  decision = "GO_PHASE3_MANUSCRIPT_REWRITE",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  workspace = phase2c_workspace,
  lineage = list(
    public_code_root = public_root,
    analysis_record_root = cfg_path("analysis_record_root"),
    controlled_package_root = cfg_path("controlled_package_root"),
    submitted_manuscript_provenance_root = cfg_path("manuscript_provenance_root"),
    historical_branch_policy = "Historical branches are provenance-only and are not evidence inputs."
  ),
  software = list(
    R = R.version$version.string,
    jsonlite = package_version_safe("jsonlite"),
    survival = package_version_safe("survival"),
    rms = package_version_safe("rms"),
    platform = paste(Sys.info()[c("sysname", "release", "machine")], collapse = "; "),
    sha256_tool = "shasum -a 256"
  ),
  reproducibility = list(
    random_seed = as.numeric(cfg$random_seed),
    bootstrap_replicates_locked_from_phase2a = 2000,
    new_inferential_analysis_in_phase2c = FALSE,
    no_new_model_fits = TRUE,
    manifest_sha256_scope = "All files listed under files; this manifest is excluded to avoid recursive hashing."
  ),
  public_repository = list(
    commit = public_commit,
    status_clean = length(public_status) == 0L,
    status_raw = public_status
  ),
  files = file_entries
)
manifest_path <- out("evidence_lock", "evidence_lock_manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, na = "null")

cat("PHASE2C_EVIDENCE_LOCK_ASSEMBLY_DONE\n")
