#!/usr/bin/env Rscript

phase2a_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2a_file, mustWork = TRUE)), "phase2a_common.R"))
suppressMessages({ library(jsonlite) })
cfg <- phase2a_load_config()
phase2a_open_log(cfg, "phase2a_04_model_specification_shift.log")
on.exit(phase2a_close_log(), add = TRUE)

if (!phase2a_phase1_gate(cfg)) stop("PHASE1_6_GATE_FAIL")

comparison_row <- function(database, comparison, effect_type, before, after, before_effect, after_effect, before_lo, before_hi, after_lo, after_hi, N_before, events_before, N_after, events_after, source_before, source_after) {
  data.frame(
    database = database,
    comparison = comparison,
    effect_type = effect_type,
    specification_before = before,
    specification_after = after,
    effect_before = before_effect,
    lower95_before = before_lo,
    upper95_before = before_hi,
    effect_after = after_effect,
    lower95_after = after_lo,
    upper95_after = after_hi,
    delta_log_effect = log(after_effect) - log(before_effect),
    null_crossing = ((before_effect - 1) * (after_effect - 1) <= 0) && before_effect != after_effect,
    N_before = N_before, events_before = events_before,
    N_after = N_after, events_after = events_after,
    source_before = source_before, source_after = source_after,
    verification_status = "LOCKED_RESULT_VERIFIED",
    interpretation_boundary = "sampling-process specification shift; not confounding attenuation",
    stringsAsFactors = FALSE
  )
}

mice <- phase2a_read(phase2a_cfg_path(cfg, "mimic_mice_results"))
eicu <- phase2a_read(phase2a_cfg_path(cfg, "eicu_locked_results"))
inspire <- phase2a_read(phase2a_cfg_path(cfg, "inspire_primary_results"))
mb <- mice[mice$model_id == "MICE_B_30d", , drop = FALSE]
mc <- mice[mice$model_id == "MICE_C_30d", , drop = FALSE]
em3 <- eicu[eicu$model_id == "HARM_eICU-CRD_M3", , drop = FALSE]
em4 <- eicu[eicu$model_id == "HARM_eICU-CRD_M4", , drop = FALSE]
ii2 <- inspire[inspire$model_id == "ADMINV5_I2_30d", , drop = FALSE]
ii3 <- inspire[inspire$model_id == "ADMINV5_I3_30d", , drop = FALSE]
if (any(vapply(list(mb, mc, em3, em4, ii2, ii3), nrow, integer(1)) != 1L)) stop("LOCKED_SPECIFICATION_RESULT_MISSING")

rows <- list(
  comparison_row("MIMIC-IV", "Model B -> Model C", "HR", "MIMIC Model B", "MIMIC Model C", mb$HR_per10, mc$HR_per10, mb$lo_per10, mb$hi_per10, mc$lo_per10, mc$hi_per10, mb$N, mb$events, mc$N, mc$events, phase2a_cfg_path(cfg, "mimic_mice_results"), phase2a_cfg_path(cfg, "mimic_mice_results")),
  comparison_row("eICU-CRD", "M3 -> M4", "RR", "eICU M3", "eICU M4", em3$RR_per10, em4$RR_per10, em3$lo, em3$hi, em4$lo, em4$hi, em3$N, em3$events, em4$N, em4$events, phase2a_cfg_path(cfg, "eicu_locked_results"), phase2a_cfg_path(cfg, "eicu_locked_results")),
  comparison_row("INSPIRE", "I2 -> I3", "HR", "INSPIRE I2", "INSPIRE I3", ii2$HR_per10, ii3$HR_per10, ii2$lo, ii2$hi, ii3$lo, ii3$hi, ii2$N, ii2$events, ii3$N, ii3$events, phase2a_cfg_path(cfg, "inspire_primary_results"), phase2a_cfg_path(cfg, "inspire_primary_results"))
)
shift <- do.call(rbind, rows)
shift$lineage_version <- "JAHA_v5_final_lineage"
phase2a_write(shift, phase2a_output_path(cfg, "machine_readable", "phase2a_specification_shift.csv"))
cat("PHASE2A_SPECIFICATION_SHIFT_DONE\n")
