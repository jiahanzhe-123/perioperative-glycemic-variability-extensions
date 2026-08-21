#!/usr/bin/env Rscript

phase2a_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2a_file, mustWork = TRUE)), "phase2a_common.R"))
suppressMessages({ library(jsonlite) })
cfg <- phase2a_load_config()
phase2a_open_log(cfg, "phase2a_01_same_patient_agreement.log")
on.exit(phase2a_close_log(), add = TRUE)

if (!phase2a_phase1_gate(cfg)) stop("PHASE1_6_GATE_FAIL")
data <- phase2a_source_data(cfg)
sp <- data$sp
target_ids <- unique(phase2a_target_ids(cfg))
primary <- sp[sp$stay_id %in% target_ids, , drop = FALSE]
candidate <- sp
if (nrow(primary) != 452L || nrow(candidate) != 453L) stop("SAME_PATIENT_DENOMINATOR_FAIL")
if (any(!is.finite(phase2a_num(sp$gv_poct))) || any(!is.finite(phase2a_num(sp$gv_lab)))) stop("SAME_PATIENT_GV_MISSING_FAIL")

agreement_row <- function(dd, cohort_id, cohort_definition, variable, poct_col, lab_col) {
  poct <- phase2a_num(dd[[poct_col]])
  lab <- phase2a_num(dd[[lab_col]])
  keep <- is.finite(poct) & is.finite(lab)
  poct <- poct[keep]; lab <- lab[keep]
  if (length(poct) < 2L) stop("Agreement cohort has fewer than two complete pairs")
  diff <- poct - lab
  pair_mean <- (poct + lab) / 2
  q_poct <- as.numeric(quantile(poct, c(.05, .25, .50, .75, .95), names = FALSE, type = 7))
  q_lab <- as.numeric(quantile(lab, c(.05, .25, .50, .75, .95), names = FALSE, type = 7))
  q_diff <- as.numeric(quantile(diff, c(.25, .50, .75), names = FALSE, type = 7))
  md <- mean(diff); sd_diff <- sd(diff)
  data.frame(
    cohort_id = cohort_id,
    cohort_definition = cohort_definition,
    N = length(poct),
    variable = variable,
    poct_mean = mean(poct), lab_mean = mean(lab),
    poct_sd = sd(poct), lab_sd = sd(lab),
    poct_median = q_poct[3], lab_median = q_lab[3],
    poct_iqr = q_poct[4] - q_poct[2], lab_iqr = q_lab[4] - q_lab[2],
    poct_iqr_q1 = q_poct[2], poct_iqr_q3 = q_poct[4],
    lab_iqr_q1 = q_lab[2], lab_iqr_q3 = q_lab[4],
    poct_q5 = q_poct[1], poct_q95 = q_poct[5],
    lab_q5 = q_lab[1], lab_q95 = q_lab[5],
    pearson_r = cor(poct, lab, method = "pearson"),
    spearman_rho = cor(poct, lab, method = "spearman"),
    mean_paired_difference = md,
    sd_paired_difference = sd_diff,
    median_paired_difference = q_diff[2],
    paired_difference_iqr = q_diff[3] - q_diff[1],
    paired_difference_iqr_q1 = q_diff[1], paired_difference_iqr_q3 = q_diff[3],
    loa_lower95 = md - 1.96 * sd_diff,
    loa_upper95 = md + 1.96 * sd_diff,
    difference_pair_mean_pearson_r = cor(diff, pair_mean, method = "pearson"),
    difference_pair_mean_spearman_rho = cor(diff, pair_mean, method = "spearman"),
    paired_difference_definition = "POCT - central-laboratory source value",
    analysis_status = "COMPLETED",
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (item in list(
  list(dd = primary, cohort = "PRIMARY_FINAL_TARGET", definition = "452 paired patients in final MIMIC target", variable = "GV_SD", poct = "gv_poct", lab = "gv_lab"),
  list(dd = candidate, cohort = "SENSITIVITY_ALL_PAIRED", definition = "453 final source-paired candidates", variable = "GV_SD", poct = "gv_poct", lab = "gv_lab"),
  list(dd = primary, cohort = "PRIMARY_FINAL_TARGET", definition = "452 paired patients in final MIMIC target", variable = "MEAN_GLUCOSE", poct = "mean_glu_poct", lab = "mean_glu_lab"),
  list(dd = candidate, cohort = "SENSITIVITY_ALL_PAIRED", definition = "453 final source-paired candidates", variable = "MEAN_GLUCOSE", poct = "mean_glu_poct", lab = "mean_glu_lab")
)) {
  rows[[length(rows) + 1L]] <- agreement_row(item$dd, item$cohort, item$definition, item$variable, item$poct, item$lab)
}
agreement <- do.call(rbind, rows)
agreement$lineage_version <- "JAHA_v5_final_lineage"
agreement$source_file <- phase2a_cfg_path(cfg, "mimic_samepatient_source")
phase2a_write(agreement, phase2a_output_path(cfg, "machine_readable", "phase2a_same_patient_agreement.csv"))

fig_dir <- phase2a_output_path(cfg, "figures", "diagnostic")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(fig_dir, "same_patient_gv_scatter_primary.png"), width = 1800, height = 1500, res = 220)
plot(primary$gv_poct, primary$gv_lab, pch = 19, col = rgb(0.12, 0.32, 0.62, 0.45),
     xlab = "POCT-derived GV SD", ylab = "Central-laboratory-derived GV SD",
     main = "Same-patient GV source comparison (primary final-target cohort)")
abline(0, 1, lty = 2, lwd = 2, col = "grey30")
abline(lm(gv_lab ~ gv_poct, data = primary), col = "firebrick", lwd = 2)
legend("topleft", legend = c("Identity line", "Descriptive linear fit"), lty = c(2, 1), col = c("grey30", "firebrick"), bty = "n")
dev.off()

gv <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "GV_SD", ]
diff <- primary$gv_poct - primary$gv_lab
pair_mean <- (primary$gv_poct + primary$gv_lab) / 2
png(file.path(fig_dir, "same_patient_gv_bland_altman_primary.png"), width = 1800, height = 1500, res = 220)
plot(pair_mean, diff, pch = 19, col = rgb(0.12, 0.32, 0.62, 0.45),
     xlab = "Pair mean GV SD", ylab = "POCT GV SD - laboratory GV SD",
     main = "Bland–Altman diagnostic (primary final-target cohort)")
abline(h = c(gv$mean_paired_difference, gv$loa_lower95, gv$loa_upper95),
       col = c("firebrick", "grey30", "grey30"), lty = c(1, 2, 2), lwd = c(2, 1, 1))
legend("topright", legend = c("Mean difference", "95% limits of agreement"), col = c("firebrick", "grey30"), lty = c(1, 2), bty = "n")
dev.off()

cat("PHASE2A_AGREEMENT_DONE\n")
