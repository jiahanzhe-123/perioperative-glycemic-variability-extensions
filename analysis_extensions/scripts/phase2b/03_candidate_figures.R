#!/usr/bin/env Rscript

phase2b_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2b_file, mustWork = TRUE)), "phase2b_common.R"))
cfg <- phase2b_load_config()
phase2a_open_log(cfg, "phase2b_03_candidate_figures.log")
on.exit(phase2a_close_log(), add = TRUE)

phase2b_require_gates(cfg)
agreement <- phase2b_read(phase2b_output_path(cfg, "machine_readable", "phase2a_same_patient_agreement.csv"))
log_scale <- phase2b_read(phase2b_output_path(cfg, "machine_readable", "phase2b_log_scale_agreement.csv"))
bootstrap <- phase2b_read(phase2b_output_path(cfg, "machine_readable", "phase2a_same_patient_bootstrap_summary.csv"))
landscape <- phase2b_read(phase2b_output_path(cfg, "machine_readable", "phase2b_analytic_context_landscape.csv"))
sp <- phase2a_source_data(cfg)$sp
target_ids <- unique(phase2a_target_ids(cfg))
primary <- sp[sp$stay_id %in% target_ids, , drop = FALSE]

fig_dir <- phase2b_output_path(cfg, "figures", "candidate")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

primary_agreement <- agreement[agreement$cohort_id == "PRIMARY_FINAL_TARGET" & agreement$variable == "GV_SD", , drop = FALSE]
primary_log <- log_scale[log_scale$cohort_id == "PRIMARY_FINAL_TARGET", , drop = FALSE]
diff <- phase2b_num(primary$gv_poct) - phase2b_num(primary$gv_lab)
pair_mean <- (phase2b_num(primary$gv_poct) + phase2b_num(primary$gv_lab)) / 2
log_keep <- is.finite(phase2b_num(primary$gv_poct)) & is.finite(phase2b_num(primary$gv_lab)) &
  phase2b_num(primary$gv_poct) > 0 & phase2b_num(primary$gv_lab) > 0
log_poct <- log(phase2b_num(primary$gv_poct)[log_keep])
log_lab <- log(phase2b_num(primary$gv_lab)[log_keep])
log_ratio <- log_poct - log_lab
log_pair_mean <- (log_poct + log_lab) / 2

png(file.path(fig_dir, "phase2b_candidate_figure_A_source_dependence.png"), width = 2600, height = 2100, res = 220)
par(mfrow = c(2, 2), mar = c(5.5, 5.5, 4.0, 1.5), oma = c(0, 0, 2.5, 0))

plot(primary$gv_poct, primary$gv_lab, pch = 19, col = rgb(0.12, 0.32, 0.62, 0.42),
     xlab = "POCT-derived GV SD", ylab = "Central-laboratory-derived GV SD",
     main = "A. Same-patient source comparison")
abline(0, 1, lty = 2, lwd = 2, col = "grey30")
abline(lm(gv_lab ~ gv_poct, data = primary), col = "firebrick", lwd = 2)
legend("topleft", legend = c("Identity line", "Descriptive linear fit"), lty = c(2, 1),
       col = c("grey30", "firebrick"), bty = "n", cex = 0.9)
mtext("Source-defined routine-care GV; neither source is a reference standard", side = 3, line = 0.2, cex = 0.78)

plot(pair_mean, diff, pch = 19, col = rgb(0.12, 0.32, 0.62, 0.42),
     xlab = "Pair mean GV SD", ylab = "POCT GV SD - laboratory GV SD",
     main = "B. Original-scale Bland–Altman")
abline(h = c(primary_agreement$mean_paired_difference, primary_agreement$loa_lower95, primary_agreement$loa_upper95),
       col = c("firebrick", "grey30", "grey30"), lty = c(1, 2, 2), lwd = c(2, 1, 1))
legend("topright", legend = c("Mean difference", "95% limits"), col = c("firebrick", "grey30"),
       lty = c(1, 2), bty = "n", cex = 0.9)

plot(log_pair_mean, log_ratio, pch = 19, col = rgb(0.12, 0.32, 0.62, 0.42),
     xlab = "Mean(log POCT GV, log laboratory GV)", ylab = "log(POCT GV / laboratory GV)",
     main = "C. Multiplicative-scale sensitivity (positive pairs only)")
abline(h = c(primary_log$mean_log_ratio, primary_log$log_loa_lower95, primary_log$log_loa_upper95),
       col = c("firebrick", "grey30", "grey30"), lty = c(1, 2, 2), lwd = c(2, 1, 1))
legend("topright", legend = c("Mean log ratio", "95% limits"), col = c("firebrick", "grey30"),
       lty = c(1, 2), bty = "n", cex = 0.9)

src <- phase2b_read(phase2a_cfg_path(cfg, "mimic_source_results"))
coef_ids <- c("SRC_samepatient_poct_30d", "SRC_samepatient_lab_30d")
coef <- src[match(coef_ids, src$model_id), , drop = FALSE]
ratio_est <- exp(phase2b_num(bootstrap$observed_delta_beta))
ratio_lo <- exp(phase2b_num(bootstrap$percentile_ci_delta_beta_lower))
ratio_hi <- exp(phase2b_num(bootstrap$percentile_ci_delta_beta_upper))
ests <- c(phase2b_num(coef$HR_per10), ratio_est)
los <- c(phase2b_num(coef$lo), ratio_lo)
his <- c(phase2b_num(coef$hi), ratio_hi)
labels <- c("POCT GV\nHR", "Laboratory GV\nHR", "Paired\nexp(delta-beta)")
ylim <- range(c(los, his, 1), finite = TRUE)
ylim <- c(max(0.45, ylim[1] * 0.8), min(1.35, ylim[2] * 1.15))
plot(seq_along(ests), ests, type = "n", log = "y", xaxt = "n", ylim = ylim,
     xlab = "Identical 409-patient / 49-event comparison", ylab = "Effect or paired ratio (log scale)",
     main = "D. Same-patient mortality coefficients")
axis(1, at = seq_along(ests), labels = labels, cex.axis = 0.78)
abline(h = 1, lty = 2, col = "grey30")
segments(seq_along(ests), los, seq_along(ests), his, lwd = 2, col = c("steelblue4", "steelblue4", "firebrick"))
points(seq_along(ests), ests, pch = 19, cex = 1.15, col = c("steelblue4", "steelblue4", "firebrick"))
legend("bottomleft", legend = c("Source-specific HR", "Paired descriptive ratio"), pch = 19,
       col = c("steelblue4", "firebrick"), bty = "n", cex = 0.82)
mtext("First two estimates use identical patients/events; no causal interaction is implied", side = 3, line = 0.2, cex = 0.72)

mtext("CANDIDATE FIGURE A — Same-patient measurement-source dependence", outer = TRUE, cex = 1.25, font = 2)
dev.off()

landscape_panel <- function(dd, title, effect_label, same_pair_legend = FALSE) {
  dd <- dd[order(dd$figure_order), , drop = FALSE]
  y <- seq_len(nrow(dd))
  y_labels <- if (same_pair_legend) {
    c("priority series", "POCT-only", "central-lab-only", "blood-gas-only", "common-source", "same-patient POCT", "same-patient laboratory")
  } else {
    dd$adjustment_model
  }
  xvals <- c(dd$lower95, dd$upper95, 1)
  xlim <- range(xvals[is.finite(xvals)], finite = TRUE)
  xlim <- c(max(0.45, xlim[1] * 0.88), min(1.45, xlim[2] * 1.10))
  plot(NA, NA, log = "x", xlim = xlim, ylim = c(nrow(dd) + 0.6, 0.4),
       yaxt = "n", xlab = paste0(effect_label, " (log scale; reference = 1)"),
       ylab = "Context order", main = title)
  axis(2, at = y, labels = y_labels, las = 1, cex.axis = 0.75)
  abline(v = 1, lty = 2, col = "grey30")
  cols <- ifelse(dd$same_patient_pair, "firebrick", "steelblue4")
  pchs <- ifelse(dd$same_patient_pair, 17, 19)
  segments(dd$lower95, y, dd$upper95, y, col = cols, lwd = 2)
  points(dd$estimate, y, col = cols, pch = pchs, cex = 1.0)
  text(dd$estimate, y, labels = paste0("  N=", dd$N, "/", dd$events), pos = 4, cex = 0.58, col = "grey25", offset = 0.25)
  if (same_pair_legend) legend("bottomright", legend = c("same-patient pair", "other source context"),
                                pch = c(17, 19), col = c("firebrick", "steelblue4"), bty = "n", cex = 0.72)
}

png(file.path(fig_dir, "phase2b_candidate_figure_B_analytic_context_landscape.png"), width = 3000, height = 2300, res = 220)
par(mfrow = c(2, 2), mar = c(5.5, 8.5, 4.2, 1.5), oma = c(0, 0, 2.5, 0))
landscape_panel(landscape[landscape$figure_panel == "MIMIC_adjustment", , drop = FALSE],
                "A. MIMIC adjustment context", "HR")
landscape_panel(landscape[landscape$figure_panel == "MIMIC_source", , drop = FALSE],
                "B. MIMIC source context", "HR", same_pair_legend = TRUE)
landscape_panel(landscape[landscape$figure_panel == "INSPIRE_timing_adjustment", , drop = FALSE],
                "C. INSPIRE timing/adjustment context", "HR")
landscape_panel(landscape[landscape$figure_panel == "eICU_adjustment", , drop = FALSE],
                "D. eICU adjustment context", "RR")
mtext("CANDIDATE FIGURE B — Analytic-context landscape (context families; not an exhaustive specification curve)",
      outer = TRUE, cex = 1.18, font = 2)
dev.off()

cat("PHASE2B_CANDIDATE_FIGURES_DONE\n")
