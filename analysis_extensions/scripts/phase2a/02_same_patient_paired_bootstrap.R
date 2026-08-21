#!/usr/bin/env Rscript

phase2a_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(phase2a_file, mustWork = TRUE)), "phase2a_common.R"))
suppressMessages({ library(survival); library(rms); library(jsonlite) })
cfg <- phase2a_load_config()
phase2a_open_log(cfg, "phase2a_02_same_patient_paired_bootstrap.log")
on.exit(phase2a_close_log(), add = TRUE)

if (!phase2a_phase1_gate(cfg)) stop("PHASE1_6_GATE_FAIL")
data <- phase2a_source_data(cfg)
base <- data$base
feat <- data$feat
sp <- data$sp
names(feat)[names(feat) != "stay_id"] <- paste0("ps_", names(feat)[names(feat) != "stay_id"])
d <- merge(base, feat, by = "stay_id")
d$mean_glu <- d$ps_mean_glucose
d$glucose_count <- d$ps_glucose_count
d$gender <- factor(d$gender)
d$procedure_cat6 <- factor(d$procedure_cat6,
  levels = c("isolated CABG", "isolated open valve", "combined CABG + open valve",
             "open aortic surgery (+/- other)", "transplant/VAD", "congenital/other open cardiac"))

K <- fromJSON(phase2a_cfg_path(cfg, "analysis_record_root") |> file.path("results", "standardization_constants.json"))
COVS <- c("age_at_admission", "gender", "bmi", "diabetes", "charlson_without_diabetes",
           "procedure_cat6", "lactate_postop_first", "creat_postop_first", "sofa_24h")
covs_fml <- paste(COVS, collapse = " + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits = 10), collapse = ","), "))")
model_formula <- paste0("Surv(t_lm_30,event_lm_30) ~ gv10 + ", rcs_mean, " + ", covs_fml)

tgt <- d[d$landmark_eligible == TRUE, , drop = FALSE]
dd <- tgt[tgt$stay_id %in% sp$stay_id, , drop = FALSE]
dd <- dd[complete.cases(dd[, c("t_lm_30", "t_lm_365", "event_lm_30", "event_lm_365", COVS)]), , drop = FALSE]
sp2 <- sp[sp$stay_id %in% dd$stay_id, , drop = FALSE]
if (nrow(dd) != 409L || nrow(sp2) != 409L || !all(dd$stay_id == sp2$stay_id)) stop("SOURCE_MODEL_ROW_ALIGNMENT_FAIL")
if (sum(dd$event_lm_30) != 49L) stop("SOURCE_MODEL_EVENT_COUNT_FAIL")

fit_source <- function(dd_in, gv_vec, mean_vec) {
  dd_fit <- dd_in
  dd_fit$gv10 <- phase2a_num(gv_vec) / 10
  dd_fit$mean_glu <- phase2a_num(mean_vec)
  keep <- complete.cases(dd_fit[, c("t_lm_30", "event_lm_30", "gv10", "mean_glu", COVS)]) & dd_fit$t_lm_30 > 0
  dd_fit <- dd_fit[keep, , drop = FALSE]
  if (nrow(dd_fit) < 100L || sum(dd_fit$event_lm_30) < 10L) {
    return(list(success = FALSE, reason = paste0("insufficient N/events: ", nrow(dd_fit), "/", sum(dd_fit$event_lm_30))))
  }
  warning_text <- character()
  fit <- tryCatch(
    withCallingHandlers(
      coxph(as.formula(model_formula), data = dd_fit),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(success = FALSE, reason = conditionMessage(fit), warnings = warning_text))
  sm <- summary(fit)
  if (!"gv10" %in% rownames(sm$coefficients)) return(list(success = FALSE, reason = "gv10 coefficient absent", warnings = warning_text))
  beta <- unname(sm$coefficients["gv10", "coef"])
  hr <- unname(sm$coefficients["gv10", "exp(coef)"])
  if (!is.finite(beta) || !is.finite(hr)) return(list(success = FALSE, reason = "non-finite gv10 coefficient", warnings = warning_text))
  list(
    success = TRUE,
    beta = beta,
    hr = hr,
    lo = unname(sm$conf.int["gv10", "lower .95"]),
    hi = unname(sm$conf.int["gv10", "upper .95"]),
    p = unname(sm$coefficients["gv10", "Pr(>|z|)"]),
    N = nrow(dd_fit),
    events = sum(dd_fit$event_lm_30),
    warnings = warning_text
  )
}

observed_poct <- fit_source(dd, sp2$gv_poct, sp2$mean_glu_poct)
observed_lab <- fit_source(dd, sp2$gv_lab, sp2$mean_glu_lab)
if (!observed_poct$success || !observed_lab$success) stop("SOURCE_MODEL_REPRODUCTION_FAIL")

locked <- phase2a_read(phase2a_cfg_path(cfg, "mimic_source_results"))
locked_poct <- locked$HR_per10[match("SRC_samepatient_poct_30d", locked$model_id)]
locked_lab <- locked$HR_per10[match("SRC_samepatient_lab_30d", locked$model_id)]
reproduction_ok <- is.finite(locked_poct) && is.finite(locked_lab) &&
  abs(observed_poct$hr - locked_poct) <= 1e-6 && abs(observed_lab$hr - locked_lab) <= 1e-6
reproduction <- data.frame(
  source = c("POCT", "central_laboratory"),
  N = c(observed_poct$N, observed_lab$N),
  events = c(observed_poct$events, observed_lab$events),
  reproduced_HR_per10 = c(observed_poct$hr, observed_lab$hr),
  locked_HR_per10 = c(locked_poct, locked_lab),
  absolute_difference = c(abs(observed_poct$hr - locked_poct), abs(observed_lab$hr - locked_lab)),
  reproduction_status = if (reproduction_ok) "PASS" else "FAIL",
  model_formula = model_formula,
  authoritative_script = phase2a_cfg_path(cfg, "public_source_model_script"),
  stringsAsFactors = FALSE
)
phase2a_write(reproduction, phase2a_output_path(cfg, "machine_readable", "phase2a_source_model_reproduction.csv"))
if (!reproduction_ok) stop("SOURCE_MODEL_REPRODUCTION_FAIL")

set.seed(as.integer(cfg$random_seed))
B <- 2000L
rep_rows <- vector("list", B)
for (b in seq_len(B)) {
  idx <- sample.int(nrow(dd), nrow(dd), replace = TRUE)
  dd_b <- dd[idx, , drop = FALSE]
  sp_b <- sp2[idx, , drop = FALSE]
  fit_p <- fit_source(dd_b, sp_b$gv_poct, sp_b$mean_glu_poct)
  fit_l <- fit_source(dd_b, sp_b$gv_lab, sp_b$mean_glu_lab)
  success <- isTRUE(fit_p$success) && isTRUE(fit_l$success)
  reason <- if (success) "" else paste(c(if (!fit_p$success) paste0("POCT: ", fit_p$reason), if (!fit_l$success) paste0("LAB: ", fit_l$reason)), collapse = " | ")
  rep_rows[[b]] <- data.frame(
    replicate = b,
    seed = as.integer(cfg$random_seed),
    fit_status = if (success) "SUCCESS" else "FAILED",
    hr_poct = if (success) fit_p$hr else NA_real_,
    hr_lab = if (success) fit_l$hr else NA_real_,
    delta_beta = if (success) log(fit_p$hr) - log(fit_l$hr) else NA_real_,
    coefficient_ratio = if (success) fit_p$hr / fit_l$hr else NA_real_,
    warning_count_poct = length(if (is.null(fit_p$warnings)) character() else fit_p$warnings),
    warning_count_lab = length(if (is.null(fit_l$warnings)) character() else fit_l$warnings),
    failure_reason = reason,
    stringsAsFactors = FALSE
  )
  if (b %% 100L == 0L) cat("bootstrap replicate", b, "of", B, "\n")
}
replicates <- do.call(rbind, rep_rows)
phase2a_write(replicates, phase2a_output_path(cfg, "machine_readable", "phase2a_same_patient_bootstrap.csv"))

ok_delta <- replicates$fit_status == "SUCCESS" & is.finite(replicates$delta_beta)
delta <- replicates$delta_beta[ok_delta]
failed <- sum(!ok_delta)
failure_pct <- 100 * failed / B
summary_row <- data.frame(
  analysis = "paired source-definition coefficient-difference sensitivity analysis",
  N = nrow(dd), events = sum(dd$event_lm_30),
  requested_replicates = B, successful_replicates = length(delta), failed_replicates = failed,
  failed_percent = failure_pct,
  observed_hr_poct = observed_poct$hr, observed_hr_lab = observed_lab$hr,
  observed_delta_beta = log(observed_poct$hr) - log(observed_lab$hr),
  observed_coefficient_ratio = observed_poct$hr / observed_lab$hr,
  bootstrap_median_delta_beta = if (length(delta)) median(delta) else NA_real_,
  bootstrap_sd_delta_beta = if (length(delta) > 1L) sd(delta) else NA_real_,
  percentile_ci_delta_beta_lower = if (length(delta)) quantile(delta, .025, names = FALSE, type = 7) else NA_real_,
  percentile_ci_delta_beta_upper = if (length(delta)) quantile(delta, .975, names = FALSE, type = 7) else NA_real_,
  failure_threshold_percent = 5,
  bootstrap_status = if (failure_pct > 5) "BOOTSTRAP_INSTABILITY" else "COMPLETED",
  seed = as.integer(cfg$random_seed),
  interpretation_boundary = "paired source-definition coefficient difference; not a causal interaction test",
  stringsAsFactors = FALSE
)
phase2a_write(summary_row, phase2a_output_path(cfg, "machine_readable", "phase2a_same_patient_bootstrap_summary.csv"))
cat("SOURCE_MODEL_REPRODUCTION_PASS\n")
cat("BOOTSTRAP_STATUS", summary_row$bootstrap_status, "failure_pct", failure_pct, "\n")
