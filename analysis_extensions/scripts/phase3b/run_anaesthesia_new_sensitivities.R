#!/usr/bin/env Rscript

# Anaesthesia candidate only: locked, post hoc sensitivity analyses.
# This script writes only to the extension workspace. It does not modify the
# public repository, controlled inputs, historical branches, or manuscript.

options(stringsAsFactors = FALSE, scipen = 999)

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
if (!length(file_arg)) stop("Run this file with Rscript so the extension root can be resolved.")
script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(extension_root)

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
  library(survival)
  library(mice)
  library(glmnet)
  library(rms)
})

cfg <- yaml::read_yaml(file.path(extension_root, "config.yaml"))
out_dir <- file.path(extension_root, "outputs/phase3b/journal_variants/anaesthesia/new_sensitivities")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

LOCKED_SEED_MICE <- 20260827L
LOCKED_SEED_PARS_BOOT <- 20260928L
LOCKED_SEED_RIDGE_FOLDS <- 20260929L
LOCKED_SEED_RIDGE_BOOT <- 20260930L
BOOT_TARGET <- 2000L

input_paths <- c(
  analysis_base = cfg$mimic_analysis_base,
  features_priority = cfg$mimic_features_priority,
  samepatient_source = cfg$mimic_samepatient_source,
  standardization_constants = file.path(cfg$analysis_record_root, "results/standardization_constants.json"),
  locked_primary_results = cfg$mimic_mice_results,
  locked_source_results = cfg$mimic_source_results
)
expected_hashes <- c(
  analysis_base = "2163168d1c2cf07828dd10d3b79392ca94545fbf36e7fd68dcd9d2675eeefba4",
  features_priority = "649dcf92d9f2fa35bab2521b123d54b1390ec508c1724ab8911ba1d0e599d906",
  samepatient_source = "9ff80386a5145f4cdc664a1ffe67db0e43a2bda08c93351bd81d6a9623b8774c",
  standardization_constants = "f67cda9ad42193e07472d979f44984284545b70fa8932b02e2fc99bae4d9e27e",
  locked_primary_results = "41c6a64f707f40518b5c809e67b26b23e9b3ade65ddea82708d5c5c1ce1c4d0a",
  locked_source_results = "ebb489c3f31bc879aa24e9ce98b6b5da66232a0312700c9ba237153c87180909"
)

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
if (any(!file.exists(input_paths))) {
  missing <- names(input_paths)[!file.exists(input_paths)]
  stop("Missing locked input(s): ", paste(missing, collapse = ", "))
}
observed_hashes <- vapply(input_paths, sha256_file, character(1))
if (!identical(unname(observed_hashes), unname(expected_hashes[names(observed_hashes)]))) {
  bad <- names(observed_hashes)[observed_hashes != expected_hashes[names(observed_hashes)]]
  stop("Locked input hash mismatch: ", paste(bad, collapse = ", "))
}

boolify <- function(x) x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
base <- read.csv(input_paths[["analysis_base"]], stringsAsFactors = FALSE, check.names = FALSE)
feat <- read.csv(input_paths[["features_priority"]], stringsAsFactors = FALSE, check.names = FALSE)
sp <- read.csv(input_paths[["samepatient_source"]], stringsAsFactors = FALSE, check.names = FALSE)
base$landmark_eligible <- boolify(base$landmark_eligible)

constants <- jsonlite::fromJSON(input_paths[["standardization_constants"]])
mean_knots <- as.numeric(constants$mean_glu_knots)
if (length(mean_knots) != 4L) stop("Expected four locked mean-glucose spline knots.")

procedure_levels <- c(
  "isolated CABG", "isolated open valve", "combined CABG + open valve",
  "open aortic surgery (+/- other)", "transplant/VAD",
  "congenital/other open cardiac"
)
clinical_covariates <- c(
  "age_at_admission", "gender", "bmi", "diabetes",
  "charlson_without_diabetes", "procedure_cat6", "lactate_postop_first",
  "creat_postop_first", "sofa_24h"
)

merged <- merge(base, feat, by = "stay_id", sort = FALSE)
# `span_hours` is present in both final-lineage frames; the feature-frame
# value is the declared retained-series span used by the locked MICE design.
if ("span_hours.y" %in% names(merged)) merged$span_hours <- merged$span_hours.y
merged$procedure_cat6 <- factor(merged$procedure_cat6, levels = procedure_levels)
merged$gender <- factor(merged$gender)
target <- merged[
  merged$landmark_eligible & !is.na(merged$gv_sd) & merged$glucose_count >= 2,
  , drop = FALSE
]
target_ge3 <- target[target$glucose_count >= 3, , drop = FALSE]

if (nrow(target) != 10561L || sum(target$event_lm_30) != 296L) {
  stop("Primary target fingerprint mismatch: observed ", nrow(target), "/", sum(target$event_lm_30))
}
if (nrow(target_ge3) != 10398L || sum(target_ge3$event_lm_30) != 272L) {
  stop(">=3 measurement fingerprint mismatch: observed ", nrow(target_ge3), "/", sum(target_ge3$event_lm_30))
}

pair_base <- base[, c(
  "stay_id", "landmark_eligible", "t_lm_30", "event_lm_30", "t_lm_365",
  "event_lm_365", clinical_covariates
)]
paired <- merge(sp, pair_base, by = "stay_id", sort = FALSE)
paired <- paired[paired$landmark_eligible & paired$t_lm_30 > 0, , drop = FALSE]
paired_cc <- paired[complete.cases(paired[, c(
  "t_lm_30", "event_lm_30", "gv_poct", "mean_glu_poct", "gv_lab",
  "mean_glu_lab", clinical_covariates
)]), , drop = FALSE]
if (nrow(paired_cc) != 409L || sum(paired_cc$event_lm_30) != 49L) {
  stop("Paired complete-case fingerprint mismatch: observed ", nrow(paired_cc), "/", sum(paired_cc$event_lm_30))
}

fmt_formula <- function(...) as.formula(paste0(...))
source_label <- c(poct = "POCT", lab = "laboratory")

fit_parsimonious <- function(dd, source) {
  gv <- dd[[paste0("gv_", source)]] / 10
  mg <- dd[[paste0("mean_glu_", source)]] / 10
  dat <- data.frame(
    time = dd$t_lm_30, event = dd$event_lm_30, gv10 = gv, mean10 = mg,
    age = dd$age_at_admission, gender = factor(dd$gender),
    diabetes = as.numeric(dd$diabetes), procedure_cat6 = factor(dd$procedure_cat6, levels = procedure_levels)
  )
  fit <- survival::coxph(
    survival::Surv(time, event) ~ gv10 + mean10 + age + gender + diabetes + procedure_cat6,
    data = dat, ties = "efron", x = TRUE, model = TRUE
  )
  sm <- summary(fit)
  co <- sm$coefficients["gv10", , drop = FALSE]
  ci <- sm$conf.int["gv10", , drop = FALSE]
  list(
    fit = fit,
    beta = unname(co[, "coef"]),
    estimate = unname(co[, "exp(coef)"]),
    lower95 = unname(ci[, "lower .95"]),
    upper95 = unname(ci[, "upper .95"]),
    p_value = unname(co[, "Pr(>|z|)"])
  )
}

fit_parsimonious_boot <- function(dd, seed, target_successes) {
  set.seed(seed)
  n <- nrow(dd); successes <- 0L; attempts <- 0L; rows <- list()
  while (successes < target_successes && attempts < 4000L) {
    attempts <- attempts + 1L
    idx <- sample.int(n, n, replace = TRUE)
    one <- tryCatch(fit_parsimonious(dd[idx, , drop = FALSE], "poct"), error = function(e) NULL)
    two <- tryCatch(fit_parsimonious(dd[idx, , drop = FALSE], "lab"), error = function(e) NULL)
    if (is.null(one) || is.null(two) || !all(is.finite(c(one$beta, two$beta)))) next
    successes <- successes + 1L
    rows[[successes]] <- data.frame(
      analysis = "parsimonious_cox", replicate = successes,
      beta_poct = one$beta, beta_lab = two$beta,
      delta_beta = one$beta - two$beta,
      stringsAsFactors = FALSE
    )
  }
  if (successes < target_successes) stop("Parsimonious bootstrap did not reach ", target_successes, " successful replicates.")
  list(data = do.call(rbind, rows), attempts = attempts, successful = successes)
}

rcs_formula <- function(source) {
  gv_name <- paste0("gv_", source)
  mg_name <- paste0("mean_glu_", source)
  # The explicit fixed knots reproduce the locked Model B mean-glucose term.
  fmt_formula(
    "~ ", gv_name, "10 + rms::rcs(", mg_name, ", c(",
    paste(format(mean_knots, digits = 10), collapse = ","), ")) + ",
    paste(clinical_covariates, collapse = " + "), " - 1"
  )
}

make_ridge_matrix <- function(dd, source) {
  d <- dd
  d[[paste0("gv_", source, "10")]] <- d[[paste0("gv_", source)]] / 10
  d[[paste0("mean_glu_", source)]] <- as.numeric(d[[paste0("mean_glu_", source)]])
  d$gender <- factor(d$gender)
  d$procedure_cat6 <- factor(d$procedure_cat6, levels = procedure_levels)
  x <- model.matrix(rcs_formula(source), data = d)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  gv_col <- paste0("gv_", source, "10")
  if (!gv_col %in% colnames(x)) stop("Ridge design does not contain ", gv_col)
  list(x = x, gv_col = gv_col)
}

fit_ridge <- function(x, y, gv_col, lambda = NULL, foldid = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cv <- NULL
  if (is.null(lambda)) {
    cv <- glmnet::cv.glmnet(
      x = x, y = y, family = "cox", alpha = 0, nfolds = 10,
      foldid = foldid, standardize = TRUE, type.measure = "deviance",
      keep = FALSE
    )
    lambda <- cv$lambda.1se
  }
  fit <- glmnet::glmnet(
    x = x, y = y, family = "cox", alpha = 0, lambda = lambda,
    standardize = TRUE
  )
  beta <- as.numeric(as.matrix(stats::coef(fit))[gv_col, 1])
  list(fit = fit, cv = cv, lambda = lambda, beta = beta, estimate = exp(beta))
}

fit_ridge_boot <- function(x_poct, x_lab, y, gv_poct, gv_lab, lambda_poct, lambda_lab, seed, target_successes) {
  set.seed(seed)
  n <- nrow(x_poct); successes <- 0L; attempts <- 0L; rows <- list()
  while (successes < target_successes && attempts < 4000L) {
    attempts <- attempts + 1L
    idx <- sample.int(n, n, replace = TRUE)
    one <- tryCatch(fit_ridge(x_poct[idx, , drop = FALSE], y[idx], gv_poct, lambda_poct), error = function(e) NULL)
    two <- tryCatch(fit_ridge(x_lab[idx, , drop = FALSE], y[idx], gv_lab, lambda_lab), error = function(e) NULL)
    if (is.null(one) || is.null(two) || !all(is.finite(c(one$beta, two$beta)))) next
    successes <- successes + 1L
    rows[[successes]] <- data.frame(
      analysis = "ridge_cox", replicate = successes,
      beta_poct = one$beta, beta_lab = two$beta,
      delta_beta = one$beta - two$beta,
      stringsAsFactors = FALSE
    )
  }
  if (successes < target_successes) stop("Ridge bootstrap did not reach ", target_successes, " successful replicates.")
  list(data = do.call(rbind, rows), attempts = attempts, successful = successes)
}

percentile <- function(x) as.numeric(stats::quantile(x, probs = c(0.025, 0.975), names = FALSE, type = 6))
safe_p <- function(x) ifelse(is.finite(x), x, NA_real_)

# --------------------------- Analysis 1 ---------------------------
pars_fit <- lapply(c("poct", "lab"), function(s) fit_parsimonious(paired_cc, s))
names(pars_fit) <- c("poct", "lab")
pars_boot <- fit_parsimonious_boot(paired_cc, LOCKED_SEED_PARS_BOOT, BOOT_TARGET)
pars_delta_ci <- percentile(pars_boot$data$delta_beta)

# --------------------------- Analysis 2 ---------------------------
y <- survival::Surv(paired_cc$t_lm_30, paired_cc$event_lm_30)
set.seed(LOCKED_SEED_RIDGE_FOLDS)
foldid <- sample(rep(seq_len(10L), length.out = nrow(paired_cc)))
ridge_design <- lapply(c("poct", "lab"), function(s) make_ridge_matrix(paired_cc, s))
names(ridge_design) <- c("poct", "lab")
ridge_fit_poct <- fit_ridge(ridge_design$poct$x, y, ridge_design$poct$gv_col, foldid = foldid)
ridge_fit_lab <- fit_ridge(ridge_design$lab$x, y, ridge_design$lab$gv_col, foldid = foldid)
ridge_boot <- fit_ridge_boot(
  ridge_design$poct$x, ridge_design$lab$x, y,
  ridge_design$poct$gv_col, ridge_design$lab$gv_col,
  ridge_fit_poct$lambda, ridge_fit_lab$lambda,
  LOCKED_SEED_RIDGE_BOOT, BOOT_TARGET
)
ridge_poct_ci <- percentile(ridge_boot$data$beta_poct)
ridge_lab_ci <- percentile(ridge_boot$data$beta_lab)
ridge_delta_ci <- percentile(ridge_boot$data$delta_beta)

# --------------------------- Analysis 3 ---------------------------
mi_data <- target_ge3[, c(
  "stay_id", "t_lm_30", "t_lm_365", "event_lm_30", "event_lm_365",
  "gv_sd", "mean_glucose", "glucose_count", "span_hours",
  "frac_central_lab", "frac_blood_gas", "frac_poct", "frac_icu_charted",
  clinical_covariates,
  "albumin_adm_first", "hgb_postop_first", "wbc_postop_first", "platelets_postop_first"
)]
names(mi_data)[names(mi_data) == "gv_sd"] <- "gv"
names(mi_data)[names(mi_data) == "mean_glucose"] <- "mean_glu"
mi_data$diabetes <- as.numeric(mi_data$diabetes)
mi_data$H365 <- mice::nelsonaalen(mi_data, "t_lm_365", "event_lm_365")
pred <- mice::make.predictorMatrix(mi_data); pred[,] <- 0L
meth <- mice::make.method(mi_data); meth[] <- ""
cont_imp <- intersect(c(
  "bmi", "lactate_postop_first", "creat_postop_first", "sofa_24h",
  "albumin_adm_first", "hgb_postop_first", "wbc_postop_first", "platelets_postop_first"
), names(mi_data))
cont_imp <- cont_imp[colSums(is.na(mi_data[, cont_imp, drop = FALSE])) > 0]
impute_preds <- setdiff(names(mi_data), c("stay_id", "t_lm_30", "frac_blood_gas"))
for (v in cont_imp) {
  meth[v] <- "pmm"
  pred[v, impute_preds] <- 1L
}
mi_fit <- mice::mice(
  mi_data, m = 50, maxit = 20, method = meth, predictorMatrix = pred,
  seed = LOCKED_SEED_MICE, printFlag = FALSE
)
mi_betas <- numeric(0); mi_vars <- numeric(0)
ge3_rcs_mean <- paste0(
  "rms::rcs(mean_glu, c(", paste(format(mean_knots, digits = 10), collapse = ","), "))"
)
ge3_formula <- as.formula(paste(
  "survival::Surv(t_lm_30, event_lm_30) ~ gv10 +", ge3_rcs_mean,
  "+ age_at_admission + gender + bmi + diabetes + charlson_without_diabetes +",
  "procedure_cat6 + lactate_postop_first + creat_postop_first + sofa_24h"
))
for (i in seq_len(50L)) {
  dd <- mice::complete(mi_fit, i)
  dd$gv10 <- dd$gv / 10
  dd$log_count <- log(dd$glucose_count)
  dd$procedure_cat6 <- factor(dd$procedure_cat6, levels = procedure_levels)
  fit <- survival::coxph(
    ge3_formula,
    data = dd, ties = "efron"
  )
  sm <- summary(fit)
  mi_betas <- c(mi_betas, sm$coefficients["gv10", "coef"])
  mi_vars <- c(mi_vars, sm$coefficients["gv10", "se(coef)"]^2)
}
rubin <- function(betas, vars) {
  m <- length(betas); qbar <- mean(betas); ubar <- mean(vars); b <- stats::var(betas)
  total <- ubar + (1 + 1 / m) * b
  if (b == 0) df <- Inf else df <- (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
  list(beta = qbar, se = sqrt(total), df = df)
}
mi_pool <- rubin(mi_betas, mi_vars)
mi_t <- mi_pool$beta / mi_pool$se
mi_p <- 2 * stats::pt(abs(mi_t), df = mi_pool$df, lower.tail = FALSE)
mi_crit <- stats::qt(0.975, df = mi_pool$df)
mi_hr <- exp(mi_pool$beta)
mi_lo <- exp(mi_pool$beta - mi_crit * mi_pool$se)
mi_hi <- exp(mi_pool$beta + mi_crit * mi_pool$se)

result_rows <- list()
add_result <- function(...) result_rows[[length(result_rows) + 1L]] <<- data.frame(...)
add_result(
  analysis_id = "ANAESTHESIA_SENS1_PARSIMONIOUS_POCT_30D", analysis_family = "parsimonious same-patient Cox",
  source = "POCT", outcome = "30-day all-cause mortality", model = "linear GV + linear mean glucose + age + sex + diabetes + procedure category",
  N = nrow(paired_cc), events = sum(paired_cc$event_lm_30), effect_type = "HR per 10 mg.dl^-1 GV",
  estimate = pars_fit$poct$estimate, lower95 = pars_fit$poct$lower95, upper95 = pars_fit$poct$upper95,
  p_value = pars_fit$poct$p_value, lambda = NA_real_, bootstrap_lower95 = NA_real_, bootstrap_upper95 = NA_real_,
  fit_status = "PASS", estimate_origin = "new_sensitivity_results.csv:parsimonious_cox", stringsAsFactors = FALSE
)
add_result(
  analysis_id = "ANAESTHESIA_SENS1_PARSIMONIOUS_LAB_30D", analysis_family = "parsimonious same-patient Cox",
  source = "laboratory", outcome = "30-day all-cause mortality", model = "linear GV + linear mean glucose + age + sex + diabetes + procedure category",
  N = nrow(paired_cc), events = sum(paired_cc$event_lm_30), effect_type = "HR per 10 mg.dl^-1 GV",
  estimate = pars_fit$lab$estimate, lower95 = pars_fit$lab$lower95, upper95 = pars_fit$lab$upper95,
  p_value = pars_fit$lab$p_value, lambda = NA_real_, bootstrap_lower95 = NA_real_, bootstrap_upper95 = NA_real_,
  fit_status = "PASS", estimate_origin = "new_sensitivity_results.csv:parsimonious_cox", stringsAsFactors = FALSE
)
add_result(
  analysis_id = "ANAESTHESIA_SENS1_PARSIMONIOUS_DELTA_BETA_30D", analysis_family = "parsimonious same-patient Cox",
  source = "POCT minus laboratory", outcome = "30-day all-cause mortality", model = "paired coefficient contrast",
  N = nrow(paired_cc), events = sum(paired_cc$event_lm_30), effect_type = "paired log-HR difference",
  estimate = pars_fit$poct$beta - pars_fit$lab$beta, lower95 = pars_delta_ci[1], upper95 = pars_delta_ci[2],
  p_value = NA_real_, lambda = NA_real_, bootstrap_lower95 = pars_delta_ci[1], bootstrap_upper95 = pars_delta_ci[2],
  fit_status = "PASS", estimate_origin = "new_sensitivity_paired_bootstrap.csv:parsimonious_cox", stringsAsFactors = FALSE
)
add_result(
  analysis_id = "ANAESTHESIA_SENS2_RIDGE_POCT_30D", analysis_family = "ridge same-patient Cox",
  source = "POCT", outcome = "30-day all-cause mortality", model = "locked Model B clinical frame; ridge alpha=0; lambda.1se",
  N = nrow(paired_cc), events = sum(paired_cc$event_lm_30), effect_type = "ridge HR per 10 mg.dl^-1 GV",
  estimate = ridge_fit_poct$estimate, lower95 = exp(ridge_poct_ci[1]), upper95 = exp(ridge_poct_ci[2]),
  p_value = NA_real_, lambda = ridge_fit_poct$lambda, bootstrap_lower95 = exp(ridge_poct_ci[1]), bootstrap_upper95 = exp(ridge_poct_ci[2]),
  fit_status = "PASS", estimate_origin = "new_sensitivity_paired_bootstrap.csv:ridge_cox", stringsAsFactors = FALSE
)
add_result(
  analysis_id = "ANAESTHESIA_SENS2_RIDGE_LAB_30D", analysis_family = "ridge same-patient Cox",
  source = "laboratory", outcome = "30-day all-cause mortality", model = "locked Model B clinical frame; ridge alpha=0; lambda.1se",
  N = nrow(paired_cc), events = sum(paired_cc$event_lm_30), effect_type = "ridge HR per 10 mg.dl^-1 GV",
  estimate = ridge_fit_lab$estimate, lower95 = exp(ridge_lab_ci[1]), upper95 = exp(ridge_lab_ci[2]),
  p_value = NA_real_, lambda = ridge_fit_lab$lambda, bootstrap_lower95 = exp(ridge_lab_ci[1]), bootstrap_upper95 = exp(ridge_lab_ci[2]),
  fit_status = "PASS", estimate_origin = "new_sensitivity_paired_bootstrap.csv:ridge_cox", stringsAsFactors = FALSE
)
add_result(
  analysis_id = "ANAESTHESIA_SENS2_RIDGE_DELTA_BETA_30D", analysis_family = "ridge same-patient Cox",
  source = "POCT minus laboratory", outcome = "30-day all-cause mortality", model = "paired ridge coefficient contrast; full-cohort lambdas fixed",
  N = nrow(paired_cc), events = sum(paired_cc$event_lm_30), effect_type = "paired log-HR difference",
  estimate = ridge_fit_poct$beta - ridge_fit_lab$beta, lower95 = ridge_delta_ci[1], upper95 = ridge_delta_ci[2],
  p_value = NA_real_, lambda = NA_real_, bootstrap_lower95 = ridge_delta_ci[1], bootstrap_upper95 = ridge_delta_ci[2],
  fit_status = "PASS", estimate_origin = "new_sensitivity_paired_bootstrap.csv:ridge_cox", stringsAsFactors = FALSE
)
add_result(
  analysis_id = "ANAESTHESIA_SENS3_MIMIC_GE3_MODEL_B_30D", analysis_family = "MIMIC >=3 retained measurements",
  source = "priority series", outcome = "30-day all-cause mortality", model = "MICE Model B; m=50; pmm; Rubin pooling",
  N = nrow(target_ge3), events = sum(target_ge3$event_lm_30), effect_type = "HR per 0.555 mmol.l^-1 (10 mg.dl^-1) GV",
  estimate = mi_hr, lower95 = mi_lo, upper95 = mi_hi, p_value = mi_p, lambda = NA_real_,
  bootstrap_lower95 = NA_real_, bootstrap_upper95 = NA_real_, fit_status = "PASS",
  estimate_origin = "new_sensitivity_mimic_ge3_mice.csv:Model_B_30d", stringsAsFactors = FALSE
)
results <- do.call(rbind, result_rows)
write.csv(results, file.path(out_dir, "new_sensitivity_results.csv"), row.names = FALSE, na = "")

boot <- rbind(pars_boot$data, ridge_boot$data)
write.csv(boot, file.path(out_dir, "new_sensitivity_paired_bootstrap.csv"), row.names = FALSE, na = "")

ge3 <- data.frame(
  analysis_id = "ANAESTHESIA_SENS3_MIMIC_GE3_MODEL_B_30D",
  cohort_definition = "Final priority-series landmark cohort with >=3 retained glucose measurements",
  N = nrow(target_ge3), events = sum(target_ge3$event_lm_30),
  HR_per10 = mi_hr, lower95 = mi_lo, upper95 = mi_hi, p_value = mi_p,
  beta = mi_pool$beta, se = mi_pool$se, df = mi_pool$df,
  m = 50L, maxit = 20L, imputation_method = "pmm", seed = LOCKED_SEED_MICE,
  fit_status = "PASS", stringsAsFactors = FALSE
)
write.csv(ge3, file.path(out_dir, "new_sensitivity_mimic_ge3_mice.csv"), row.names = FALSE, na = "")

fmt3 <- function(x) sprintf("%.3f", as.numeric(x))
fmt_p <- function(x) {
  if (is.na(x)) return("not estimated")
  if (x < 0.001) return("<0.001")
  sprintf("%.3f", x)
}
get_result <- function(id) results[results$analysis_id == id, , drop = FALSE]
display_row <- function(id, label, interpretation) {
  r <- get_result(id)
  data.frame(
    `Sensitivity analysis` = label,
    `N/events` = paste0(r$N, "/", r$events),
    `Estimate (95% CI)` = paste0(fmt3(r$estimate), " (", fmt3(r$lower95), " to ", fmt3(r$upper95), ")"),
    `p value` = fmt_p(r$p_value),
    Interpretation = interpretation,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
display <- rbind(
  display_row("ANAESTHESIA_SENS1_PARSIMONIOUS_POCT_30D", "Parsimonious Cox: POCT-derived GV HR", "Reduced covariate complexity"),
  display_row("ANAESTHESIA_SENS1_PARSIMONIOUS_LAB_30D", "Parsimonious Cox: laboratory-derived GV HR", "Reduced covariate complexity"),
  display_row("ANAESTHESIA_SENS1_PARSIMONIOUS_DELTA_BETA_30D", "Parsimonious Cox: paired Δβ (POCT minus laboratory)", "Bootstrap percentile interval; not an HR"),
  display_row("ANAESTHESIA_SENS2_RIDGE_POCT_30D", "Ridge Cox: POCT-derived GV HR", "λ.1se penalised coefficient; bootstrap interval"),
  display_row("ANAESTHESIA_SENS2_RIDGE_LAB_30D", "Ridge Cox: laboratory-derived GV HR", "λ.1se penalised coefficient; bootstrap interval"),
  display_row("ANAESTHESIA_SENS2_RIDGE_DELTA_BETA_30D", "Ridge Cox: paired Δβ (POCT minus laboratory)", "Bootstrap percentile interval; not an HR"),
  display_row("ANAESTHESIA_SENS3_MIMIC_GE3_MODEL_B_30D", "MIMIC Model B: at least 3 retained measurements", "MICE m=50; Rubin pooling")
)
write.csv(display, file.path(out_dir, "new_sensitivity_table.csv"), row.names = FALSE, na = "")

manifest <- list(
  run_status = "COMPLETED",
  run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  extension_root = extension_root,
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  input_paths = unname(as.list(input_paths)),
  input_sha256 = unname(as.list(observed_hashes)),
  seeds = list(mice_ge3 = LOCKED_SEED_MICE, parsimonious_bootstrap = LOCKED_SEED_PARS_BOOT,
               ridge_cv_folds = LOCKED_SEED_RIDGE_FOLDS, ridge_bootstrap = LOCKED_SEED_RIDGE_BOOT),
  bootstrap = list(requested_successful = BOOT_TARGET,
                   parsimonious_attempts = pars_boot$attempts, parsimonious_successful = pars_boot$successful,
                   ridge_attempts = ridge_boot$attempts, ridge_successful = ridge_boot$successful),
  fingerprints = list(primary_target_N = nrow(target), primary_target_events = sum(target$event_lm_30),
                      paired_complete_case_N = nrow(paired_cc), paired_complete_case_events = sum(paired_cc$event_lm_30),
                      ge3_target_N = nrow(target_ge3), ge3_target_events = sum(target_ge3$event_lm_30)),
  output_files = c("new_sensitivity_results.csv", "new_sensitivity_paired_bootstrap.csv",
                   "new_sensitivity_mimic_ge3_mice.csv", "new_sensitivity_table.csv",
                   "new_sensitivity_qc.md"),
  locked_reference = list(primary_model_b_hr_per10 = 0.979380602603274,
                          samepatient_poct_hr_per10 = 0.758857245048576,
                          samepatient_lab_hr_per10 = 0.970897006154936),
  session = paste(capture.output(sessionInfo()), collapse = "\n")
)
manifest_path <- file.path(out_dir, "new_sensitivity_run_manifest.json")
write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"), manifest_path)

qc <- c(
  "# Anaesthesia new sensitivity-analysis QC", "",
  paste0("Run status: **", manifest$run_status, "**"), "",
  "## Passed checks", "",
  paste0("- Locked input SHA-256 values matched before fitting."),
  paste0("- Primary target fingerprint verified: N=", nrow(target), ", events=", sum(target$event_lm_30), "."),
  paste0("- Same-patient complete-case fingerprint verified: N=", nrow(paired_cc), ", events=", sum(paired_cc$event_lm_30), "."),
  paste0("- At-least-three-measurement fingerprint verified: N=", nrow(target_ge3), ", events=", sum(target_ge3$event_lm_30), "."),
  paste0("- Parsimonious Cox bootstrap successful replicates: ", pars_boot$successful, "/", BOOT_TARGET, "."),
  paste0("- Ridge Cox bootstrap successful replicates: ", ridge_boot$successful, "/", BOOT_TARGET, "."),
  "- MICE Model B fit completed for all 50 imputed datasets and pooled with Rubin's rules.",
  "- No manuscript or authoritative input file was written by the runner.", "",
  "## Interpretation boundary", "",
  "These are post hoc sensitivity analyses. Their purpose is to assess robustness to covariate complexity and a minimum measurement-count boundary; they do not convert the primary null result into a confirmatory positive result and do not establish a causal source interaction."
)
writeLines(qc, file.path(out_dir, "new_sensitivity_qc.md"))

cat("ANAESTHESIA_NEW_SENSITIVITIES_COMPLETED\n")
print(results[, c("analysis_id", "N", "events", "estimate", "lower95", "upper95", "p_value", "lambda", "fit_status")], row.names = FALSE)
