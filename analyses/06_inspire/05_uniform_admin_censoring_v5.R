# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 05_uniform_admin_censoring_v5.R — INSPIRE 30日全因死亡的统一行政删失(analysis of record)。
#
# 迁移自 v5 复现包 code/current_v5/18_uniform_admin_censoring_rerun_v5.R(逻辑不变,
# 路径改由 config/paths.yml 驱动)。这是 INSPIRE 模块的 analysis-of-record 脚本:
# 稿件报告值(I2 30d HR 0.905;48h HR 1.104)由本脚本产生。
#
# 时间规则:死亡者用复合死亡时间;所有非事件在统一的术后第30日行政终点删失
# (研究团队对死亡登记覆盖至30日的溯源声明);绝不使用出院删失。
# 365日 Cox 有意不重跑:无同等文件证的统一登记覆盖终点(见 365d_withdrawal_v5.csv)。

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
# rm(list=ls()) 会清除页首 source 的配置,须重新加载:
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])

SEED <- 20260726L
set.seed(SEED)
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({
  library(survival)
  library(jsonlite)
  library(rms)
})

SOURCE_ROOT <- PGV("inspire_work")
OUT_ROOT <- PGV("inspire_record_work")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT_ROOT, "18_uniform_admin_censoring_v5.log"), split = TRUE)

OUTCOME_VERSION <- "30d_reconciled_composite_v3 + uniform_admin_registry_censoring_v5"
TIME_RULE <- "Composite death time for deaths; common postoperative day-30 administrative censoring for all non-events (study-team attestation of linked-mortality coverage through day 30); never discharge censoring"

bool_int <- function(x) as.integer(x %in% c(TRUE, "t", "True", "TRUE", "true", 1, "1"))
rbind_fill <- function(x) {
  nm <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(z) {
    for (j in setdiff(nm, names(z))) z[[j]] <- NA
    z[, nm, drop = FALSE]
  }))
}

# Apply the fixed, external end of mortality-registry coverage.  `day30_cutoff`
# is used only because coverage through this point was confirmed outside the raw
# export; the raw extract does not contain a per-patient coverage-end field.
admin_day30 <- function(dd, landmark_col) {
  lm <- dd[[landmark_col]]
  death <- dd$death_time_composite
  cutoff <- dd$day30_cutoff
  dd$event_30d <- as.integer(!is.na(death) & death > lm & death <= cutoff)
  dd$followup_end_30 <- pmin(ifelse(is.na(death), Inf, death), cutoff)
  dd$t30 <- (dd$followup_end_30 - lm) / 1440.0
  stopifnot(all(is.finite(dd$t30)), all(dd$t30 > 0))
  stopifnot(all(abs(dd$t30[dd$event_30d == 1L] -
                    (death[dd$event_30d == 1L] - lm[dd$event_30d == 1L]) / 1440.0) < 1e-10))
  stopifnot(all(abs(dd$t30[dd$event_30d == 0L] -
                    (cutoff[dd$event_30d == 0L] - lm[dd$event_30d == 0L]) / 1440.0) < 1e-10))
  dd
}

r30 <- read.csv(file.path(SOURCE_ROOT, "data", "outcome_30d_reconciled.csv"))
base <- read.csv(file.path(SOURCE_ROOT, "data", "inspire_base.csv"))
comorb <- read.csv(file.path(SOURCE_ROOT, "data", "comorbidity.csv"))
anchor_features <- read.csv(file.path(SOURCE_ROOT, "data", "glucose_features_anchor.csv"))
K <- fromJSON(file.path(SOURCE_ROOT, "results", "standardization_constants.json"))

# ---- 24-hour landmark frame ----
d <- merge(r30[, c("subject_id", "op_id", "landmark_24h", "day30_cutoff",
                   "death_time_composite", "event_30d_reconciled")],
           base, by = c("subject_id", "op_id"))
d <- merge(d, comorb, by = "subject_id", all.x = TRUE)
d <- d[is.na(d$death_time_composite) | d$death_time_composite > d$landmark_24h, ]
d <- admin_day30(d, "landmark_24h")
d$gv <- d$gv_sd
d$mean_glu <- d$mean_glucose
d$n_gv <- d$n_glucose_0_24h
d$gv10 <- d$gv / 10
d$log_count <- log(d$n_gv)
d$bmi <- d$weight / (d$height / 100)^2
d$diabetes <- bool_int(d$diabetes)
d$sex <- factor(d$sex)
d$asa_f <- factor(d$asa)
d$emop_f <- bool_int(d$emop)
d$surgery_group <- factor(d$surgery_group)

stopifnot(nrow(d) == 1353L, sum(d$event_30d) == 27L)
stopifnot(identical(d$event_30d, bool_int(d$event_30d_reconciled)))
rcs_mean <- paste0("rms::rcs(mean_glu, c(",
                   paste(format(K$mean_glu_knots, digits = 10), collapse = ","), "))")
kSD <- K$gv_sd / 10
mask30 <- c("t30", "event_30d", "gv", "gv10", "mean_glu", "age", "log_count", "span_hours")
stopifnot(all(colSums(is.na(d[, mask30, drop = FALSE])) == 0L))
cc30 <- d[complete.cases(d[, mask30, drop = FALSE]), ]
stopifnot(nrow(cc30) == 1353L, sum(cc30$event_30d) == 27L)

rep_fit <- function(fit, model_id, model, dd, note) {
  s <- summary(fit)
  data.frame(
    model_id = model_id, frame = "INSPIRE_OPEND_24H_LANDMARK_V2",
    outcome_version = OUTCOME_VERSION, time_rule = TIME_RULE, model = model,
    N = nrow(dd), events = sum(dd$event_30d),
    HR_per10 = s$coefficients["gv10", "exp(coef)"],
    lo = s$conf.int["gv10", "lower .95"], hi = s$conf.int["gv10", "upper .95"],
    P = s$coefficients["gv10", "Pr(>|z|)"],
    HR_perSD = exp(s$coefficients["gv10", "coef"] * kSD),
    ph_global = tryCatch(cox.zph(fit)$table["GLOBAL", "p"], error = function(e) NA_real_),
    note = note, stringsAsFactors = FALSE
  )
}

f_i1 <- coxph(Surv(t30, event_30d) ~ gv10 + age, data = cc30, x = TRUE, y = TRUE)
f_i2 <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean, " + age")),
              data = cc30, x = TRUE, y = TRUE)
f_i3 <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean,
                                " + age + log_count + span_hours")),
              data = cc30, x = TRUE, y = TRUE)
res30 <- rbind(
  rep_fit(f_i1, "ADMINV5_I1_30d", "Model I1", cc30,
          "Formula-variable-only frame; uniform day-30 administrative censoring"),
  rep_fit(f_i2, "ADMINV5_I2_30d", "Model I2 (designated INSPIRE extension)", cc30,
          "Formula-variable-only frame; uniform day-30 administrative censoring"),
  rep_fit(f_i3, "ADMINV5_I3_30d", "Model I3", cc30,
          "Formula-variable-only frame; uniform day-30 administrative censoring")
)
write.csv(res30, file.path(OUT_ROOT, "30d_primary_results_v5.csv"), row.names = FALSE)
print(res30[, c("model_id", "N", "events", "HR_per10", "lo", "hi", "P", "ph_global")])

# ---- 30-day standardized risk contrast, I2 ----
risk_contrast <- function(dd) {
  fit <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean, " + age")), data = dd)
  bh <- basehaz(fit, centered = FALSE)
  H <- max(bh$hazard[bh$time <= max(dd$t30)])
  q25 <- dd; q25$gv <- K$gv_q25; q25$gv10 <- q25$gv / 10
  q75 <- dd; q75$gv <- K$gv_q75; q75$gv10 <- q75$gv / 10
  p25 <- mean(1 - exp(-H * exp(predict(fit, q25, type = "lp", reference = "zero"))))
  p75 <- mean(1 - exp(-H * exp(predict(fit, q75, type = "lp", reference = "zero"))))
  c(rd = p75 - p25, rr = p75 / p25, risk_q25 = p25, risk_q75 = p75)
}

point <- risk_contrast(cc30)
B <- 1000L
boot <- matrix(NA_real_, nrow = B, ncol = 4L,
               dimnames = list(NULL, c("rd", "rr", "risk_q25", "risk_q75")))
fails <- 0L
set.seed(SEED + 1L)
for (i in seq_len(B)) {
  ans <- tryCatch(risk_contrast(cc30[sample.int(nrow(cc30), nrow(cc30), replace = TRUE), ]),
                  error = function(e) NULL)
  if (is.null(ans) || any(!is.finite(ans))) { fails <- fails + 1L } else { boot[i, ] <- ans }
}
ci <- apply(boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
absrisk <- data.frame(
  model_id = "ADMINV5_ABSRISK_I2_30d", horizon = "postoperative day 30",
  N = nrow(cc30), events = sum(cc30$event_30d), gv_q25 = K$gv_q25, gv_q75 = K$gv_q75,
  risk_q25 = point["risk_q25"], risk_q75 = point["risk_q75"],
  rd = point["rd"], rd_lo = ci[1, "rd"], rd_hi = ci[2, "rd"],
  rr = point["rr"], rr_lo = ci[1, "rr"], rr_hi = ci[2, "rr"],
  n_boot = B, n_boot_failed = fails, outcome_version = OUTCOME_VERSION,
  time_rule = TIME_RULE, stringsAsFactors = FALSE
)
write.csv(absrisk, file.path(OUT_ROOT, "30d_absolute_risk_v5.csv"), row.names = FALSE)
print(absrisk)

# ---- Correct 48-hour landmark sensitivity ----
d48 <- read.csv(file.path(SOURCE_ROOT, "data", "cohort_48h_landmark.csv"))
d48 <- merge(d48, comorb, by = "subject_id", all.x = TRUE)
d48 <- d48[is.na(d48$death_time_composite) | d48$death_time_composite > d48$landmark_48h, ]
d48 <- admin_day30(d48, "landmark_48h")
d48$gv <- d48$gv_sd
d48$mean_glu <- d48$mean_glucose
d48$gv10 <- d48$gv / 10
d48$log_count <- log(d48$n_glucose_0_48h)
K48 <- list(
  gv_sd = sd(d48$gv),
  mean_glu_knots = unname(quantile(d48$mean_glu, c(0.05, 0.35, 0.65, 0.95))),
  N = nrow(d48), events_30d = sum(d48$event_30d), seed = SEED,
  frame = "INSPIRE_OPEND_48H_LANDMARK_V3"
)
write_json(K48, file.path(OUT_ROOT, "inspire_48h_constants_v5.json"), pretty = TRUE, auto_unbox = TRUE)
rcs48 <- paste0("rms::rcs(mean_glu, c(", paste(format(K48$mean_glu_knots, digits = 10), collapse = ","), "))")
cc48 <- d48[complete.cases(d48[, c("t30", "event_30d", "gv", "gv10", "mean_glu", "age")]), ]
f48 <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs48, " + age")),
              data = cc48, x = TRUE, y = TRUE)
s48 <- summary(f48)
res48 <- data.frame(
  model_id = "ADMINV5_I2_48h_landmark_30d", frame = "INSPIRE_OPEND_48H_LANDMARK_V3",
  outcome_version = OUTCOME_VERSION, time_rule = TIME_RULE,
  N = nrow(cc48), events = sum(cc48$event_30d),
  HR_per10 = s48$coefficients["gv10", "exp(coef)"],
  lo = s48$conf.int["gv10", "lower .95"], hi = s48$conf.int["gv10", "upper .95"],
  P = s48$coefficients["gv10", "Pr(>|z|)"],
  HR_perSD = exp(s48$coefficients["gv10", "coef"] * K48$gv_sd / 10),
  ph_global = tryCatch(cox.zph(f48)$table["GLOBAL", "p"], error = function(e) NA_real_),
  note = "0–48-hour exposure with risk beginning at the 48-hour landmark",
  stringsAsFactors = FALSE
)
write.csv(res48, file.path(OUT_ROOT, "48h_landmark_result_v5.csv"), row.names = FALSE)
print(res48)

# ---- 30-day resolution, anchor, and limited joint-model checks ----
sens_rows <- list()
add_sens <- function(...) sens_rows[[length(sens_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
fit_simple <- function(dd, xvar, divisor, id, label) {
  dd$x <- dd[[xvar]] / divisor
  dd <- dd[complete.cases(dd[, c("t30", "event_30d", "x", "mean_glu", "age")]), ]
  if (nrow(dd) < 100L || sum(dd$event_30d) < 8L) {
    add_sens(model_id = id, cohort = label, N = nrow(dd), events = sum(dd$event_30d),
             outcome_version = OUTCOME_VERSION, time_rule = TIME_RULE,
             note = "Not reliably estimable (<100 patients or <8 events)")
    return(invisible(NULL))
  }
  fit <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ x + ", rcs_mean, " + age")), data = dd)
  s <- summary(fit)
  add_sens(model_id = id, cohort = label, N = nrow(dd), events = sum(dd$event_30d),
           HR = s$coefficients["x", "exp(coef)"], lo = s$conf.int["x", "lower .95"],
           hi = s$conf.int["x", "upper .95"], P = s$coefficients["x", "Pr(>|z|)"],
           outcome_version = OUTCOME_VERSION, time_rule = TIME_RULE, note = "")
}

for (rule in list(
  list(name = "ge2", keep = d$n_gv >= 2L),
  list(name = "ge3tp", keep = d$n_gv >= 3L),
  list(name = "ge3tp_ge3val", keep = d$n_gv >= 3L & d$n_distinct_values >= 3L)
)) {
  for (metric in list(
    list(var = "gv_sd", divisor = 10, tag = "gv"),
    list(var = "arv", divisor = 10, tag = "arv"),
    list(var = "mad_glucose", divisor = 10, tag = "mad"),
    list(var = "iqr_glucose", divisor = 10, tag = "iqr")
  )) {
    fit_simple(d[rule$keep, ], metric$var, metric$divisor,
               paste0("ADMINV5_RES_", rule$name, "_", metric$tag, "_30d"),
               paste0(rule$name, "; ", metric$var))
  }
}

# The 48-hour window is intentionally excluded from this 24-hour-risk-origin
# anchor panel; its valid 48-hour-landmark analysis is written above.
for (a in c("opend_0_24h", "opstart_0_24h", "orout_0_24h", "icuin_0_24h", "opend_0_12h")) {
  af <- anchor_features[anchor_features$anchor == a & anchor_features$n >= 2L, ]
  dd <- merge(d, af[, c("subject_id", "mean_glucose", "gv_sd")], by = "subject_id",
              suffixes = c("", "_anchor"))
  dd$gv_anchor <- dd$gv_sd_anchor
  dd$mean_glu <- dd$mean_glucose_anchor
  fit_simple(dd, "gv_anchor", 10, paste0("ADMINV5_ANCH_", a, "_30d"), paste0("anchor=", a))
}

joint <- d[!is.na(d$hba1c_pct) & !is.na(d$eag_mg_dl) & d$eag_mg_dl > 0, ]
joint$shr <- joint$mean_glu / joint$eag_mg_dl
JK <- fromJSON(file.path(SOURCE_ROOT, "results", "inspire_joint_constants.json"))
joint$shr_zf <- (joint$shr - JK$shr_mean) / JK$shr_sd
joint$gv_zf <- (joint$gv - JK$gv_mean) / JK$gv_sd
joint <- joint[complete.cases(joint[, c("t30", "event_30d", "shr_zf", "gv_zf", "mean_glu", "age", "diabetes")]), ]
if (nrow(joint) >= 100L && sum(joint$event_30d) >= 8L) {
  fj <- coxph(Surv(t30, event_30d) ~ shr_zf + gv_zf + age + diabetes, data = joint)
  sj <- summary(fj)
  for (term in c("shr_zf", "gv_zf")) {
    add_sens(model_id = paste0("ADMINV5_JOINT_", term, "_30d"), cohort = "limited SHR-GV joint model",
             N = nrow(joint), events = sum(joint$event_30d), term = term,
             HR = sj$coefficients[term, "exp(coef)"], lo = sj$conf.int[term, "lower .95"],
             hi = sj$conf.int[term, "upper .95"], P = sj$coefficients[term, "Pr(>|z|)"],
             outcome_version = OUTCOME_VERSION, time_rule = TIME_RULE,
             note = "Event-limited; not part of the GV main hierarchy")
  }
} else {
  add_sens(model_id = "ADMINV5_JOINT_30d", cohort = "limited SHR-GV joint model",
           N = nrow(joint), events = sum(joint$event_30d), outcome_version = OUTCOME_VERSION,
           time_rule = TIME_RULE, note = "Not reliably estimable (<100 patients or <8 events)")
}

sens <- rbind_fill(sens_rows)
write.csv(sens, file.path(OUT_ROOT, "30d_sensitivity_joint_results_v5.csv"), row.names = FALSE)
print(sens[, intersect(c("model_id", "N", "events", "HR", "lo", "hi", "P", "note"), names(sens))])

# ---- Auditable time-rule and coverage QC ----
v4_time <- (pmin(ifelse(d$event_30d == 1L, d$death_time_composite, d$discharge_time),
                d$day30_cutoff) - d$landmark_24h) / 1440.0
postdischarge_non_events <- d$event_30d == 0L & d$discharge_time < d$day30_cutoff
qc <- data.frame(
  check = c(
    "coverage_assumption_recorded", "postlandmark_frame", "reconciled_30d_events",
    "all_event_times_equal_composite_death_time", "all_non_event_times_equal_day30_admin_end",
    "postdischarge_deaths_retained_to_death_time", "postdischarge_non_events_retained_to_day30",
    "no_discharge_censoring_in_30d_models", "formula_variable_cc_mask",
    "48h_landmark_frame", "48h_reconciled_30d_events",
    "48h_event_times_equal_composite_death_time", "48h_non_event_times_equal_day30_admin_end",
    "invalid_365d_cox_withdrawn"
  ),
  status = c("ATTESTED", rep("PASS", 13L)),
  value = c(
    "Linked mortality coverage through postoperative day 30 attested by the study team/data custodian",
    nrow(d), sum(d$event_30d), sum(d$event_30d), sum(d$event_30d == 0L),
    sum(d$event_30d == 1L & d$death_time_composite > d$discharge_time, na.rm = TRUE),
    sum(postdischarge_non_events), sum(abs(d$t30 - v4_time) >= 1e-10 & d$event_30d == 0L),
    nrow(cc30), nrow(d48), sum(d48$event_30d), sum(d48$event_30d), sum(d48$event_30d == 0L),
    "No 365-day Cox/absolute-risk result carried forward"
  ),
  detail = c(
    "Raw export lacks registry_coverage_end; this is an explicit source-governance attestation, not an inferred field",
    "Composite death later than the 24-hour landmark or no recorded death", "Deaths by postoperative day 30",
    "Every event uses reconciled composite death time", "Every non-event is administratively censored at day 30",
    "Post-discharge deaths are not moved to discharge", "No-death patients discharged before day 30 remain in the risk set",
    "Count of non-events whose v5 time correctly differs from legacy discharge-censoring time",
    "All stated I1-I3 variables complete; ASA/BMI not used in 30-day core formulas",
    "Correct 48-hour landmark analysis frame", "Deaths by postoperative day 30 after the 48-hour landmark",
    "Every 48-hour-landmark event uses reconciled composite death time",
    "Every 48-hour-landmark non-event is administratively censored at day 30",
    "No common 365-day registry coverage end documented"
  ),
  stringsAsFactors = FALSE
)
write.csv(qc, file.path(OUT_ROOT, "30d_time_rule_qc_v5.csv"), row.names = FALSE)

withdrawal <- data.frame(
  analysis = "INSPIRE 365-day Cox and derived measures",
  status = "WITHDRAWN_FROM_MANUSCRIPT_V5",
  rationale = "v4 used event-dependent discharge censoring and no common 365-day registry coverage end is documented; MICE for ASA/BMI cannot repair follow-up bias",
  replacement = "No 365-day INSPIRE HR, absolute risk, RMST, calibration, or performance claim",
  stringsAsFactors = FALSE
)
write.csv(withdrawal, file.path(OUT_ROOT, "365d_withdrawal_v5.csv"), row.names = FALSE)

manifest <- data.frame(
  version = "v5", outcome_version = OUTCOME_VERSION, time_rule = TIME_RULE,
  primary_N = nrow(cc30), primary_events = sum(cc30$event_30d),
  registry_coverage_note = "Declared adequate through postoperative day 30; no raw coverage-end field exported",
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(OUT_ROOT, "v5_analysis_manifest.csv"), row.names = FALSE)

print(qc)
cat("V5_ADMIN_CENSORING_RERUN_DONE\n")
sink()
