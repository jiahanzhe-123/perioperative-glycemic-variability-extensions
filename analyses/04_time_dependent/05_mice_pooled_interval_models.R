# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 05_mice_pooled_interval_models.R — 365d PH 处理(MICE m=50 池化)。
# 迁移自 stats_fix_phase1/scripts/01_mimic_ph_repairs_mice_pooled.R(逻辑不变,路径配置化)。
#   A) 365d interval-specific GV estimates (d1-7 / d8-30 / d31-365), MICE-pooled, with 95% CI.
#   B) 365d GV x log(time) continuous time interaction, MICE-pooled (sensitivity).
#   C) 30d Model B + prespecified log(time) terms for violating covariates (age, Charlson), MICE-pooled.
#   S) Sanity refit: MICE Model B 30d/365d must reproduce frozen primary values.
# 输入:mimic_record_work 的 mice_m50_object.rds 与 standardization_constants.json
# (由 analyses/03_primary_mimic/03_run_primary_models_mice.R 在 BMI 修复帧上产生)。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
# rm(list=ls()) 会清除页首 source 的配置,须重新加载:
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
P35  <- PGV("mimic_record_work")
dir.create(file.path(ROOT,"logs"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(ROOT,"logs","01_mimic_ph_repairs.log"), split=TRUE)

K <- fromJSON(file.path(P35,"results","standardization_constants.json"))
imp <- readRDS(file.path(P35,"results","mice_m50_object.rds"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")

rubin <- function(betas, vars){
  m <- length(betas); qb <- mean(unlist(betas)); U <- mean(unlist(vars))
  B <- var(unlist(betas)); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, fmi=(1+1/m)*B/Tm)
}
grab_pool <- function(bs, vs){
  rb <- rubin(bs, vs); q <- qt(.975, rb$df)
  list(HR=exp(rb$beta), lo=exp(rb$beta-q*rb$se), hi=exp(rb$beta+q*rb$se),
       P=2*pt(abs(rb$beta/rb$se), df=rb$df, lower.tail=FALSE), fmi=rb$fmi, m=length(bs))
}

# ---------- S) sanity refit: frozen primary must reproduce ----------
cat("== sanity: MICE Model B refit ==\n")
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  fml <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml))
  bs <- c(); vs <- c()
  for (i in 1:50) {
    dd <- complete(imp, i); dd$gv10 <- dd$gv/10
    s <- summary(coxph(fml, data=dd))
    bs <- c(bs, s$coefficients["gv10","coef"]); vs <- c(vs, s$coefficients["gv10","se(coef)"]^2)
  }
  r <- grab_pool(bs, vs)
  cat(sprintf("ModelB_%sd: HR=%.6f CI=(%.6f, %.6f) P=%.6f [frozen 30d: 0.979381 (0.922463, 1.039810) P=0.495213; 365d: 0.984981 (0.941996, 1.029927) P=0.506233]\n",
              hz, r$HR, r$lo, r$hi, r$P))
}

# ---------- A) 365d interval-specific GV estimates, MICE-pooled ----------
cat("== A) 365d GV intervals ==\n")
tm <- "t_lm_365"; ev <- "event_lm_365"
INTS <- c("d1-7","d8-30","d31-365")
pool_int <- setNames(vector("list",3), INTS)
for (lv in INTS) pool_int[[lv]] <- list(b=c(), v=c())
n_at_risk <- setNames(rep(NA_integer_,3), INTS); ev_int <- setNames(rep(NA_integer_,3), INTS)
fml_int <- as.formula(paste0("Surv(tstart,tstop,evsplit) ~ gv10:interval + ", rcs_mean, " + ", covs_fml))
for (i in 1:50) {
  dd <- complete(imp, i); dd$gv10 <- dd$gv/10
  sp <- survSplit(as.formula(paste0("Surv(",tm,",",ev,") ~ .")),
                  data=dd[, c("stay_id",tm,ev,"gv10","mean_glu",COVS)],
                  cut=c(7,30), zero=-0.5, episode="interval", start="tstart", end="tstop", event="evsplit")
  sp$interval <- factor(sp$interval, labels=INTS)
  sp <- sp[!is.na(sp$evsplit),]
  if (i==1) {
    ev_int <- tapply(sp$evsplit, sp$interval, sum)[INTS]
    n_at_risk <- tapply(seq_len(nrow(sp)), sp$interval, function(ix) length(unique(sp$stay_id[ix])))[INTS]
  }
  fit <- tryCatch(coxph(fml_int, data=sp), error=function(e) NULL)
  if (is.null(fit)) { cat("imputation", i, "interval fit failed\n"); next }
  si <- summary(fit)$coefficients
  for (lv in INTS) {
    term <- paste0("gv10:interval", lv)
    pool_int[[lv]]$b <- c(pool_int[[lv]]$b, si[term,"coef"])
    pool_int[[lv]]$v <- c(pool_int[[lv]]$v, si[term,"se(coef)"]^2)
  }
}
rowsA <- list()
for (lv in INTS) {
  r <- grab_pool(pool_int[[lv]]$b, pool_int[[lv]]$v)
  rowsA[[lv]] <- data.frame(model_id="MICE_B_365d_GV_intervals", interval=lv,
    n_at_risk_start=as.integer(n_at_risk[lv]), interval_events=as.integer(ev_int[lv]),
    HR_per10=r$HR, lo=r$lo, hi=r$hi, P=r$P, fmi=r$fmi, n_imputations=r$m, stringsAsFactors=FALSE)
  cat(sprintf("%s: N_at_risk=%d events=%d HR=%.4f (%.4f-%.4f) P=%.4f fmi=%.4f\n",
              lv, n_at_risk[lv], ev_int[lv], r$HR, r$lo, r$hi, r$P, r$fmi))
}
tabA <- do.call(rbind, rowsA)
write.csv(tabA, file.path(ROOT,"results","PH365_INTERVAL_MICE_POOLED.csv"), row.names=FALSE)

# ---------- B) 365d GV x log(time) sensitivity, MICE-pooled ----------
cat("== B) 365d GV x log(time) ==\n")
fml_tt365 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + tt(gv10) + ", rcs_mean, " + ", covs_fml))
poolB <- list(gv10=list(b=c(),v=c()), `tt(gv10)`=list(b=c(),v=c()))
for (i in 1:50) {
  dd <- complete(imp, i); dd$gv10 <- dd$gv/10
  fit <- tryCatch(coxph(fml_tt365, data=dd, tt=function(x,t,...) x*log(pmax(t,0.5))), error=function(e) NULL)
  if (is.null(fit)) { cat("imputation", i, "tt365 fit failed\n"); next }
  si <- summary(fit)$coefficients
  for (tn in names(poolB)) { poolB[[tn]]$b <- c(poolB[[tn]]$b, si[tn,"coef"]); poolB[[tn]]$v <- c(poolB[[tn]]$v, si[tn,"se(coef)"]^2) }
}
rowsB <- list()
for (tn in names(poolB)) {
  rb <- rubin(poolB[[tn]]$b, poolB[[tn]]$v); q <- qt(.975, rb$df)
  rowsB[[tn]] <- data.frame(model_id="MICE_B_365d_GV_logtime", term=tn,
    beta=rb$beta, se=rb$se, HR=exp(rb$beta), lo=exp(rb$beta-q*rb$se), hi=exp(rb$beta+q*rb$se),
    P=2*pt(abs(rb$beta/rb$se), df=rb$df, lower.tail=FALSE), fmi=rb$fmi, n_imputations=length(poolB[[tn]]$b), stringsAsFactors=FALSE)
  cat(sprintf("%s: beta=%.5f HR=%.4f (%.4f-%.4f) P=%.4f\n", tn, rb$beta, exp(rb$beta), exp(rb$beta-q*rb$se), exp(rb$beta+q*rb$se), rowsB[[tn]]$P))
}
tabB <- do.call(rbind, rowsB)
write.csv(tabB, file.path(ROOT,"results","PH365_GV_LOGTIME_SENSITIVITY_MICE_POOLED.csv"), row.names=FALSE)

# ---------- C) 30d Model B + prespecified log(time) terms for age + Charlson, MICE-pooled ----------
cat("== C) 30d tt-corrected (age, charlson) ==\n")
tm <- "t_lm_30"; ev <- "event_lm_30"
viol <- c("age_at_admission","charlson_without_diabetes")  # post-repair CC PH: P=0.0451 / 0.0415
fml_tt30 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml, " + ",
                               paste(paste0("tt(", viol, ")"), collapse=" + ")))
bs <- c(); vs <- c()
for (i in 1:50) {
  dd <- complete(imp, i); dd$gv10 <- dd$gv/10
  fit <- tryCatch(coxph(fml_tt30, data=dd, tt=function(x,t,...) x*log(pmax(t,0.5))), error=function(e) NULL)
  if (is.null(fit)) { cat("imputation", i, "tt30 fit failed\n"); next }
  si <- summary(fit)$coefficients
  bs <- c(bs, si["gv10","coef"]); vs <- c(vs, si["gv10","se(coef)"]^2)
}
r <- grab_pool(bs, vs)
tabC <- data.frame(model_id="MICE_B_30d_tt_corrected",
  term="age_at_admission+charlson_without_diabetes",
  N=10561, events=296, HR_per10=r$HR, lo=r$lo, hi=r$hi, P=r$P, fmi=r$fmi, n_imputations=r$m,
  note="Model B with prespecified time-varying log(time) terms for PH-violating covariates (age_at_admission, charlson_without_diabetes); GV effect pooled over m=50",
  stringsAsFactors=FALSE)
write.csv(tabC, file.path(ROOT,"results","PH30_MODELB_TT_CORRECTED_MICE_POOLED.csv"), row.names=FALSE)
cat(sprintf("30d tt-corrected GV: HR=%.4f (%.4f-%.4f) P=%.4f fmi=%.4f m=%d\n", r$HR, r$lo, r$hi, r$P, r$fmi, r$m))
cat("PHASE1_MIMIC_DONE\n")
sink()
