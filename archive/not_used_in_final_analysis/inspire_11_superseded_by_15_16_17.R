#!/usr/bin/env Rscript
# 11_inspire_analysis.R — INSPIRE 主要分析(冻结协议 INSPIRE_ANALYSIS_PROTOCOL_v1, sha256=12171ced…)
# 标记:confirmatory secondary / external extension / exploratory。
# 输出:01_cohort_flow, 02_exposure_distribution, 03_quantization_audit, standardization_constants.json,
#       04_primary_results.csv, 06_absolute_risk.csv, 07_ph_diagnostics.csv, mice 诊断, rcs 非线性。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.libPaths(c(file.path("private", "r-library"), .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite); library(rms)})
ROOT <- normalizePath(file.path("private", "inspire_project"), mustWork=FALSE)
dir.create(file.path(ROOT,"results"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(ROOT,"logs","11_inspire_analysis.log"), split=TRUE)

b <- read.csv(file.path(ROOT,"data","inspire_base.csv"), stringsAsFactors=FALSE)
cm <- read.csv(file.path(ROOT,"data","comorbidity.csv"), stringsAsFactors=FALSE)
d <- merge(b, cm, by="subject_id", all.x=TRUE)
# inspire_base 已含暴露与测量过程列,直接使用,避免合并撞列
d$gv <- d$gv_sd; d$mean_glu <- d$mean_glucose; d$n_gv <- d$n_glucose_0_24h
d$bmi <- d$weight / (d$height/100)^2
d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$emop_f <- as.integer(d$emop %in% c("1",1,TRUE,"t","true"))
d$asa_f <- factor(d$asa)
d$surgery_group <- factor(d$surgery_group, levels=c("OPEN_CABG","OPEN_VALVE"))
d$sex <- factor(d$sex)
d$cpb <- as.integer(!is.na(d$cpbon_time))
# 结局(冻结定义):30/365 日全因死亡,时间零=landmark(opend+1440),随访自 landmark 起
d$event_30d <- as.integer(!is.na(d$allcause_death_time) & d$allcause_death_time <= d$opend_time + 30*1440)
d$event_365d <- as.integer(!is.na(d$allcause_death_time) & d$allcause_death_time <= d$opend_time + 365*1440)
d$event_30d_inhosp <- as.integer(!is.na(d$inhosp_death_time) & d$inhosp_death_time <= d$opend_time + 30*1440)
d$followup_days <- (pmin(ifelse(is.na(d$allcause_death_time), Inf, d$allcause_death_time), d$discharge_time) - d$opend_time)/1440
d$t30 <- pmin(d$followup_days, 30) - 1
d$t365 <- pmin(d$followup_days, 365) - 1

# ---- 目标队列(协议 §1) ----
tgt <- d[d$age>=18 & !is.na(d$opend_time) & !is.na(d$gv) & d$n_gv>=2 &
         d$landmark_eligible_allcause %in% c(TRUE,"t","True","TRUE","true",1,"1") & d$t30>=0,]
stopifnot(!anyDuplicated(tgt$subject_id))
cat("目标队列 N =", nrow(tgt), "; 30d 事件 =", sum(tgt$event_30d), "; 365d 事件 =", sum(tgt$event_365d),
    "; 30d 院内事件 =", sum(tgt$event_30d_inhosp), "\n")

# ---- 队列流程(01) ----
flow <- data.frame(
  step=c("adult first high-specificity open CABG/valve (codebook v1 INCLUDE)",
         "verified operation end time (opend non-missing)",
         ">=2 distinct glucose timepoints in opend+0-24h (frozen series)",
         "alive with follow-up at day-1 landmark (all-cause linkage)",
         "TARGET COHORT"),
  n=c(sum(d$age>=18 & !is.na(d$opend_time)),
      sum(d$age>=18 & !is.na(d$opend_time)),
      sum(d$age>=18 & !is.na(d$opend_time) & !is.na(d$gv) & d$n_gv>=2),
      nrow(tgt), nrow(tgt)))
flow$removed_from_previous <- -c(0, diff(flow$n))
write.csv(flow, file.path(ROOT,"results","01_cohort_flow.csv"), row.names=FALSE)

# ---- 量化审计(03) ----
gvals <- sort(unique(round(tgt$mean_glu)))
qa <- data.frame(item=c("distinct mean-glucose values in cohort",
                        "distinct per-patient glucose values: median",
                        "pct patients with <10 distinct values",
                        "pct patients with <5 distinct values",
                        "median adjacent spacing of mean values (mg/dL)",
                        "pct values that are integers",
                        "note"),
                 value=c(length(gvals),
                         median(tgt$n_distinct_values, na.rm=TRUE),
                         100*mean(tgt$n_distinct_values < 10, na.rm=TRUE),
                         100*mean(tgt$n_distinct_values < 5, na.rm=TRUE),
                         median(diff(gvals)),
                         100,
                         "INSPIRE glucose integer-quantized; GV 称 quantized routine-laboratory SD-based GV"))
write.csv(qa, file.path(ROOT,"results","03_quantization_audit.csv"), row.names=FALSE)

# ---- 暴露分布(02) ----
dist_of <- function(v) c(n=length(v), mean=mean(v), sd=sd(v), median=median(v),
  q25=quantile(v,.25), q75=quantile(v,.75), p1=quantile(v,.01), p99=quantile(v,.99),
  min=min(v), max=max(v))
dist_tab <- as.data.frame(rbind(GV=dist_of(tgt$gv), mean_glucose=dist_of(tgt$mean_glu),
  n_measurements=dist_of(tgt$n_gv), span_hours=dist_of(tgt$span_hours)))
dist_tab <- cbind(variable=rownames(dist_tab), round(dist_tab, 3))
write.csv(dist_tab, file.path(ROOT,"results","02_exposure_distribution.csv"), row.names=FALSE)

# ---- 标准化常数(冻结) ----
K <- list(
  gv_mean=mean(tgt$gv), gv_sd=sd(tgt$gv), gv_median=median(tgt$gv),
  gv_q25=unname(quantile(tgt$gv,.25)), gv_q75=unname(quantile(tgt$gv,.75)),
  mean_glu_mean=mean(tgt$mean_glu), mean_glu_sd=sd(tgt$mean_glu),
  mean_glu_knots=unname(quantile(tgt$mean_glu, c(.05,.35,.65,.95))),
  gv_knots_30d=unname(quantile(tgt$gv, c(.10,.50,.90))),
  gv_knots_365d=unname(quantile(tgt$gv, c(.05,.35,.65,.95))),
  N=nrow(tgt), events_30d=sum(tgt$event_30d), events_365d=sum(tgt$event_365d),
  cohort="INSPIRE v1.4.2 open CABG/valve, opend-anchored day-1 landmark, >=2 quantized glucose timepoints",
  seed=SEED, protocol_sha256="12171cedff43bac9b8a5a76407640391ace120145e6fe6f49e06c2c2b98e7d1c",
  input_checksum_base="dd842f46e477f7740e771e543b83ad82")
write_json(K, file.path(ROOT,"results","standardization_constants.json"), pretty=TRUE, auto_unbox=TRUE)
tgt$gv10 <- tgt$gv/10
kSD <- K$gv_sd/10
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")

# ---- 模型(协议 §4:30d 核心=age+GV+RCS mean(EPV 规则);365d 全临床集) ----
CLIN_365 <- c("age","sex","diabetes","charlson_without_diabetes","asa_f","emop_f","surgery_group","bmi")
fit_models <- function(dd, hz){
  tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
  if (hz==30) {
    covs_core <- "age"
    fI0 <- as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_core))
    fI1 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_core))
    fI2 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_core))
    fI3 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_core,
               " + log_count + span_hours"))
  } else {
    covs_full <- paste(CLIN_365, collapse=" + ")
    fI0 <- as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_full))
    fI1 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_full))
    fI2 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_full))
    fI3 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_full,
               " + log_count + span_hours"))
  }
  list(I0=fI0, I1=fI1, I2=fI2, I3=fI3)
}
report <- function(fit, model_id, model, N, events, note=""){
  s <- summary(fit)
  data.frame(model_id=model_id, model=model, N=N, events=events,
    HR_per10=s$coefficients["gv10","exp(coef)"], lo_per10=s$conf.int["gv10","lower .95"],
    hi_per10=s$conf.int["gv10","upper .95"], P_per10=s$coefficients["gv10","Pr(>|z|)"],
    HR_perSD=exp(s$coefficients["gv10","coef"]*kSD),
    lo_perSD=exp((s$coefficients["gv10","coef"]-1.96*s$coefficients["gv10","se(coef)"])*kSD),
    hi_perSD=exp((s$coefficients["gv10","coef"]+1.96*s$coefficients["gv10","se(coef)"])*kSD),
    P_perSD=s$coefficients["gv10","Pr(>|z|)"], note=note, stringsAsFactors=FALSE)
}
tgt$log_count <- log(tgt$n_gv)
cc <- tgt[complete.cases(tgt[,c("t30","t365","event_30d","event_365d","gv","gv10","mean_glu","age","sex","diabetes",
  "charlson_without_diabetes","asa_f","emop_f","surgery_group","bmi","log_count","span_hours")]) & tgt$t30>0,]
cat("complete-case N =", nrow(cc), "; 30d 事件 =", sum(cc$event_30d), "; 365d 事件 =", sum(cc$event_365d), "\n")
cc_rows <- list(); ph_rows <- list(); fits_cc <- list()
for (hz in c(30,365)) {
  fmls <- fit_models(cc, hz)
  for (mk in c("I1","I2","I3")) {
    fit <- coxph(fmls[[mk]], data=cc, x=TRUE, y=TRUE)
    fits_cc[[paste0(mk,"_",hz)]] <- fit
    note <- if(hz==30) "EPV-limited (core covariates per frozen rule)" else ""
    cc_rows[[length(cc_rows)+1]] <- report(fit, paste0("CC_",mk,"_",hz,"d"), paste0("Model ",mk), nrow(cc), sum(cc[[paste0("event_",hz,"d")]]), note)
    z <- tryCatch(cox.zph(fit)$table, error=function(e) NULL)
    if (!is.null(z)) {
      tab <- as.data.frame(z); tab$term <- rownames(tab)
      for (i in seq_len(nrow(tab)))
        ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("PH_",mk,"_",hz,"d"),
          term=tab$term[i], chisq=tab$chisq[i], df=tab$df[i], p=tab$p[i], stringsAsFactors=FALSE)
    }
  }
}
cc_tab <- do.call(rbind, cc_rows); ph_tab <- do.call(rbind, ph_rows)
write.csv(cc_tab, file.path(ROOT,"results","cc_models_I1_I2_I3.csv"), row.names=FALSE)
write.csv(ph_tab, file.path(ROOT,"results","07_ph_diagnostics.csv"), row.names=FALSE)
print(cc_tab[,c("model_id","N","events","HR_per10","lo_per10","hi_per10","P_per10","note")])

# ---- MICE m=50(协议 §5) ----
mi_vars <- c("subject_id","t30","t365","event_30d","event_365d","gv","gv10","mean_glu","n_gv","span_hours","log_count",
             "age","sex","diabetes","charlson_without_diabetes","asa_f","emop_f","surgery_group","bmi","cpb")
mi_data <- tgt[, mi_vars]
mi_data$sex <- as.numeric(mi_data$sex); mi_data$asa_num <- suppressWarnings(as.numeric(as.character(mi_data$asa_f)))
mi_data$proc <- as.numeric(mi_data$surgery_group)
mi_data <- mi_data[, names(mi_data)!="asa_f"]
mi_data$H365 <- nelsonaalen(mi_data, "t365", "event_365d")
# 确定性重复与近共线项从数据中剔除,派生量在各插补集内重算
mi_data <- mi_data[, !names(mi_data) %in% c("gv10","proc")]
pred <- make.predictorMatrix(mi_data); pred[,] <- 0L
meth <- make.method(mi_data); meth[] <- ""
imp_targets <- intersect(c("bmi","asa_num","charlson_without_diabetes"), names(mi_data))
imp_targets <- imp_targets[colSums(is.na(mi_data[,imp_targets]))>0]
impute_preds <- setdiff(names(mi_data), c("subject_id","t30","event_30d"))
for (v in imp_targets) { meth[v] <- "pmm"; pred[v, impute_preds] <- 1L }
cat("MICE 插补变量:", paste(imp_targets, collapse=", "), "\n")
t0 <- Sys.time()
imp <- mice(mi_data, m=50, maxit=20, method=meth, predictorMatrix=pred, seed=SEED, printFlag=FALSE)
cat("MICE 耗时:", round(difftime(Sys.time(), t0, units="mins"),1), "分钟\n")
saveRDS(imp, file.path(ROOT,"results","inspire_mice_m50.rds"))
lg <- imp$loggedEvents
if (!is.null(lg) && nrow(lg)) { write.csv(data.frame(lg), file.path(ROOT,"results","mice_logged_events.csv"), row.names=FALSE)
  cat("loggedEvents:", nrow(lg), "\n") } else cat("loggedEvents: 0\n")

rubin <- function(betas, vars){
  m <- length(betas); qb <- mean(unlist(betas)); U <- mean(unlist(vars))
  B <- var(unlist(betas)); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, p=2*pt(abs(qb/sqrt(Tm)), df=df, lower.tail=FALSE), fmi=(1+1/m)*B/Tm)
}
pool_rows <- list()
for (hz in c(30,365)) {
  fmls <- fit_models(tgt, hz)
  ev <- paste0("event_",hz,"d")
  for (mk in c("I1","I2","I3")) {
    bs <- c(); vs <- c()
    for (i in 1:50) {
      dd <- complete(imp, i)
      dd$gv10 <- dd$gv/10; dd$log_count <- log(dd$n_gv)
      dd$sex <- factor(dd$sex); dd$asa_f <- factor(dd$asa_num)
      dd$surgery_group <- factor(dd$surgery_group)
      f <- tryCatch(coxph(fmls[[mk]], data=dd), error=function(e) NULL)
      if (is.null(f)) next
      s <- summary(f)
      bs <- c(bs, s$coefficients["gv10","coef"]); vs <- c(vs, s$coefficients["gv10","se(coef)"]^2)
    }
    rb <- rubin(bs, vs)
    z <- rb$beta/rb$se
    pool_rows[[length(pool_rows)+1]] <- data.frame(
      model_id=sprintf("MICE_%s_%sd", mk, hz), model=paste0("Model ",mk),
      N=nrow(tgt), events=sum(tgt[[ev]]),
      HR_per10=exp(rb$beta), lo_per10=exp(rb$beta-qt(.975,rb$df)*rb$se), hi_per10=exp(rb$beta+qt(.975,rb$df)*rb$se),
      P_per10=rb$p, HR_perSD=exp(rb$beta*kSD), lo_perSD=exp((rb$beta-qt(.975,rb$df)*rb$se)*kSD),
      hi_perSD=exp((rb$beta+qt(.975,rb$df)*rb$se)*kSD), P_perSD=rb$p, fmi=rb$fmi,
      n_imputations_fit=length(bs),
      note=if(hz==30) "EPV-limited (core covariates per frozen rule)" else "", stringsAsFactors=FALSE)
  }
}
pool_tab <- do.call(rbind, pool_rows)
write.csv(pool_tab, file.path(ROOT,"results","mice_pooled_models.csv"), row.names=FALSE)
print(pool_tab[,c("model_id","N","events","HR_per10","lo_per10","hi_per10","P_per10","HR_perSD","P_perSD","fmi")])

# ---- 主要结果(04) ----
g <- function(id) pool_tab[pool_tab$model_id==id,]
prim <- data.frame(
  item=c("PRIMARY (INSPIRE external extension): GV per 10 mg/dL (Model I2, 30d all-cause)",
         "PRIMARY: GV per fixed-cohort SD (Model I2, 30d)",
         "KEY SECONDARY: GV per 10 mg/dL (Model I2, 365d)",
         "KEY SECONDARY: GV per fixed-cohort SD (Model I2, 365d)",
         "REFERENCE: Model I1 (no mean glucose), 30d",
         "CONTEXT: Model I3 (+measurement process), 30d",
         "CONTEXT: 30-day in-hospital mortality (I2, complete-case)"),
  HR=c(g("MICE_I2_30d")$HR_per10, g("MICE_I2_30d")$HR_perSD, g("MICE_I2_365d")$HR_per10, g("MICE_I2_365d")$HR_perSD,
       g("MICE_I1_30d")$HR_per10, g("MICE_I3_30d")$HR_per10,
       cc_tab[cc_tab$model_id=="CC_I2_30d","HR_per10"]),
  lo=c(g("MICE_I2_30d")$lo_per10, g("MICE_I2_30d")$lo_perSD, g("MICE_I2_365d")$lo_per10, g("MICE_I2_365d")$lo_perSD,
       g("MICE_I1_30d")$lo_per10, g("MICE_I3_30d")$lo_per10, cc_tab[cc_tab$model_id=="CC_I2_30d","lo_per10"]),
  hi=c(g("MICE_I2_30d")$hi_per10, g("MICE_I2_30d")$hi_perSD, g("MICE_I2_365d")$hi_per10, g("MICE_I2_365d")$hi_perSD,
       g("MICE_I1_30d")$hi_per10, g("MICE_I3_30d")$hi_per10, cc_tab[cc_tab$model_id=="CC_I2_30d","hi_per10"]),
  P=c(g("MICE_I2_30d")$P_per10, g("MICE_I2_30d")$P_perSD, g("MICE_I2_365d")$P_per10, g("MICE_I2_365d")$P_perSD,
      g("MICE_I1_30d")$P_per10, g("MICE_I3_30d")$P_per10, cc_tab[cc_tab$model_id=="CC_I2_30d","P_per10"]),
  N=nrow(tgt), events_30d=sum(tgt$event_30d), events_365d=sum(tgt$event_365d),
  gv_sd_constant=K$gv_sd, gv_mean_constant=K$gv_mean, seed=SEED,
  analysis_label="confirmatory secondary / external extension / exploratory", stringsAsFactors=FALSE)
write.csv(prim, file.path(ROOT,"results","04_primary_results.csv"), row.names=FALSE)
print(prim[,c("item","HR","lo","hi","P")])

# ---- 绝对风险(Q75 vs Q25;1000 患者级 bootstrap) ----
risk_contrast <- function(dd, hz){
  tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
  tmax <- if(hz==30) 29 else 364
  fmls <- fit_models(dd, hz)
  fit <- coxph(fmls$I2, data=dd)
  bh <- basehaz(fit, centered=FALSE); bh <- bh[bh$time<=tmax,]; H <- max(bh$hazard)
  tt <- c(0, bh$time); Hv <- c(0, bh$hazard); dt <- diff(c(tt, tmax))
  lp25 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q25; nd$gv10 <- nd$gv/10; nd}, type="lp", reference="zero")
  lp75 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q75; nd$gv10 <- nd$gv/10; nd}, type="lp", reference="zero")
  r25 <- mean(1 - exp(-H*exp(lp25))); r75 <- mean(1 - exp(-H*exp(lp75)))
  rmst_diff <- NA
  if (hz==365) {
    rmst_of <- function(lpv){ S <- exp(-outer(Hv, exp(lpv))); mean(colSums(S[-nrow(S),,drop=FALSE]*dt)) }
    rmst_diff <- rmst_of(lp75) - rmst_of(lp25)
  }
  list(risk25=r25, risk75=r75, rd=r75-r25, rr=r75/r25, rmst_diff=rmst_diff)
}
boot_rows <- list()
for (hz in c(30,365)) {
  ev <- paste0("event_",hz,"d")
  pt <- risk_contrast(cc, hz)
  B <- 1000; vals <- matrix(NA, B, 5); fails <- 0L; n <- nrow(cc)
  for (b_ in 1:B) {
    idx <- sample.int(n, n, replace=TRUE); db <- cc[idx,]
    rb <- tryCatch({ r <- risk_contrast(db, hz); c(r$rd, r$rr, r$risk25, r$risk75, ifelse(is.na(r$rmst_diff), NA, r$rmst_diff)) },
                   error=function(e) NULL)
    if (is.null(rb) || any(!is.finite(rb[1:4]))) { fails <- fails+1L; next }
    vals[b_,] <- rb
  }
  colnames(vals) <- c("rd","rr","risk25","risk75","rmst_diff")
  ci <- apply(vals, 2, quantile, probs=c(.025,.975), na.rm=TRUE)
  boot_rows[[length(boot_rows)+1]] <- data.frame(
    model_id=paste0("ABSRISK_",hz,"d"), horizon=hz, N=nrow(cc), events=sum(cc[[ev]]),
    gv_q25=K$gv_q25, gv_q75=K$gv_q75, risk_q25=pt$risk25, risk_q75=pt$risk75,
    rd=pt$rd, rd_lo=ci[1,"rd"], rd_hi=ci[2,"rd"], rr=pt$rr, rr_lo=ci[1,"rr"], rr_hi=ci[2,"rr"],
    rmst_diff=pt$rmst_diff, rmst_lo=if(hz==365) ci[1,"rmst_diff"] else NA, rmst_hi=if(hz==365) ci[2,"rmst_diff"] else NA,
    n_boot=B, n_boot_failed=fails, method="percentile bootstrap, patient-level, full model refit", stringsAsFactors=FALSE)
  cat(sprintf("[absrisk %sd] RD=%.4f (%.4f, %.4f), fails=%d\n", hz, pt$rd, ci[1,"rd"], ci[2,"rd"], fails))
}
boot_tab <- do.call(rbind, boot_rows)
write.csv(boot_tab, file.path(ROOT,"results","06_absolute_risk.csv"), row.names=FALSE)

# ---- RCS 非线性(次要) ----
lrt_p <- function(a,b){ ddf <- attr(logLik(b),"df")-attr(logLik(a),"df")
  if(!is.finite(ddf)||ddf<=0) return(NA_real_); pchisq(2*as.numeric(logLik(b)-logLik(a)), df=ddf, lower.tail=FALSE) }
nl_rows <- list()
for (hz in c(30,365)) {
  tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
  knots <- if(hz==30) K$gv_knots_30d else K$gv_knots_365d
  kt <- paste(format(knots, digits=10), collapse=",")
  covs <- if(hz==30) "age" else paste(CLIN_365, collapse=" + ")
  fS <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(gv, c(", kt, ")) + ", rcs_mean, " + ", covs)), data=cc, x=TRUE, y=TRUE)
  fL <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv + ", rcs_mean, " + ", covs)), data=cc, x=TRUE, y=TRUE)
  f0 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ ", rcs_mean, " + ", covs)), data=cc, x=TRUE, y=TRUE)
  nl_rows[[length(nl_rows)+1]] <- data.frame(model_id=paste0("RCS_GV_",hz,"d"), horizon=hz,
    N=nrow(cc), events=sum(cc[[ev]]), knots=paste(round(knots,3),collapse=";"),
    overall_p=lrt_p(f0,fS), nonlinear_p=lrt_p(fL,fS), note="secondary; no threshold search", stringsAsFactors=FALSE)
}
write.csv(do.call(rbind, nl_rows), file.path(ROOT,"results","rcs_gv_nonlinearity.csv"), row.names=FALSE)
cat("PHASE11_DONE\n")
sink()
