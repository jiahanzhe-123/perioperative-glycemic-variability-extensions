# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 03_alt_gv_bloodonly_framefix.R
# 修复 reviewer 问题①:替代 GV 指标分析
# 修复点:
#   1. 完整病例 frame 与主模型完全一致(补回 frame_extra:GV-only 两个 horizon 及
#      SHR-GV 30d 额外要求 bmi+diabetes 完整)——原 03_alt_gv_bloodonly.R 缺失,
#      导致 N=5013/10325 而非主模型的 4869/9786。
#   2. 非 GV-only 队列所有模型联合纳入 shr_z(原脚本未纳入,无法复现主模型)。
#   3. 每个 (cohort,horizon) 第一行为主模型复现行,并与
#      final_primary_models.csv 自动核对(N/events 精确一致,HR 容差 1e-6),不一致即 stop。
#   4. 每行输出元数据:精确队列/协变量规格/SHR是否同调/标准化样本/缺失规则。
# 输入与暴露版本不变(final blood-only sequence),不覆盖原脚本与原输出。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr)})
ROOT <- PGV("mimic_work")
OUT  <- file.path(ROOT,"reviewer_issue_verification_20260728/01_frame_repair/alt_gv")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sink(file.path(ROOT,"reviewer_issue_verification_20260728/logs/03_alt_gv_framefix.log"), split=TRUE)

dat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA","NULL"))
num <- c("shr","gv","age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first",
         "creat_postop_first","wbc_postop_first","albumin_adm_first","hgb_postop_first",
         "platelets_postop_first","sofa_24h","survival_time_days","postop_30d_death_flag","one_year_death_flag",
         "diabetes","glucose_cv_postop_24h","glucose_mean_postop_24h","arv","tw_arv","n_glucose_postop_24h")
for (nm in num) dat[[nm]] <- suppressWarnings(as.numeric(as.character(dat[[nm]])))
dat$cv_pct <- dat$glucose_cv_postop_24h * 100
dat$gender <- factor(dat$gender)
pg <- tolower(trimws(dat$procedure_group_main))
dat$procedure_group_main_model <- factor(ifelse(pg=="cabg","cabg",ifelse(pg=="open_valve","open_valve",
  ifelse(pg=="transplant_vad","transplant_vad","aortic_congenital_other"))),
  levels=c("cabg","open_valve","aortic_congenital_other","transplant_vad"))
dat$event_30d <- as.integer(dat$postop_30d_death_flag==1)
dat$event_1y <- as.integer(dat$one_year_death_flag==1)
dat$t30 <- pmin(dat$survival_time_days,30); dat$t365 <- pmin(dat$survival_time_days,365)
reduced <- "age_at_admission + gender + charlson_comorbidity_index + procedure_group_main_model + lactate_postop_first + creat_postop_first + sofa_24h"
full <- paste(reduced, "+ bmi + wbc_postop_first + albumin_adm_first + hgb_postop_first + platelets_postop_first")
safe_z <- function(x) as.numeric(scale(suppressWarnings(as.numeric(x))))

# ---- 与 01_refit_primary_models.R 完全一致的 frame 规则 ----
build_primary_frame <- function(dat, cohort_flag, cohort, horizon){
  d <- dat[dat[[cohort_flag]]==1,]
  ev <- if(horizon==30) "event_30d" else "event_1y"; tm <- if(horizon==30) "t30" else "t365"
  covs <- if(cohort=="GV-only") reduced else if(horizon==30) reduced else full
  fe <- if(cohort=="GV-only" || (cohort=="SHR-GV" && horizon==30)) c("bmi","diabetes") else NULL
  expo_vars <- if(cohort=="GV-only") "gv" else c("gv","shr")
  keep <- complete.cases(d[,c(tm,ev,expo_vars,strsplit(covs," \\+ ")[[1]],fe)]) & d[[tm]]>0
  list(d=d[keep,], ev=ev, tm=tm, covs=covs, fe=fe)
}
cov_label <- function(covs) if(covs==reduced) paste0("reduced: ",reduced) else paste0("full: ",full)
MISS_RULE <- "complete-case on outcome, exposure(s), model covariates and frame_extra (bmi+diabetes for GV-only both horizons and SHR-GV 30d); survival time > 0"

# ---- 主模型复现核对 ----
prim <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/02_primary_models_refit/final_primary_models.csv"), stringsAsFactors=FALSE)
repro_rows <- list()
assert_repro <- function(fr, cohort, horizon){
  d <- fr$d; d$gv_z <- safe_z(d$gv)
  if (cohort!="GV-only") d$shr_z <- safe_z(d$shr)
  expo <- if(cohort=="GV-only") "gv_z" else "shr_z + gv_z"
  fit <- coxph(as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",expo," + ",fr$covs)), data=d)
  s <- summary(fit)
  ref <- prim[prim$cohort==cohort & prim$horizon_days==horizon,]
  stopifnot(nrow(ref)==1)
  hr_gv <- s$coefficients["gv_z","exp(coef)"]
  ok_N <- nrow(d)==ref$N; ok_ev <- sum(d[[fr$ev]])==ref$events; ok_hr <- abs(hr_gv-ref$HR_gv_z)<1e-6
  if (cohort!="GV-only") ok_hr <- ok_hr && abs(s$coefficients["shr_z","exp(coef)"]-ref$HR_shr_z)<1e-6
  repro_rows[[length(repro_rows)+1]] <<- data.frame(
    module="alt_gv", cohort=cohort, horizon_days=horizon,
    N=nrow(d), primary_N=ref$N, N_match=ok_N,
    events=sum(d[[fr$ev]]), primary_events=ref$events, events_match=ok_ev,
    HR_gv=hr_gv, primary_HR_gv=ref$HR_gv_z, HR_match=ok_hr)
  if(!(ok_N && ok_ev && ok_hr))
    stop(paste("REPRODUCTION FAILED:", cohort, horizon, "N", nrow(d), "vs", ref$N,
               "events", sum(d[[fr$ev]]), "vs", ref$events, "HR", hr_gv, "vs", ref$HR_gv_z))
  invisible(TRUE)
}

metrics <- data.frame(key=c("sd","cv","arv","tw_arv"),
                      var=c("gv","cv_pct","arv","tw_arv"),
                      label=c("SD (mg/dL)","CV (%)","ARV (mg/dL)","time-weighted ARV (mg/dL)"),
                      min_n=c(2,2,3,3))

run_metric <- function(fr, cohort, horizon, mvar, mlabel, min_n){
  d <- fr$d
  n_frame <- nrow(d)
  d <- d[!is.na(d[[mvar]]) & d$n_glucose_postop_24h>=min_n,]
  if(nrow(d)<100 || sum(d[[fr$ev]])<20)
    return(list(data.frame(cohort=cohort,horizon_days=horizon,metric=mlabel,min_n=min_n,
      N=nrow(d),events=sum(d[[fr$ev]]),frame_N=n_frame,sample_loss_vs_frame=n_frame-nrow(d),
      note="NOT ESTIMABLE (<100 N or <20 events)", stringsAsFactors=FALSE)))
  # 在该指标可计算的同一患者子集上拟合指定指标;paired_role 标记替代指标及其配对 SD 对照
  fit_on_subset <- function(fit_var, fit_label, role){
    dd <- d
    dd$x_z <- safe_z(dd[[fit_var]]); dd$mean_z <- safe_z(dd$glucose_mean_postop_24h)
    if (cohort!="GV-only") dd$shr_z <- safe_z(dd$shr)
    shr_part <- if(cohort!="GV-only") "shr_z + " else ""
    f1 <- as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",shr_part,"x_z + ",fr$covs))
    f2 <- as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",shr_part,"x_z + mean_z + ",fr$covs))
    fit1 <- coxph(f1, data=dd); fit2 <- coxph(f2, data=dd)
    s1 <- summary(fit1); s2 <- summary(fit2)
    # 共线性诊断(线性近似)
    vif_x <- vif_mean <- NA
    v <- tryCatch({
      numvars <- c("x_z","mean_z","age_at_admission","charlson_comorbidity_index","sofa_24h","creat_postop_first","lactate_postop_first")
      if (horizon!=30 && cohort!="GV-only") numvars <- c(numvars,"bmi","wbc_postop_first","albumin_adm_first","hgb_postop_first","platelets_postop_first")
      dq <- dd[, numvars]; dq <- dq[complete.cases(dq),]
      vif_x <- 1/(1-summary(lm(x_z ~ ., data=dq))$r.squared)
      vif_mean <- 1/(1-summary(lm(mean_z ~ . - mean_z + x_z, data=dq))$r.squared)
      c(vif_x, vif_mean)
    }, error=function(e) NULL)
    if(!is.null(v)) { vif_x <- v[1]; vif_mean <- v[2] }
    data.frame(cohort=cohort, horizon_days=horizon, metric=fit_label, min_n=min_n, paired_role=role,
      N=nrow(dd), events=sum(dd[[fr$ev]]), frame_N=n_frame, sample_loss_vs_frame=n_frame-nrow(dd),
      HR_clinical=s1$coefficients["x_z","exp(coef)"], lo1=s1$conf.int["x_z","lower .95"],
      hi1=s1$conf.int["x_z","upper .95"], P_clinical=s1$coefficients["x_z","Pr(>|z|)"],
      HR_plus_mean=s2$coefficients["x_z","exp(coef)"], lo2=s2$conf.int["x_z","lower .95"],
      hi2=s2$conf.int["x_z","upper .95"], P_plus_mean=s2$coefficients["x_z","Pr(>|z|)"],
      mean_HR=s2$coefficients["mean_z","exp(coef)"], mean_P=s2$coefficients["mean_z","Pr(>|z|)"],
      vif_x=vif_x, vif_mean=vif_mean,
      metric_mean=mean(dd[[fit_var]],na.rm=TRUE), metric_sd=sd(dd[[fit_var]],na.rm=TRUE),
      analytic_cohort=paste0(if(cohort=="SHR-GV") "final_v2_A" else if(cohort=="open-core") "final_v2_open_core" else "final_v2_C"," (",cohort,"), primary complete-case frame"),
      covariate_specification=cov_label(fr$covs),
      shr_jointly_adjusted=cohort!="GV-only",
      standardization_sample=paste0("z within this model sample (N=",nrow(dd),")"),
      missing_data_rule=MISS_RULE, note="", stringsAsFactors=FALSE)
  }
  rows <- list(fit_on_subset(mvar, mlabel, "alternative_metric"))
  # 配对 SD 对照:min_n>2 的指标(ARV/tw-ARV)必须在同一子集上同时给出 SD 估计,
  # 使指标间差异不混入样本差异(reviewer 要求)
  if (min_n>2)
    rows[[length(rows)+1]] <- fit_on_subset("gv", paste0("SD (paired comparator on ", mlabel, " subset)"), "paired_SD_comparator")
  rows
}

res <- list()
for (co in c("A","open_core","C")) {
  flag <- paste0("final_v2_", co)
  cname <- c("A"="SHR-GV","open_core"="open-core","C"="GV-only")[co]
  for (h in c(30,365)) {
    fr <- build_primary_frame(dat, flag, cname, h)
    assert_repro(fr, cname, h)   # 第一行 = 主模型复现(失败即停)
    for (i in 1:nrow(metrics)) {
      rr <- tryCatch(run_metric(fr, cname, h, metrics$var[i], metrics$label[i], metrics$min_n[i]),
                     error=function(e) list(data.frame(cohort=cname,horizon_days=h,metric=metrics$label[i],error=conditionMessage(e))))
      for (r in rr) res[[length(res)+1]] <- r
    }
  }
}
res_tab <- dplyr::bind_rows(res)
write_csv(res_tab, file.path(OUT,"alt_gv_models_framefix.csv"))

# 分布表与相关矩阵(GV-only 1y 主 frame 内,描述性)
frC <- build_primary_frame(dat,"final_v2_C","GV-only",365)
dist_rows <- lapply(1:nrow(metrics), function(i){
  v <- frC$d[[metrics$var[i]]]
  data.frame(metric=metrics$label[i], n_valid=sum(!is.na(v)),
    median=median(v,na.rm=TRUE), q1=quantile(v,.25,na.rm=TRUE), q3=quantile(v,.75,na.rm=TRUE),
    min=min(v,na.rm=TRUE), max=max(v,na.rm=TRUE))
})
write_csv(do.call(rbind,dist_rows), file.path(OUT,"alt_gv_metrics_distribution_framefix.csv"))
cm <- cor(frC$d[,c("gv","cv_pct","arv","tw_arv","glucose_mean_postop_24h")], use="pairwise.complete.obs")
write.csv(round(cm,3), file.path(OUT,"alt_gv_correlation_matrix_framefix.csv"))

# ARV 样本损失(相对各自主 frame)
arv_loss <- dplyr::bind_rows(lapply(c("A","open_core","C"), function(co){
  cname <- c("A"="SHR-GV","open_core"="open-core","C"="GV-only")[co]
  dplyr::bind_rows(lapply(c(30,365), function(h){
    fr <- build_primary_frame(dat, paste0("final_v2_",co), cname, h)
    data.frame(cohort=cname, horizon_days=h, frame_N=nrow(fr$d),
      n_arv_available=sum(!is.na(fr$d$arv) & fr$d$n_glucose_postop_24h>=3),
      loss=nrow(fr$d)-sum(!is.na(fr$d$arv) & fr$d$n_glucose_postop_24h>=3))
  }))
}))
write_csv(arv_loss, file.path(OUT,"alt_gv_arv_sample_loss_framefix.csv"))
write_csv(dplyr::bind_rows(repro_rows), file.path(OUT,"reproduction_check_alt_gv.csv"))
print(dplyr::bind_rows(repro_rows))
print(res_tab %>% filter(cohort=="GV-only") %>% select(cohort,horizon_days,metric,N,events,HR_clinical,P_clinical,HR_plus_mean,P_plus_mean))
cat("PART2_FRAMEFIX_DONE\n")
sink()
