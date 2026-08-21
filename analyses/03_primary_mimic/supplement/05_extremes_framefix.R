# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 05_extremes_bloodonly_framefix.R
# 修复 reviewer 问题①:低血糖/高血糖/血糖负荷分析
# 修复点:
#   1. 完整病例 frame 与主模型完全一致(补回 frame_extra:bmi+diabetes)。
#   2. 非 GV-only 队列 keep 中显式纳入 shr(口径修正,数值不变)。
#   3. 【数据源修复】负荷变量(any_lt70 等)改自 final_glucose_features.csv——
#      原脚本误用 00_glucose_series/series_features.csv,其测量集合与主 GV 不一致
#      (6415/10581 个 stay 的 SD 不同,max diff 46.7 mg/dL),违反"所有指标与主 GV
#      同一观测集合"的预设规则。
#   4. adjustment="none" 行即主模型,与 final_primary_models.csv 自动核对,不一致即 stop。
#   5. 每行输出元数据:精确队列/协变量规格/SHR是否同调/标准化样本/缺失规则。
# 输入与暴露版本不变(final blood-only sequence),不覆盖原脚本与原输出。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr)})
ROOT <- PGV("mimic_work")
OUT  <- file.path(ROOT,"reviewer_issue_verification_20260728/01_frame_repair/extremes")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sink(file.path(ROOT,"reviewer_issue_verification_20260728/logs/05_extremes_framefix.log"), split=TRUE)

dat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA","NULL"))
# 修复:负荷变量必须来自与主 GV 同一测量集合的 final_glucose_features.csv
feat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_glucose_features.csv"), stringsAsFactors=FALSE)
dat <- merge(dat, feat[,c("stay_id","any_lt70","any_lt54","any_gt180","any_gt250","prop_lt70","prop_gt180","prop_70_180")],
             by="stay_id", all.x=TRUE)
num <- c("shr","gv","age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first",
         "creat_postop_first","wbc_postop_first","albumin_adm_first","hgb_postop_first",
         "platelets_postop_first","sofa_24h","survival_time_days","postop_30d_death_flag","one_year_death_flag",
         "diabetes","glucose_mean_postop_24h","glucose_min_postop_24h","glucose_max_postop_24h","glucose_range_postop_24h",
         "prop_lt70","prop_gt180","prop_70_180","n_glucose_postop_24h")
for (nm in num) dat[[nm]] <- suppressWarnings(as.numeric(as.character(dat[[nm]])))
for (nm in c("any_lt70","any_lt54","any_gt180","any_gt250"))
  dat[[nm]] <- as.integer(dat[[nm]] %in% c(TRUE,"TRUE","True","true",1,"1"))
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
cohort_flag_of <- function(cohort) c("SHR-GV"="final_v2_A","open-core"="final_v2_open_core","GV-only"="final_v2_C")[cohort]

prim <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/02_primary_models_refit/final_primary_models.csv"), stringsAsFactors=FALSE)
repro_rows <- list()
assert_repro <- function(fr, cohort, horizon, module){
  d <- fr$d; d$gv_z <- safe_z(d$gv)
  if (cohort!="GV-only") d$shr_z <- safe_z(d$shr)
  expo <- if(cohort=="GV-only") "gv_z" else "shr_z + gv_z"
  fit <- coxph(as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",expo," + ",fr$covs)), data=d)
  s <- summary(fit)
  ref <- prim[prim$cohort==cohort & prim$horizon_days==horizon,]; stopifnot(nrow(ref)==1)
  hr_gv <- s$coefficients["gv_z","exp(coef)"]
  ok_N <- nrow(d)==ref$N; ok_ev <- sum(d[[fr$ev]])==ref$events; ok_hr <- abs(hr_gv-ref$HR_gv_z)<1e-6
  if (cohort!="GV-only") ok_hr <- ok_hr && abs(s$coefficients["shr_z","exp(coef)"]-ref$HR_shr_z)<1e-6
  repro_rows[[length(repro_rows)+1]] <<- data.frame(
    module=module, cohort=cohort, horizon_days=horizon,
    N=nrow(d), primary_N=ref$N, N_match=ok_N,
    events=sum(d[[fr$ev]]), primary_events=ref$events, events_match=ok_ev,
    HR_gv=hr_gv, primary_HR_gv=ref$HR_gv_z, HR_match=ok_hr)
  if(!(ok_N && ok_ev && ok_hr))
    stop(paste("REPRODUCTION FAILED:", module, cohort, horizon))
  invisible(TRUE)
}

# ---- 1. 描述:GV 三分位中的负荷分布(GV-only 1y 主 frame 内) ----
frC1y <- build_primary_frame(dat,"final_v2_C","GV-only",365)
dC <- frC1y$d
dC$gv_tertile <- cut(dC$gv, breaks=quantile(dC$gv, c(0,1/3,2/3,1), na.rm=TRUE), labels=c("T1","T2","T3"), include.lowest=TRUE)
desc <- dC %>% group_by(gv_tertile) %>%
  summarise(n=n(),
    any_lt70_pct=round(100*mean(any_lt70,na.rm=TRUE),1),
    any_lt54_pct=round(100*mean(any_lt54,na.rm=TRUE),2),
    any_gt180_pct=round(100*mean(any_gt180,na.rm=TRUE),1),
    any_gt250_pct=round(100*mean(any_gt250,na.rm=TRUE),1),
    prop_lt70_med=round(median(prop_lt70,na.rm=TRUE),3),
    prop_gt180_med=round(median(prop_gt180,na.rm=TRUE),3),
    prop_70_180_med=round(median(prop_70_180,na.rm=TRUE),3),
    min_glucose_med=round(median(glucose_min_postop_24h,na.rm=TRUE),0),
    max_glucose_med=round(median(glucose_max_postop_24h,na.rm=TRUE),0),
    range_med=round(median(glucose_range_postop_24h,na.rm=TRUE),0), .groups="drop")
write_csv(desc, file.path(OUT,"extremes_burden_by_gv_tertile_framefix.csv"))
print(desc)

cormat <- cor(dC[,c("gv","any_lt70","any_gt180","prop_lt70","prop_gt180","glucose_min_postop_24h","glucose_max_postop_24h","glucose_range_postop_24h","glucose_mean_postop_24h")], use="pairwise.complete.obs")
write.csv(round(cormat,3), file.path(OUT,"extremes_extreme_variable_correlations_framefix.csv"))

# ---- 2. 序列调整(none 行 = 主模型复现) ----
run_ext <- function(fr, cohort, horizon, adj){
  d <- fr$d
  extra <- switch(adj,
    none="", hypo=" + any_lt70", hyper=" + any_gt180",
    minmax=" + glucose_min_postop_24h + glucose_max_postop_24h",
    burden=" + prop_lt70 + prop_gt180")
  if (adj=="hypo" && sum(d$any_lt70==1)<10)
    return(data.frame(cohort=cohort,horizon_days=horizon,adjustment=adj,N=nrow(d),events=sum(d[[fr$ev]]),
      note=paste0("低血糖事件过少(n=",sum(d$any_lt70==1),"),仅描述"), stringsAsFactors=FALSE))
  d$gv_z <- safe_z(d$gv)
  expo <- if(cohort=="GV-only") "gv_z" else "shr_z + gv_z"
  if (cohort!="GV-only") d$shr_z <- safe_z(d$shr)
  fit <- coxph(as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",expo," + ",fr$covs,extra)), data=d)
  s <- summary(fit)
  data.frame(cohort=cohort, horizon_days=horizon, adjustment=adj, N=nrow(d), events=sum(d[[fr$ev]]),
    HR_gv=s$coefficients["gv_z","exp(coef)"], lo=s$conf.int["gv_z","lower .95"],
    hi=s$conf.int["gv_z","upper .95"], P_gv=s$coefficients["gv_z","Pr(>|z|)"],
    gv_mean=mean(d$gv,na.rm=TRUE), gv_sd=sd(d$gv,na.rm=TRUE),
    analytic_cohort=paste0(cohort_flag_of(cohort)," (",cohort,"), primary complete-case frame"),
    covariate_specification=cov_label(fr$covs),
    shr_jointly_adjusted=cohort!="GV-only",
    standardization_sample=paste0("z within this model sample (N=",nrow(d),")"),
    missing_data_rule=MISS_RULE, note="", stringsAsFactors=FALSE)
}
ext_list <- list()
for (cf in c("final_v2_A","final_v2_C")) {
  cname <- c(final_v2_A="SHR-GV",final_v2_C="GV-only")[cf]
  for (h in c(30,365)) {
    fr <- build_primary_frame(dat, cf, cname, h)
    assert_repro(fr, cname, h, "extremes")
    for (a in c("none","hypo","hyper","minmax","burden"))
      ext_list[[length(ext_list)+1]] <- run_ext(fr, cname, h, a)
  }
}
ext_tab <- dplyr::bind_rows(ext_list)
write_csv(ext_tab, file.path(OUT,"extremes_extreme_adjustments_framefix.csv"))
write_csv(dplyr::bind_rows(repro_rows), file.path(OUT,"reproduction_check_extremes.csv"))
print(dplyr::bind_rows(repro_rows))
print(ext_tab %>% select(cohort,horizon_days,adjustment,N,events,HR_gv,P_gv,note))
cat("PART5_FRAMEFIX_DONE\n")
sink()
