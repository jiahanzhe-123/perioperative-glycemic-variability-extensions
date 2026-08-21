# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 01_refit_primary_models: 用 final blood-only 暴露重跑 6 个核心模型(规则与冻结完全一致)
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr)})
ROOT <- PGV("mimic_work")
OUT <- file.path(ROOT,"final_statistical_freeze_20260727/02_primary_models_refit")
sink(file.path(ROOT,"final_statistical_freeze_20260727/logs/01_primary_refit.log"), split=TRUE)

fz <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/06_release_objects/final_analysis_dataset.csv"), stringsAsFactors=FALSE)
num <- c("age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first","creat_postop_first",
         "wbc_postop_first","albumin_adm_first","hgb_postop_first","platelets_postop_first","sofa_24h",
         "survival_time_days","postop_30d_death_flag","one_year_death_flag","diabetes",
         "gv_bloodonly","shr_bloodonly","mean_glucose_bloodonly")
for (nm in num) fz[[nm]] <- suppressWarnings(as.numeric(as.character(fz[[nm]])))
fz$gender <- factor(fz$gender)
pg <- tolower(trimws(fz$procedure_group_main))
fz$procedure_group_main_model <- factor(ifelse(pg=="cabg","cabg",ifelse(pg=="open_valve","open_valve",
  ifelse(pg=="transplant_vad","transplant_vad","aortic_congenital_other"))),
  levels=c("cabg","open_valve","aortic_congenital_other","transplant_vad"))
fz$event_30d <- as.integer(fz$postop_30d_death_flag==1)
fz$event_1y <- as.integer(fz$one_year_death_flag==1)
fz$t30 <- pmin(fz$survival_time_days,30); fz$t365 <- pmin(fz$survival_time_days,365)
reduced <- "age_at_admission + gender + charlson_comorbidity_index + procedure_group_main_model + lactate_postop_first + creat_postop_first + sofa_24h"
full <- paste(reduced, "+ bmi + wbc_postop_first + albumin_adm_first + hgb_postop_first + platelets_postop_first")
safe_z <- function(x) as.numeric(scale(suppressWarnings(as.numeric(x))))

run_model <- function(cohort_flag, cohort, horizon){
  d <- fz[fz[[cohort_flag]]==1,]
  ev <- if(horizon==30) "event_30d" else "event_1y"; tm <- if(horizon==30) "t30" else "t365"
  covs <- if(cohort=="GV-only") reduced else if(horizon==30) reduced else full
  fe <- if(cohort=="GV-only" || (cohort=="SHR-GV" && horizon==30)) c("bmi","diabetes") else NULL
  expo_vars <- if(cohort=="GV-only") "gv_bloodonly" else c("gv_bloodonly","shr_bloodonly")
  keep <- complete.cases(d[,c(tm,ev,expo_vars,strsplit(covs," \\+ ")[[1]],fe)]) & d[[tm]]>0
  d <- d[keep,]; d$gv_z <- safe_z(d$gv_bloodonly); d$shr_z <- safe_z(d$shr_bloodonly)
  expo <- if(cohort=="GV-only") "gv_z" else "shr_z + gv_z"
  f <- as.formula(paste0("Surv(",tm,",",ev,") ~ ",expo," + ",covs))
  fit <- coxph(f, data=d)
  s <- summary(fit); z <- cox.zph(fit)
  out <- data.frame(cohort=cohort, horizon_days=horizon, N=nrow(d), events=sum(d[[ev]]),
    exposure_mean=mean(d$gv_bloodonly,na.rm=TRUE), exposure_sd=sd(d$gv_bloodonly,na.rm=TRUE),
    model_type=if(cohort=="GV-only") "GV-only Cox" else "joint SHR+GV Cox",
    covariate_specification=if(horizon==30) "reduced" else "full",
    ph_global_p=z$table["GLOBAL","p"], stringsAsFactors=FALSE)
  for (e in if(cohort=="GV-only") "gv_z" else c("shr_z","gv_z")) {
    out[[paste0("HR_",e)]] <- s$coefficients[e,"exp(coef)"]
    out[[paste0("lo_",e)]] <- s$conf.int[e,"lower .95"]
    out[[paste0("hi_",e)]] <- s$conf.int[e,"upper .95"]
    out[[paste0("P_",e)]] <- s$coefficients[e,"Pr(>|z|)"]
    out[[paste0("ph_p_",e)]] <- z$table[e,"p"]
  }
  attr(out,"fit") <- NULL
  list(tab=out, fit=fit)
}
res <- list(); fits <- list()
for (co in c("final_v2_A","final_v2_open_core","final_v2_C")) {
  cn <- c(final_v2_A="SHR-GV", final_v2_open_core="open-core", final_v2_C="GV-only")[co]
  for (h in c(30,365)) {
    r <- run_model(co, cn, h)
    res[[length(res)+1]] <- r$tab
    fits[[paste(cn,h,sep="_")]] <- r$fit
  }
}
tab <- dplyr::bind_rows(res)
write_csv(tab, file.path(OUT,"final_primary_models.csv"))
saveRDS(fits, file.path(OUT,"final_primary_model_objects.rds"))
ph <- tab %>% select(cohort,horizon_days,ph_global_p,starts_with("ph_p_"))
write_csv(ph, file.path(OUT,"final_primary_ph_tests.csv"))
write_csv(tab %>% select(cohort,horizon_days,N,events,model_type,covariate_specification), file.path(OUT,"final_primary_model_diagnostics.csv"))
print(tab %>% select(cohort,horizon_days,N,events,starts_with("HR_"),starts_with("P_")) %>%
        mutate(across(where(is.numeric), ~round(.,4))))

# 与冻结结果比较
frozen_ref <- data.frame(
  cohort=c("SHR-GV","SHR-GV","open-core","open-core","GV-only","GV-only"),
  horizon_days=c(30,365,30,365,30,365),
  frozen_HR=c(1.090,1.122,1.278,1.152,1.074,1.070),
  frozen_lo=c(0.978,1.027,1.108,1.038,1.001,1.013),
  frozen_hi=c(1.216,1.225,1.475,1.279,1.153,1.129),
  frozen_P=c(0.119,0.011,0.000,0.008,0.046,0.015),
  frozen_N=c(4869,4381,4757,4162,9787,9787), frozen_events=c(114,279,99,254,256,650))
cmp <- merge(tab %>% mutate(HR=coalesce(HR_gv_z), P=coalesce(P_gv_z), lo=coalesce(lo_gv_z), hi=coalesce(hi_gv_z)),
             frozen_ref, by=c("cohort","horizon_days"))
cmp$abs_diff <- round(cmp$HR - cmp$frozen_HR, 4)
cmp$rel_diff_pct <- round(100*(cmp$HR - cmp$frozen_HR)/cmp$frozen_HR, 2)
cmp$conclusion_changed <- (cmp$P<0.05) != (cmp$frozen_P<0.05)
cmp$interpretation_changed <- cmp$conclusion_changed
write_csv(cmp, file.path(ROOT,"final_statistical_freeze_20260727/05_comparison_tables/frozen_vs_bloodonly_primary_models.csv"))
print(cmp %>% select(cohort,horizon_days,frozen_HR,HR,frozen_P,P,abs_diff,rel_diff_pct,conclusion_changed))
cat("PRIMARY_REFIT_DONE\n")
sink()
