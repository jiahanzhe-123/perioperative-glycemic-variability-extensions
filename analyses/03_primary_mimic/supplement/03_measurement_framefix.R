# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 03_measurement_bloodonly_framefix.R
# 修复 reviewer 问题①:测量过程敏感性分析
# 修复点:
#   1. 完整病例 frame 与主模型完全一致(补回 frame_extra:bmi+diabetes)。
#   2. 非 GV-only 队列 keep 中显式纳入 shr(原脚本依赖 coxph 静默删行;
#      实测 open-core/SHR-GV 池 shr 均完整,数值不变,仅为口径正确)。
#   3. adjustment="none" 行即主模型,与 final_primary_models.csv 自动核对,不一致即 stop。
#   4. 每行输出元数据:精确队列/协变量规格/SHR是否同调/标准化样本/缺失规则。
# 输入与暴露版本不变(final blood-only sequence),不覆盖原脚本与原输出。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr); library(splines)})
ROOT <- PGV("mimic_work")
OUT  <- file.path(ROOT,"reviewer_issue_verification_20260728/01_frame_repair/measurement")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sink(file.path(ROOT,"reviewer_issue_verification_20260728/logs/03_measurement_framefix.log"), split=TRUE)

dat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA","NULL"))
feat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_glucose_features.csv"), stringsAsFactors=FALSE)
dat <- merge(dat, feat[, setdiff(names(feat), c("arv","tw_arv","span_hours"))], by="stay_id", all.x=TRUE)
num <- c("shr","gv","age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first",
         "creat_postop_first","wbc_postop_first","albumin_adm_first","hgb_postop_first",
         "platelets_postop_first","sofa_24h","survival_time_days","postop_30d_death_flag","one_year_death_flag",
         "diabetes","glucose_mean_postop_24h","n_glucose_postop_24h","span_hours","density_per_hour",
         "frac_poct","frac_central_lab","frac_blood_gas","frac_icu_charted")
for (nm in num) dat[[nm]] <- suppressWarnings(as.numeric(as.character(dat[[nm]])))
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
meta_cols <- function(d, cohort, covs) list(
  gv_mean=mean(d$gv,na.rm=TRUE), gv_sd=sd(d$gv,na.rm=TRUE),
  analytic_cohort=paste0(cohort_flag_of(cohort)," (",cohort,"), primary complete-case frame"),
  covariate_specification=cov_label(covs),
  shr_jointly_adjusted=cohort!="GV-only",
  standardization_sample=paste0("z within this model sample (N=",nrow(d),")"),
  missing_data_rule=MISS_RULE)

# ---- 1. 描述:GV 与测量过程的相关(GV-only 1y 主 frame 内) ----
frC1y <- build_primary_frame(dat,"final_v2_C","GV-only",365)
d_c <- frC1y$d
stopifnot(sum(is.na(d_c$span_hours))==0, sum(is.na(d_c$frac_poct))==0)  # frame 内测量过程变量应完整
desc <- data.frame(
  variable=c("glucose count","span (hours)","density (/h)","POCT fraction","central-lab fraction","blood-gas fraction"),
  cor_with_gv=c(cor(d_c$gv,d_c$n_glucose_postop_24h), cor(d_c$gv,d_c$span_hours),
                cor(d_c$gv,d_c$density_per_hour), cor(d_c$gv,d_c$frac_poct),
                cor(d_c$gv,d_c$frac_central_lab), cor(d_c$gv,d_c$frac_blood_gas)),
  cor_with_mean=c(cor(d_c$glucose_mean_postop_24h,d_c$n_glucose_postop_24h), cor(d_c$glucose_mean_postop_24h,d_c$span_hours),
                cor(d_c$glucose_mean_postop_24h,d_c$density_per_hour), cor(d_c$glucose_mean_postop_24h,d_c$frac_poct),
                cor(d_c$glucose_mean_postop_24h,d_c$frac_central_lab), cor(d_c$glucose_mean_postop_24h,d_c$frac_blood_gas))
)
write_csv(desc, file.path(OUT,"measurement_process_correlations_framefix.csv"))

# ---- 2. 序列调整(none 行 = 主模型复现) ----
run_adj <- function(fr, cohort, horizon, adj){
  d <- fr$d
  extra <- switch(adj,
    none="", count=" + n_glucose_postop_24h", span=" + span_hours",
    count_span=" + n_glucose_postop_24h + span_hours",
    count_span_source=" + n_glucose_postop_24h + span_hours + frac_poct")
  d$gv_z <- safe_z(d$gv)
  expo <- if(cohort=="GV-only") "gv_z" else "shr_z + gv_z"
  if (cohort!="GV-only") d$shr_z <- safe_z(d$shr)
  fit <- coxph(as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",expo," + ",fr$covs,extra)), data=d)
  s <- summary(fit)
  data.frame(cohort=cohort, horizon_days=horizon, adjustment=adj, N=nrow(d), events=sum(d[[fr$ev]]),
    HR_gv=s$coefficients["gv_z","exp(coef)"], lo=s$conf.int["gv_z","lower .95"],
    hi=s$conf.int["gv_z","upper .95"], P=s$coefficients["gv_z","Pr(>|z|)"],
    meta_cols(d, cohort, fr$covs), stringsAsFactors=FALSE)
}
adj_list <- list()
for (cf in c("final_v2_A","final_v2_open_core","final_v2_C")) {
  cname <- c(final_v2_A="SHR-GV",final_v2_open_core="open-core",final_v2_C="GV-only")[cf]
  for (h in c(30,365)) {
    fr <- build_primary_frame(dat, cf, cname, h)
    assert_repro(fr, cname, h, "measurement")
    for (a in c("none","count","span","count_span","count_span_source"))
      adj_list[[length(adj_list)+1]] <- run_adj(fr, cname, h, a)
  }
}
adj_tab <- dplyr::bind_rows(adj_list)
write_csv(adj_tab, file.path(OUT,"measurement_sequential_adjustment_framefix.csv"))

# ---- 3. 非线性检查(count/span RCS, GV-only 1y 主 frame 内) ----
d <- frC1y$d; d$gv_z <- safe_z(d$gv)
nl <- list()
for (v in c("n_glucose_postop_24h","span_hours")) {
  f_lin <- as.formula(paste0("Surv(t365,event_1y) ~ gv_z + ", v, " + ", reduced))
  f_rcs <- as.formula(paste0("Surv(t365,event_1y) ~ gv_z + ns(", v, ", df=3) + ", reduced))
  fit1 <- coxph(f_lin, data=d); fit2 <- coxph(f_rcs, data=d)
  p <- anova(fit1,fit2)$`Pr(>|Chi|)`[2]
  nl[[v]] <- data.frame(variable=v, N=nrow(d), events=sum(d$event_1y),
    linear_term_p=summary(fit1)$coefficients[v,"Pr(>|z|)"], rcs_vs_linear_p=p,
    analytic_cohort="final_v2_C (GV-only) 1y primary frame", covariate_specification=cov_label(reduced),
    shr_jointly_adjusted=FALSE, missing_data_rule=MISS_RULE)
}
write_csv(dplyr::bind_rows(nl), file.path(OUT,"measurement_nonlinearity_checks_framefix.csv"))

# ---- 4. 限制性分析(自主 frame 起,样本损失透明) ----
series <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/glucose_minute_series_bloodonly.csv"), stringsAsFactors=FALSE)
series$minute <- as.POSIXct(series$minute)
poct_gv <- series %>% filter(grepl("poct", sources)) %>% group_by(stay_id) %>%
  summarise(gv_poct=if(n()>=2) sd(value,na.rm=TRUE) else NA_real_, .groups="drop")
lab_gv <- series %>% filter(grepl("central_lab", sources)) %>% group_by(stay_id) %>%
  summarise(gv_lab=if(n()>=2) sd(value,na.rm=TRUE) else NA_real_, .groups="drop")
dat <- dat %>% left_join(poct_gv, by="stay_id") %>% left_join(lab_gv, by="stay_id")

run_restrict <- function(cohort_flag, cohort, horizon, label, mask, gvcol="gv"){
  fr <- build_primary_frame(dat, cohort_flag, cohort, horizon)
  d <- fr$d[fr$d$stay_id %in% dat$stay_id[dat[[cohort_flag]]==1 & mask],]
  n_frame <- nrow(fr$d)
  if(nrow(d)<100 || sum(d[[fr$ev]])<20)
    return(data.frame(cohort=cohort,horizon_days=horizon,restriction=label,N=nrow(d),events=sum(d[[fr$ev]]),
      frame_N=n_frame,sample_loss_vs_frame=n_frame-nrow(d),
      note="NOT ESTIMABLE (<100 N or <20 events)", stringsAsFactors=FALSE))
  d$gv_z <- safe_z(d[[gvcol]])
  expo <- if(cohort=="GV-only") "gv_z" else "shr_z + gv_z"
  if (cohort!="GV-only") d$shr_z <- safe_z(d$shr)
  fit <- coxph(as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",expo," + ",fr$covs)), data=d)
  s <- summary(fit)
  data.frame(cohort=cohort, horizon_days=horizon, restriction=label, N=nrow(d), events=sum(d[[fr$ev]]),
    frame_N=n_frame, sample_loss_vs_frame=n_frame-nrow(d),
    HR_gv=s$coefficients["gv_z","exp(coef)"], lo=s$conf.int["gv_z","lower .95"],
    hi=s$conf.int["gv_z","upper .95"], P=s$coefficients["gv_z","Pr(>|z|)"],
    gv_mean=mean(d[[gvcol]],na.rm=TRUE), gv_sd=sd(d[[gvcol]],na.rm=TRUE),
    analytic_cohort=paste0(cohort_flag_of(cohort)," (",cohort,"), primary frame + restriction"),
    covariate_specification=cov_label(fr$covs), shr_jointly_adjusted=cohort!="GV-only",
    standardization_sample=paste0("z within restricted sample (N=",nrow(d),")"),
    missing_data_rule=MISS_RULE, note="", stringsAsFactors=FALSE)
}
rests <- list(
  list(">=3 glucose", dat$n_glucose_postop_24h>=3, "gv"),
  list(">=4 glucose", dat$n_glucose_postop_24h>=4, "gv"),
  list("span>=12h", dat$span_hours>=12 & !is.na(dat$span_hours), "gv"),
  list("POCT-dominant (>=50% POCT)", !is.na(dat$frac_poct) & dat$frac_poct>=0.5, "gv"),
  list("POCT-only GV (>=2 POCT)", !is.na(dat$gv_poct), "gv_poct"),
  list("central-lab-only GV (>=2 lab)", !is.na(dat$gv_lab), "gv_lab"))
rest_tab <- dplyr::bind_rows(lapply(rests, function(r)
  dplyr::bind_rows(lapply(c("final_v2_A","final_v2_C"), function(cf)
    dplyr::bind_rows(lapply(c(30,365), function(h)
      tryCatch(run_restrict(cf, c(final_v2_A="SHR-GV",final_v2_C="GV-only")[cf], h, r[[1]], r[[2]], r[[3]]),
               error=function(e) data.frame(cohort=c(final_v2_A="SHR-GV",final_v2_C="GV-only")[cf],horizon_days=h,restriction=r[[1]],error=conditionMessage(e)))))))))
write_csv(rest_tab, file.path(OUT,"measurement_restriction_analyses_framefix.csv"))
write_csv(dplyr::bind_rows(repro_rows), file.path(OUT,"reproduction_check_measurement.csv"))
print(dplyr::bind_rows(repro_rows))
print(adj_tab %>% filter(horizon_days==365) %>% select(cohort,adjustment,N,events,HR_gv,P))
cat("PART3_FRAMEFIX_DONE\n")
sink()
