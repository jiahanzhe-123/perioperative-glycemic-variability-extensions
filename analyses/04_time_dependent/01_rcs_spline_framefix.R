# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 06_rcs_spline_framefix.R
# 修复 reviewer 补查项:Supplementary Table S2(RCS)仍使用宽 frame。
# 根因:02_run_final_v3_bloodonly.R 中 fit_rcs_one() 的 frame 只要求
#   exposure+协变量+sofa 完整,未使用 fit_primary() 的 main_sequence common frame
#   (额外要求 bmi+diabetes) → S2 出现 N=5013/10325,与主模型 4869/9786 不一致。
# 本脚本:
#   1. 用与主模型完全一致的 primary frame 重跑全部 RCS;
#   2. 线性版拟合必须精确复现主模型 HR(逐模型断言,失败即停);
#   3. 用旧(宽)frame 复制一遍,验证能复现当前 S2 数字(证明 frame 是唯一差异源);
#   4. 每行输出完整元数据。
# 暴露版本不变(final blood-only);冻结的模型规格选择不变(30d reduced / 1y full);
# knots 在各自分析样本内重取(0.05/0.35/0.65/0.95),与主分析"样本内标准化"口径一致。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr); library(rms)})
ROOT <- PGV("mimic_work")
OUT  <- file.path(ROOT,"reviewer_issue_verification_20260728/03_rcs_framefix")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sink(file.path(ROOT,"reviewer_issue_verification_20260728/logs/06_rcs_framefix.log"), split=TRUE)

dat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA","NULL"))
num <- c("shr","gv","age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first",
         "creat_postop_first","wbc_postop_first","albumin_adm_first","hgb_postop_first",
         "platelets_postop_first","sofa_24h","survival_time_days","postop_30d_death_flag","one_year_death_flag",
         "diabetes","glucose_mean_postop_24h")
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
# 旧(宽)frame:仅 exposure+协变量+sofa 完整 —— 复刻 fit_rcs_one 原行为,仅用于对照复现
build_legacy_rcs_frame <- function(dat, cohort_flag, cohort, horizon){
  d <- dat[dat[[cohort_flag]]==1,]
  ev <- if(horizon==30) "event_30d" else "event_1y"; tm <- if(horizon==30) "t30" else "t365"
  covs <- if(cohort=="GV-only") reduced else if(horizon==30) reduced else full
  expo_vars <- if(cohort=="GV-only") "gv" else c("gv","shr")
  keep <- complete.cases(d[,c(tm,ev,expo_vars,strsplit(covs," \\+ ")[[1]])]) & d[[tm]]>0
  list(d=d[keep,], ev=ev, tm=tm, covs=covs, fe=NULL)
}
cov_label <- function(covs) if(covs==reduced) "M3_reduced" else "M3_full"
MISS_RULE <- "complete-case on outcome, exposure(s), model covariates and frame_extra (bmi+diabetes for GV-only both horizons and SHR-GV 30d); survival time > 0"

lrt_p <- function(a, b){
  la <- tryCatch(logLik(a), error=function(e) NULL); lb <- tryCatch(logLik(b), error=function(e) NULL)
  if (is.null(la) || is.null(lb)) return(NA_real_)
  ddf <- attr(lb,"df") - attr(la,"df")
  if (!is.finite(ddf) || ddf<=0) return(NA_real_)
  pchisq(2*(as.numeric(lb)-as.numeric(la)), df=ddf, lower.tail=FALSE)
}

prim <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/02_primary_models_refit/final_primary_models.csv"), stringsAsFactors=FALSE)
repro_rows <- list(); legacy_rows <- list(); rcs_rows <- list()

run_rcs <- function(fr, cohort, horizon, exposure, frame_label, check_primary=FALSE){
  d <- fr$d
  other <- if(cohort=="GV-only") character(0) else setdiff(c("shr","gv"), exposure)
  d[[paste0(exposure,"_z")]] <- safe_z(d[[exposure]])
  for (o in other) d[[paste0(o,"_z")]] <- safe_z(d[[o]])
  ez <- paste0(exposure,"_z"); oz <- if(length(other)) paste0(other,"_z") else NULL
  knots <- as.numeric(quantile(d[[ez]], probs=c(.05,.35,.65,.95), names=FALSE, type=7))
  terms_base <- c(oz, strsplit(fr$covs," \\+ ")[[1]])
  knot_text <- paste(format(knots, scientific=FALSE, digits=12), collapse=",")
  spline_fml <- as.formula(paste("Surv(",fr$tm,",",fr$ev,") ~ rcs(",ez,", c(",knot_text,")) +", paste(terms_base, collapse=" + ")))
  linear_fml <- as.formula(paste("Surv(",fr$tm,",",fr$ev,") ~ ",ez," + ", paste(terms_base, collapse=" + ")))
  noexp_fml  <- as.formula(paste("Surv(",fr$tm,",",fr$ev,") ~ ", paste(terms_base, collapse=" + ")))
  spline_fit <- coxph(spline_fml, data=d, x=TRUE, y=TRUE, model=TRUE, ties="efron")
  linear_fit <- coxph(linear_fml, data=d, x=TRUE, y=TRUE, model=TRUE, ties="efron")
  noexp_fit  <- coxph(noexp_fml,  data=d, x=TRUE, y=TRUE, model=TRUE, ties="efron")
  # 线性版 = 主模型:必须精确复现主 HR
  if (check_primary) {
    ref <- prim[prim$cohort==cohort & prim$horizon_days==horizon,]; stopifnot(nrow(ref)==1)
    hr_lin <- summary(linear_fit)$coefficients[ez,"exp(coef)"]
    ref_hr <- if(exposure=="gv") ref$HR_gv_z else ref$HR_shr_z
    ok_N <- nrow(d)==ref$N; ok_ev <- sum(d[[fr$ev]])==ref$events; ok_hr <- abs(hr_lin-ref_hr)<1e-6
    repro_rows[[length(repro_rows)+1]] <<- data.frame(
      module="rcs", cohort=cohort, horizon_days=horizon, exposure=exposure,
      N=nrow(d), primary_N=ref$N, N_match=ok_N,
      events=sum(d[[fr$ev]]), primary_events=ref$events, events_match=ok_ev,
      HR_linear=hr_lin, primary_HR=ref_hr, HR_match=ok_hr)
    if(!(ok_N && ok_ev && ok_hr))
      stop(paste("REPRODUCTION FAILED: rcs", cohort, horizon, exposure))
  }
  data.frame(cohort=cohort, horizon_days=horizon,
    endpoint=if(horizon==30) "30-day mortality" else "1-year mortality",
    exposure=if(exposure=="shr") "SHR" else "GV",
    N=nrow(d), events=sum(d[[fr$ev]]),
    knot_probabilities="0.05;0.35;0.65;0.95",
    knot_values=paste(format(knots, digits=10), collapse=";"),
    reference_median=median(d[[ez]], na.rm=TRUE),
    overall_p=lrt_p(noexp_fit, spline_fit),
    nonlinear_p=lrt_p(linear_fit, spline_fit),
    implementation="rms::rcs basis + survival::coxph",
    model_specification=cov_label(fr$covs),
    analytic_cohort=paste0(if(cohort=="SHR-GV") "final_v2_A" else if(cohort=="open-core") "final_v2_open_core" else "final_v2_C"," (",cohort,"), primary complete-case frame"),
    covariate_specification=paste0(cov_label(fr$covs),": ",fr$covs),
    shr_jointly_adjusted=cohort!="GV-only",
    standardization_sample=paste0("z within this analytic sample (N=",nrow(d),"); knots within same sample"),
    missing_data_rule=if(frame_label=="primary") MISS_RULE else "LEGACY wide frame (exposure+covariates+sofa complete only) — reproduction check only",
    frame=frame_label, stringsAsFactors=FALSE)
}

for (cf in c("final_v2_A","final_v2_open_core","final_v2_C")) {
  cohort <- c(final_v2_A="SHR-GV",final_v2_open_core="open-core",final_v2_C="GV-only")[cf]
  exposures <- if(cohort=="GV-only") "gv" else c("shr","gv")
  for (h in c(30,365)) {
    fr_new <- build_primary_frame(dat, cf, cohort, h)
    fr_old <- build_legacy_rcs_frame(dat, cf, cohort, h)
    for (ex in exposures) {
      rcs_rows[[length(rcs_rows)+1]] <- run_rcs(fr_new, cohort, h, ex, "primary", check_primary=TRUE)
      legacy_rows[[length(legacy_rows)+1]] <- run_rcs(fr_old, cohort, h, ex, "legacy", check_primary=FALSE)
    }
  }
}
new_tab <- dplyr::bind_rows(rcs_rows)
old_tab <- dplyr::bind_rows(legacy_rows)
chk <- dplyr::bind_rows(repro_rows)
write_csv(new_tab, file.path(OUT,"rcs_primary_frame_framefix.csv"))
write_csv(old_tab, file.path(OUT,"rcs_legacy_frame_replication.csv"))
write_csv(chk, file.path(OUT,"reproduction_check_rcs.csv"))

# 与当前 S2(freeze v3 RCS 表)对照
s2 <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/03_dependent_analyses_refit/v3_pipeline/tables/final_v3_rcs.csv"), stringsAsFactors=FALSE)
s2 <- s2[s2$cohort %in% c("A","open_core","C"),]
s2$cohort <- c("A"="SHR-GV","open_core"="open-core","C"="GV-only")[s2$cohort]
cmp <- merge(old_tab %>% select(cohort,horizon_days,exposure,N,events,overall_p,nonlinear_p),
             s2 %>% select(cohort,horizon_days,exposure,N,events,overall_p,nonlinear_p),
             by=c("cohort","horizon_days","exposure"), suffixes=c("_legacy_rerun","_S2_current"))
cmp$N_match <- cmp$N_legacy_rerun==cmp$N_S2_current
cmp$events_match <- cmp$events_legacy_rerun==cmp$events_S2_current
cmp$overall_p_absdiff <- abs(cmp$overall_p_legacy_rerun - cmp$overall_p_S2_current)
cmp$nonlinear_p_absdiff <- abs(cmp$nonlinear_p_legacy_rerun - cmp$nonlinear_p_S2_current)
write_csv(cmp, file.path(OUT,"legacy_vs_current_S2_comparison.csv"))

cat("=== 主模型复现断言(线性版) ===\n"); print(chk)
cat("\n=== 旧 frame 复现当前 S2 对照 ===\n")
print(cmp %>% select(cohort,horizon_days,exposure,N_legacy_rerun,N_S2_current,N_match,overall_p_absdiff,nonlinear_p_absdiff))
cat("\n=== 新(primary frame)RCS ===\n")
print(new_tab %>% select(cohort,horizon_days,exposure,N,events,overall_p,nonlinear_p) %>%
        mutate(across(where(is.numeric), ~round(.,4))))
cat("RCS_FRAMEFIX_DONE\n")
sink()
