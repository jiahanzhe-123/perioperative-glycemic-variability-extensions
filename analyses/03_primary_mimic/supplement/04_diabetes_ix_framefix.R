# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 04_diabetes_ix_bloodonly_framefix.R
# 修复:糖尿病交互分析同样缺失 frame_extra(SHR-GV 30d 出现 5013、GV-only 出现 10325)。
# 修复:所有 MIMIC 主队列交互模型固定在主模型 frame 上,并对 frame N/events 做断言。
# 设计说明(保持原预设不变):非 GV-only 队列的交互模型为单暴露设计(逐个 gv/shr/mean × diabetes),
# 不同时纳入另一暴露;harmonized MIMIC/eICU 部分使用既有 harmonized 数据集,不受影响,原样保留。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr); library(sandwich); library(lmtest)})
ROOT <- PGV("mimic_work")
OUT  <- file.path(ROOT,"reviewer_issue_verification_20260728/01_frame_repair/diabetes_ix")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sink(file.path(ROOT,"reviewer_issue_verification_20260728/logs/04_diabetes_ix_framefix.log"), split=TRUE)

dat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA","NULL"))
num <- c("shr","gv","age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first",
         "creat_postop_first","wbc_postop_first","albumin_adm_first","hgb_postop_first",
         "platelets_postop_first","sofa_24h","survival_time_days","postop_30d_death_flag","one_year_death_flag",
         "glucose_mean_postop_24h","diabetes")
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
MISS_RULE <- "complete-case on outcome, exposure(s), model covariates and frame_extra (bmi+diabetes for GV-only both horizons and SHR-GV 30d); survival time > 0"

build_primary_frame <- function(dat, cohort_flag, cohort, horizon){
  d <- dat[dat[[cohort_flag]]==1,]
  ev <- if(horizon==30) "event_30d" else "event_1y"; tm <- if(horizon==30) "t30" else "t365"
  covs <- if(cohort=="GV-only") reduced else if(horizon==30) reduced else full
  fe <- if(cohort=="GV-only" || (cohort=="SHR-GV" && horizon==30)) c("bmi","diabetes") else NULL
  expo_vars <- if(cohort=="GV-only") "gv" else c("gv","shr")
  keep <- complete.cases(d[,c(tm,ev,expo_vars,strsplit(covs," \\+ ")[[1]],fe)]) & d[[tm]]>0
  list(d=d[keep,], ev=ev, tm=tm, covs=covs, fe=fe)
}
prim <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/02_primary_models_refit/final_primary_models.csv"), stringsAsFactors=FALSE)
assert_frame <- function(fr, cohort, horizon){
  ref <- prim[prim$cohort==cohort & prim$horizon_days==horizon,]; stopifnot(nrow(ref)==1)
  if(!(nrow(fr$d)==ref$N && sum(fr$d[[fr$ev]])==ref$events))
    stop(paste("FRAME MISMATCH:", cohort, horizon, nrow(fr$d), "vs", ref$N))
  invisible(TRUE)
}

run_ix <- function(fr, cohort, horizon, exposure_raw, zname){
  d <- fr$d
  d[[zname]] <- safe_z(d[[exposure_raw]])
  dm_n <- sum(d$diabetes==0); dm_y <- sum(d$diabetes==1)
  ev_y <- sum(d[[fr$ev]][d$diabetes==1]); ev_n <- sum(d[[fr$ev]][d$diabetes==0])
  f_int <- as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",zname,"*diabetes + ",fr$covs))
  f_add <- as.formula(paste0("Surv(",fr$tm,",",fr$ev,") ~ ",zname," + diabetes + ",fr$covs))
  fi <- coxph(f_int, data=d); fa <- coxph(f_add, data=d)
  wald_p <- summary(fi)$coefficients[paste0(zname,":diabetes"),"Pr(>|z|)"]
  lrt_p <- anova(fa,fi)$`Pr(>|Chi|)`[2]
  b <- coef(fi); V <- vcov(fi); ix <- paste0(zname,":diabetes")
  hr_ndm <- exp(b[zname]); se_ndm <- sqrt(V[zname,zname])
  lo_ndm <- exp(b[zname]-1.96*se_ndm); hi_ndm <- exp(b[zname]+1.96*se_ndm)
  b1 <- b[zname] + b[ix]; se1 <- sqrt(V[zname,zname] + V[ix,ix] + 2*V[zname,ix])
  hr_dm <- exp(b1); lo_dm <- exp(b1-1.96*se1); hi_dm <- exp(b1+1.96*se1)
  data.frame(db="MIMIC-IV", cohort=cohort, horizon_days=horizon, exposure=zname,
    N=nrow(d), N_nondiabetic=dm_n, N_diabetic=dm_y, events_nondiabetic=ev_n, events_diabetic=ev_y,
    HR_nondiabetic=unname(hr_ndm), lo_ndm=unname(lo_ndm), hi_ndm=unname(hi_ndm),
    HR_diabetic=unname(hr_dm), lo_dm=unname(lo_dm), hi_dm=unname(hi_dm),
    wald_p=wald_p, lrt_p=lrt_p,
    analytic_cohort=paste0(if(cohort=="SHR-GV") "final_v2_A" else if(cohort=="open-core") "final_v2_open_core" else "final_v2_C"," (",cohort,"), primary complete-case frame"),
    covariate_specification=if(fr$covs==reduced) paste0("reduced: ",reduced) else paste0("full: ",full),
    shr_jointly_adjusted=FALSE,
    design_note="single-exposure interaction model (prespecified); other exposure not co-adjusted",
    standardization_sample=paste0("z within this analytic sample (N=",nrow(d),")"),
    missing_data_rule=MISS_RULE, stringsAsFactors=FALSE)
}
jobs <- list(
  list("final_v2_A","SHR-GV",365,"gv","gv_z"), list("final_v2_A","SHR-GV",365,"shr","shr_z"),
  list("final_v2_A","SHR-GV",365,"glucose_mean_postop_24h","mean_z"),
  list("final_v2_A","SHR-GV",30,"gv","gv_z"), list("final_v2_A","SHR-GV",30,"shr","shr_z"),
  list("final_v2_C","GV-only",365,"gv","gv_z"), list("final_v2_C","GV-only",30,"gv","gv_z"),
  list("final_v2_open_core","open-core",365,"gv","gv_z"), list("final_v2_open_core","open-core",365,"shr","shr_z"))
mimic_ix <- dplyr::bind_rows(lapply(jobs, function(j){
  fr <- build_primary_frame(dat, j[[1]], j[[2]], j[[3]])
  assert_frame(fr, j[[2]], j[[3]])
  run_ix(fr, j[[2]], j[[3]], j[[4]], j[[5]])
}))

# ============ harmonized MIMIC + eICU(不受影响,原样保留) ============
mm <- read.csv(file.path(PGV("method_audit_work"),"00_audit","mimic_model_data_charlson_wo_diabetes.csv"))
ee <- read.csv(file.path(PGV("replication_work"),"03_eicu_harmonized","eicu_model_data.csv"))
run_hix <- function(df, db, sev, cluster=NULL){
  d <- df[df$in_main %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
  d <- d[d$in_landmark %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
  d <- d[!is.na(d$diabetes) & !is.na(d$glucose_sd_24h) & !is.na(d$glucose_mean_24h),]
  d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
  d$post_landmark_hosp_mortality <- as.integer(d$post_landmark_hosp_mortality %in% c(TRUE,"t","True","TRUE","true",1,"1"))
  cc <- c("post_landmark_hosp_mortality","glucose_sd_24h","glucose_mean_24h","age","sex","procedure_category","creatinine", sev)
  d <- d[complete.cases(d[,cc]),]
  d$gv_z <- as.numeric(scale(d$glucose_sd_24h)); d$mean_z <- as.numeric(scale(d$glucose_mean_24h))
  d$sex <- factor(d$sex); d$procedure_category <- factor(d$procedure_category)
  sev_f <- paste(sev, collapse=" + ")
  f_int <- as.formula(paste0("post_landmark_hosp_mortality ~ gv_z*diabetes + age + sex + procedure_category + creatinine + ", sev_f))
  f_add <- as.formula(paste0("post_landmark_hosp_mortality ~ gv_z + diabetes + age + sex + procedure_category + creatinine + ", sev_f))
  fi <- glm(f_int, data=d, family=binomial); fa <- glm(f_add, data=d, family=binomial)
  V_i <- if(is.null(cluster)) vcovHC(fi, type="HC1") else vcovCL(fi, cluster=d[[cluster]])
  wald_p <- coeftest(fi, vcov.=V_i)["gv_z:diabetes","Pr(>|z|)"]
  lrt_p <- anova(fa,fi,test="Chisq")$`Pr(>Chi)`[2]
  ci_i <- coeftest(fi, vcov.=V_i)
  b <- ci_i[,"Estimate"]; Sv <- V_i
  b0 <- b["gv_z"]; se0 <- sqrt(Sv["gv_z","gv_z"])
  b1 <- b["gv_z"] + b["gv_z:diabetes"]
  se1 <- sqrt(Sv["gv_z","gv_z"] + Sv["gv_z:diabetes","gv_z:diabetes"] + 2*Sv["gv_z","gv_z:diabetes"])
  data.frame(db=db, cohort="harmonized-main", horizon_days=NA, exposure="gv_z",
    N=nrow(d), N_nondiabetic=sum(d$diabetes==0), N_diabetic=sum(d$diabetes==1),
    events_nondiabetic=sum(d$post_landmark_hosp_mortality[d$diabetes==0]),
    events_diabetic=sum(d$post_landmark_hosp_mortality[d$diabetes==1]),
    HR_nondiabetic=exp(b0), lo_ndm=exp(b0-1.96*se0), hi_ndm=exp(b0+1.96*se0),
    HR_diabetic=exp(b1), lo_dm=exp(b1-1.96*se1), hi_dm=exp(b1+1.96*se1),
    wald_p=wald_p, lrt_p=lrt_p,
    analytic_cohort="harmonized analysis set (separate pipeline; unaffected by MIMIC primary frame)",
    covariate_specification=paste0("age + sex + procedure_category + creatinine + ", sev_f),
    shr_jointly_adjusted=FALSE, design_note="harmonized module; in_main & in_landmark flags",
    standardization_sample="z within harmonized sample", missing_data_rule="complete-case on harmonized variables",
    stringsAsFactors=FALSE)
}
hix <- dplyr::bind_rows(
  run_hix(mm, "MIMIC-IV (harmonized)", c("charlson_without_diabetes","sofa")),
  run_hix(ee, "eICU-CRD (harmonized)", "apachescore", cluster="hospitalid"))
ix_all <- dplyr::bind_rows(mimic_ix, hix)
write_csv(ix_all, file.path(OUT,"diabetes_ix_interaction_models_framefix.csv"))
print(ix_all %>% select(db,cohort,horizon_days,exposure,N,HR_nondiabetic,HR_diabetic,wald_p,lrt_p) %>% mutate(across(where(is.numeric), ~round(.,3))))
cat("PART4_FRAMEFIX_DONE\n")
sink()
