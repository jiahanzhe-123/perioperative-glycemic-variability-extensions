# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# Harmonized cross-database replication: Models 0-3 + prespecified sensitivities
# MIMIC: sandwich HC SE; eICU: hospital-clustered robust SE (+glmer sensitivity)
set.seed(20260724)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(sandwich); library(lmtest); library(lme4); library(dplyr); library(readr); library(broom)})

ROOT <- PGV("replication_work")
mimic <- read_csv(file.path(ROOT,"02_mimic_harmonized/mimic_model_data.csv"), show_col_types=FALSE)
eicu  <- read_csv(file.path(ROOT,"03_eicu_harmonized/eicu_model_data.csv"), show_col_types=FALSE)
sink(file.path(ROOT,"logs/models_run.log"), split=TRUE)

# ---- BMI 冻结规则(protocol §6): 两库 BMI 完整率均 ≥85% 且加入 BMI 后两库保留率均 ≥85% ----
mm0 <- mimic %>% filter(in_main, in_landmark, glucose_n>=2)
ee0 <- eicu %>% filter(in_main, in_landmark, glucose_n>=2)
bmi_comp <- c(MIMIC=mean(!is.na(mm0$bmi)), eICU=mean(!is.na(ee0$bmi)))
cc_m <- sum(complete.cases(mm0[,c("post_landmark_hosp_mortality","glucose_sd_24h","glucose_mean_24h","age","sex","procedure_category","diabetes","creatinine","bmi","charlson","sofa")]))/nrow(mm0)
cc_e <- sum(complete.cases(ee0[,c("post_landmark_hosp_mortality","glucose_sd_24h","glucose_mean_24h","age","sex","procedure_category","diabetes","creatinine","bmi","apachescore")]))/nrow(ee0)
BMI_INCLUDE <- all(bmi_comp >= 0.85) && all(c(cc_m, cc_e) >= 0.85)
cat("BMI rule: completeness MIMIC", round(bmi_comp[1],3), "eICU", round(bmi_comp[2],3),
    "; retention with BMI MIMIC", round(cc_m,3), "eICU", round(cc_e,3),
    "-> BMI_INCLUDE =", BMI_INCLUDE, "\n")

fit_one <- function(df, db, cohort_label, outcome, severity, cluster_var=NULL, exposure="gv_z", extra="") {
  # 构造公式
  f1 <- paste0(outcome, " ~ ", exposure, " + age + sex + procedure_category + diabetes + creatinine", if (BMI_INCLUDE) " + bmi" else "")
  f0 <- paste0(outcome, " ~ ", exposure)
  f2 <- paste0(f1, " + ", severity)
  f3 <- paste0(f2, " + mean_z")
  if (extra != "") { f2 <- paste0(f2, " + ", extra); f3 <- paste0(f3, " + ", extra) }
  # z 标准化(本库、本分析集内): 用 M2 的完整病例集
  sev_vars <- strsplit(severity, " \\+ ")[[1]]
  cc_vars <- unique(c(outcome, "glucose_sd_24h","glucose_cv_24h","glucose_mean_24h","age","sex","procedure_category","diabetes","creatinine", if (BMI_INCLUDE) "bmi", sev_vars,
                      if(extra!="") strsplit(extra," \\+ ")[[1]] else NULL))
  d <- df[complete.cases(df[, cc_vars]), ]
  # 暴露已在外部 z 化? 未——在此按 cc 集 z 化
  d$gv_z <- as.numeric(scale(d$glucose_sd_24h)); d$mean_z <- as.numeric(scale(d$glucose_mean_24h)); d$cv_z <- as.numeric(scale(d$glucose_cv_24h))
  if (exposure=="cv_z") { f0 <- sub("gv_z","cv_z",f0); f1 <- sub("gv_z","cv_z",f1); f2 <- sub("gv_z","cv_z",f2); f3 <- sub("gv_z","cv_z",f3) }
  out <- list()
  for (mn in c("M0","M1","M2","M3")) {
    fm <- switch(mn, M0=f0, M1=f1, M2=f2, M3=f3)
    fit <- glm(as.formula(fm), data=d, family=binomial)
    if (!is.null(cluster_var)) {
      V <- vcovCL(fit, cluster=d[[cluster_var]])
    } else {
      V <- vcovHC(fit, type="HC1")
    }
    ct <- coeftest(fit, vcov.=V)
    rn <- rownames(ct)
    i <- which(rn==exposure)
    out[[mn]] <- data.frame(db=db, cohort=cohort_label, outcome=outcome, model=mn,
      n=nrow(d), events=sum(d[[outcome]]==1),
      exposure=exposure, extra=extra, beta=ct[i,1], se=ct[i,2], or=exp(ct[i,1]),
      ci_lo=exp(ct[i,1]-1.96*ct[i,2]), ci_hi=exp(ct[i,1]+1.96*ct[i,2]), p=ct[i,4],
      formula=fm, stringsAsFactors=FALSE)
    assign(paste0("fit_",mn), fit)
  }
  list(res=do.call(rbind,out), fits=list(M0=fit_M0,M1=fit_M1,M2=fit_M2,M3=fit_M3), data=d, vcov_note=ifelse(is.null(cluster_var),"HC1","cluster"))
}

res_all <- list()

# ---------- 主分析 ----------
# MIMIC main
mm <- mimic %>% filter(in_main, in_landmark, glucose_n>=2) %>%
  mutate(sex=factor(sex), procedure_category=ifelse(is.na(procedure_category),"unknown",procedure_category) %>% factor() %>% relevel(ref="cabg"))
r <- fit_one(mm, "MIMIC-IV","main","post_landmark_hosp_mortality","charlson + sofa")
res_all[["mimic_main"]] <- r
saveRDS(r$fits, file.path(ROOT,"05_models/mimic_main_models.rds"))

# eICU main
ee <- eicu %>% filter(in_main, in_landmark, glucose_n>=2) %>%
  mutate(sex=factor(sex), procedure_category=ifelse(is.na(procedure_category),"unknown",procedure_category) %>% factor() %>% relevel(ref="cabg"))
r <- fit_one(ee, "eICU-CRD","main","post_landmark_hosp_mortality","apachescore", cluster_var="hospitalid")
res_all[["eicu_main"]] <- r
saveRDS(r$fits, file.path(ROOT,"05_models/eicu_main_models.rds"))

primary <- do.call(rbind, lapply(res_all, `[[`, "res"))
primary$analysis <- "primary"

# ---------- 敏感性 ----------
sens <- list()
run_sens <- function(df, db, label, outcome, severity, cluster_var=NULL, exposure="gv_z", extra="", fname) {
  d <- df %>% mutate(sex=factor(sex), procedure_category=ifelse(is.na(procedure_category),"unknown",procedure_category) %>% factor() %>% relevel(ref="cabg"))
  r <- tryCatch(fit_one(d, db, label, outcome, severity, cluster_var, exposure, extra),
                error=function(e){cat("SENS FAIL", db, label, outcome, extra, ":", conditionMessage(e), "\n"); NULL})
  if (!is.null(r)) saveRDS(r$fits, file.path(ROOT,paste0("05_models/",fname,".rds")))
  r
}

# S1 eICU broad
sens[["s1_eicu_broad"]] <- run_sens(eicu %>% filter(in_sensc, in_landmark, glucose_n>=2), "eICU-CRD","sensC_broad","post_landmark_hosp_mortality","apachescore","hospitalid",fname="s1_eicu_broad")
# S2 +aortic (both)
sens[["s2a_mimic"]] <- run_sens(mimic %>% filter(in_sensA, in_landmark, glucose_n>=2), "MIMIC-IV","sensA_aortic","post_landmark_hosp_mortality","charlson + sofa",fname="s2a_mimic")
sens[["s2b_eicu"]] <- run_sens(eicu %>% filter(in_sensa, in_landmark, glucose_n>=2), "eICU-CRD","sensA_aortic","post_landmark_hosp_mortality","apachescore","hospitalid",fname="s2b_eicu")
# S3 all high-spec (both)
sens[["s3a_mimic"]] <- run_sens(mimic %>% filter(in_sensB, in_landmark, glucose_n>=2), "MIMIC-IV","sensB_all_highspec","post_landmark_hosp_mortality","charlson + sofa",fname="s3a_mimic")
sens[["s3b_eicu"]] <- run_sens(eicu %>% filter(in_sensb, in_landmark, glucose_n>=2), "eICU-CRD","sensB_all_highspec","post_landmark_hosp_mortality","apachescore","hospitalid",fname="s3b_eicu")
# S4 non-landmark hospital mortality
sens[["s4a_mimic"]] <- run_sens(mimic %>% filter(in_main, glucose_n>=2), "MIMIC-IV","main","hosp_mortality","charlson + sofa",fname="s4a_mimic")
sens[["s4b_eicu"]] <- run_sens(eicu %>% filter(in_main, glucose_n>=2), "eICU-CRD","main","hosp_mortality","apachescore","hospitalid",fname="s4b_eicu")
# S5 ICU mortality
sens[["s5a_mimic"]] <- run_sens(mimic %>% filter(in_main, glucose_n>=2), "MIMIC-IV","main","icu_mortality","charlson + sofa",fname="s5a_mimic")
sens[["s5b_eicu"]] <- run_sens(eicu %>% filter(in_main, glucose_n>=2), "eICU-CRD","main","icu_mortality","apachescore","hospitalid",fname="s5b_eicu")
# S6 CV 替代 SD
sens[["s6a_mimic"]] <- run_sens(mimic %>% filter(in_main, in_landmark, glucose_n>=2), "MIMIC-IV","main","post_landmark_hosp_mortality","charlson + sofa",exposure="cv_z",fname="s6a_mimic")
sens[["s6b_eicu"]] <- run_sens(eicu %>% filter(in_main, in_landmark, glucose_n>=2), "eICU-CRD","main","post_landmark_hosp_mortality","apachescore","hospitalid",exposure="cv_z",fname="s6b_eicu")
# S7 + glucose_n + measurement_span
sens[["s7a_mimic"]] <- run_sens(mimic %>% filter(in_main, in_landmark, glucose_n>=2), "MIMIC-IV","main","post_landmark_hosp_mortality","charlson + sofa",extra="glucose_n + measurement_span_minutes",fname="s7a_mimic")
sens[["s7b_eicu"]] <- run_sens(eicu %>% filter(in_main, in_landmark, glucose_n>=2), "eICU-CRD","main","post_landmark_hosp_mortality","apachescore","hospitalid",extra="glucose_n + measurement_span_minutes",fname="s7b_eicu")
# S9 eICU 医院>=10 例
hosp_ok <- eicu %>% filter(in_main, in_landmark, glucose_n>=2) %>% count(hospitalid) %>% filter(n>=10) %>% pull(hospitalid)
sens[["s9_eicu"]] <- run_sens(eicu %>% filter(in_main, in_landmark, glucose_n>=2, hospitalid %in% hosp_ok), "eICU-CRD","main_hospge10","post_landmark_hosp_mortality","apachescore","hospitalid",fname="s9_eicu")
# S10 排除 glucose_n==2
sens[["s10a_mimic"]] <- run_sens(mimic %>% filter(in_main, in_landmark, glucose_n>=3), "MIMIC-IV","main_ge3","post_landmark_hosp_mortality","charlson + sofa",fname="s10a_mimic")
sens[["s10b_eicu"]] <- run_sens(eicu %>% filter(in_main, in_landmark, glucose_n>=3), "eICU-CRD","main_ge3","post_landmark_hosp_mortality","apachescore","hospitalid",fname="s10b_eicu")

sens_res <- do.call(rbind, lapply(sens[!sapply(sens,is.null)], `[[`, "res"))

# S8 eICU random intercept (主队列, M2/M3)
d8 <- res_all[["eicu_main"]]$data
ri_out <- data.frame()
ri_status <- "ok"
for (mn in c("M2","M3")) {
  fm <- res_all[["eicu_main"]]$res[res_all[["eicu_main"]]$res$model==mn,"formula"]
  ri <- tryCatch(glmer(as.formula(paste0(fm," + (1|hospitalid)")), data=d8, family=binomial,
                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))),
                 error=function(e){ri_status<<-paste("glmer",mn,"failed:",conditionMessage(e)); NULL})
  if (is.null(ri)) next
  s <- summary(ri)$coefficients["gv_z",]
  ri_out <- rbind(ri_out, data.frame(db="eICU-CRD", cohort="main", outcome="post_landmark_hosp_mortality",
    model=paste0(mn,"_random_intercept"), n=nobs(ri), events=sum(d8$post_landmark_hosp_mortality==1),
    exposure="gv_z", beta=s[1], se=s[2], or=exp(s[1]), ci_lo=exp(s[1]-1.96*s[2]), ci_hi=exp(s[1]+1.96*s[2]),
    p=s[4], formula=paste0(fm," + (1|hospitalid)")))
  saveRDS(ri, file.path(ROOT,paste0("05_models/s8_eicu_ri_",mn,".rds")))
}
cat("S8 random intercept status:", ri_status, "\n")

# ---------- 汇总输出 ----------
all_res <- bind_rows(primary, sens_res, ri_out)
# attenuation
att <- all_res %>% filter(model %in% c("M2","M3")) %>%
  select(db,cohort,outcome,exposure,extra,model,beta) %>%
  tidyr::pivot_wider(names_from=model, values_from=beta) %>%
  mutate(attenuation_pct = ifelse(abs(M2)>1e-6 & M2>0, 100*(M2-M3)/M2, NA_real_))
write_csv(all_res, file.path(ROOT,"05_models/all_model_results.csv"))
write_csv(att, file.path(ROOT,"05_models/attenuation.csv"))

# 标准化参数
std_params <- rbind(
  data.frame(db="MIMIC-IV", cohort="main", var="glucose_sd_24h", mean=mean(res_all[["mimic_main"]]$data$glucose_sd_24h,na.rm=T), sd=sd(res_all[["mimic_main"]]$data$glucose_sd_24h,na.rm=T)),
  data.frame(db="MIMIC-IV", cohort="main", var="glucose_mean_24h", mean=mean(res_all[["mimic_main"]]$data$glucose_mean_24h,na.rm=T), sd=sd(res_all[["mimic_main"]]$data$glucose_mean_24h,na.rm=T)),
  data.frame(db="eICU-CRD", cohort="main", var="glucose_sd_24h", mean=mean(res_all[["eicu_main"]]$data$glucose_sd_24h,na.rm=T), sd=sd(res_all[["eicu_main"]]$data$glucose_sd_24h,na.rm=T)),
  data.frame(db="eICU-CRD", cohort="main", var="glucose_mean_24h", mean=mean(res_all[["eicu_main"]]$data$glucose_mean_24h,na.rm=T), sd=sd(res_all[["eicu_main"]]$data$glucose_mean_24h,na.rm=T)))
write_csv(std_params, file.path(ROOT,"05_models/standardization_params.csv"))

# 完整系数(主模型 M2/M3)
full_coef <- function(fits, db){
  out <- data.frame()
  for (mn in names(fits)) {
    t <- tidy(fits[[mn]], conf.int=TRUE, exponentiate=FALSE)
    t$model <- mn; t$db <- db; out <- rbind(out,t)
  }
  out
}
fc <- rbind(full_coef(res_all[["mimic_main"]]$fits,"MIMIC-IV"), full_coef(res_all[["eicu_main"]]$fits,"eICU-CRD"))
write_csv(fc, file.path(ROOT,"05_models/primary_full_coefficients.csv"))

print(primary)
print(att)
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()
cat("MODELS_DONE\n")
