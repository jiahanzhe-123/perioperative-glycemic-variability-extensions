# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 03_run_primary_models_mice.R — 主分析:MICE(m=50, maxit=20, pmm, seed=20260726)+ Models A/B/C。
# 修复:移除确定性重复(gv10/log_count 改为各插补集内重算);
#       移除线性互补 frac_blood_gas 与近常数 t_lm_30 的预测角色;
#       Nelson–Aalen 仅保留 H365(含全部随访信息,H30 为其截断所致共线的根源)。
# 数据:analysis_base_bmi_repaired.csv(BMI-DQ-1 规则已由
#       analyses/01_cohort_construction/03_apply_bmi_plausibility.py 应用;
#       bmi<10 或 >80 -> NA -> MICE)。这是 analysis of record 的输入帧;
#       BMI 规则前的旧估计(0.977/0.992)不保留,仅存档。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
# rm(list=ls()) 会清除页首 source 的配置,须重新加载:
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
SEED <- 20260726L; set.seed(SEED)
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
dir.create(file.path(ROOT,"logs"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(ROOT,"results"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(ROOT,"logs","03b_mice_rerun.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$span_hours <- d$ps_span_hours; d$frac_central_lab <- d$ps_frac_central_lab
d$frac_blood_gas <- d$ps_frac_blood_gas; d$frac_poct <- d$ps_frac_poct; d$frac_icu_charted <- d$ps_frac_icu_charted
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
kSD <- K$gv_sd / 10
tgt$gender <- factor(tgt$gender)
tgt$procedure_cat6 <- factor(tgt$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
EXPAND <- c("albumin_adm_first","hgb_postop_first","wbc_postop_first","platelets_postop_first")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")

mi_data <- tgt[, c("stay_id","t_lm_30","t_lm_365","event_lm_30","event_lm_365",
                   "gv","mean_glu","glucose_count","span_hours",
                   "frac_central_lab","frac_blood_gas","frac_poct","frac_icu_charted",
                   COVS, EXPAND)]
mi_data$diabetes <- as.numeric(mi_data$diabetes)
mi_data$H365 <- nelsonaalen(mi_data, "t_lm_365", "event_lm_365")
pred <- make.predictorMatrix(mi_data); pred[,] <- 0L
meth <- make.method(mi_data); meth[] <- ""
cont_imp <- intersect(c("bmi","lactate_postop_first","creat_postop_first","sofa_24h", EXPAND), names(mi_data))
cont_imp <- cont_imp[colSums(is.na(mi_data[,cont_imp]))>0]
impute_preds <- setdiff(names(mi_data), c("stay_id","t_lm_30","frac_blood_gas"))
for (v in cont_imp) { meth[v] <- "pmm"; pred[v, impute_preds] <- 1L }
cat("MICE 插补变量:", paste(cont_imp, collapse=", "), "\n")
cat("预测集:", paste(impute_preds, collapse=", "), "\n")
t0 <- Sys.time()
imp <- mice(mi_data, m=50, maxit=20, method=meth, predictorMatrix=pred, seed=SEED, printFlag=FALSE)
cat("MICE 耗时:", round(difftime(Sys.time(), t0, units="mins"),1), "分钟\n")
saveRDS(imp, file.path(ROOT,"results","mice_m50_object.rds"))
lg <- imp$loggedEvents
if (!is.null(lg) && nrow(lg)) { write.csv(data.frame(lg), file.path(ROOT,"results","mice_logged_events.csv"), row.names=FALSE)
  cat("loggedEvents:", nrow(lg), "\n"); print(data.frame(lg)) } else cat("loggedEvents: 0\n")

rubin <- function(betas, vars){
  m <- length(betas); qb <- mean(unlist(betas)); U <- mean(unlist(vars))
  B <- var(unlist(betas)); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, fmi=(1+1/m)*B/Tm)
}
fml_A <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_fml))
fml_B <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml))
fml_C <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml,
           " + log_count + span_hours + frac_central_lab + frac_poct + frac_icu_charted"))
pool_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  betas <- list(A=list(b=c(),v=c()), B=list(b=c(),v=c()), C=list(b=c(),v=c()))
  for (i in 1:50) {
    dd <- complete(imp, i)
    dd$gv10 <- dd$gv/10; dd$log_count <- log(dd$glucose_count)
    for (mn in c("A","B","C")) {
      fit <- tryCatch(coxph(get(paste0("fml_",mn))(tm,ev), data=dd), error=function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit)
      betas[[mn]]$b <- c(betas[[mn]]$b, s$coefficients["gv10","coef"])
      betas[[mn]]$v <- c(betas[[mn]]$v, s$coefficients["gv10","se(coef)"]^2)
    }
  }
  for (mn in c("A","B","C")) {
    rb <- rubin(betas[[mn]]$b, betas[[mn]]$v)
    tv <- rb$beta/rb$se; pv <- 2*pt(abs(tv), df=rb$df, lower.tail=FALSE)
    pool_rows[[length(pool_rows)+1]] <- data.frame(
      model_id=sprintf("MICE_%s_%sd", mn, hz), analysis="primary (MICE m=50)",
      cohort="GV-only landmark", outcome=paste0("mortality by index day ",hz," among day-1 landmark survivors"),
      model=paste0("Model ",mn), N=nrow(tgt), events=sum(tgt[[ev]]),
      HR_per10=exp(rb$beta), lo_per10=exp(rb$beta-qt(.975,rb$df)*rb$se), hi_per10=exp(rb$beta+qt(.975,rb$df)*rb$se),
      P_per10=pv, HR_perSD=exp(rb$beta*kSD), lo_perSD=exp((rb$beta-qt(.975,rb$df)*rb$se)*kSD),
      hi_perSD=exp((rb$beta+qt(.975,rb$df)*rb$se)*kSD), P_perSD=pv, fmi=rb$fmi,
      n_imputations_fit=length(betas[[mn]]$b), stringsAsFactors=FALSE)
  }
}
pool_tab <- do.call(rbind, pool_rows)
write.csv(pool_tab, file.path(ROOT,"results","mice_pooled_models.csv"), row.names=FALSE)
print(pool_tab[,c("model_id","HR_per10","lo_per10","hi_per10","P_per10","HR_perSD","P_perSD","fmi")])

g <- function(id) pool_tab[pool_tab$model_id==id,]
prim <- data.frame(
  item=c("PRIMARY TEST: GV per 10 mg/dL (Model B, 30d)","PRIMARY: GV per fixed-cohort-SD (Model B, 30d)",
         "KEY SECONDARY: GV per 10 mg/dL (Model B, 365d)","KEY SECONDARY: GV per fixed-cohort-SD (Model B, 365d)"),
  HR=c(g("MICE_B_30d")$HR_per10, g("MICE_B_30d")$HR_perSD, g("MICE_B_365d")$HR_per10, g("MICE_B_365d")$HR_perSD),
  lo=c(g("MICE_B_30d")$lo_per10, g("MICE_B_30d")$lo_perSD, g("MICE_B_365d")$lo_per10, g("MICE_B_365d")$lo_perSD),
  hi=c(g("MICE_B_30d")$hi_per10, g("MICE_B_30d")$hi_perSD, g("MICE_B_365d")$hi_per10, g("MICE_B_365d")$hi_perSD),
  P=c(g("MICE_B_30d")$P_per10, g("MICE_B_30d")$P_perSD, g("MICE_B_365d")$P_per10, g("MICE_B_365d")$P_perSD),
  N=nrow(tgt), events_30d=sum(tgt$event_lm_30), events_365d=sum(tgt$event_lm_365),
  gv_sd_constant=K$gv_sd, gv_mean_constant=K$gv_mean, stringsAsFactors=FALSE)
write.csv(prim, file.path(ROOT,"results","04_primary_results.csv"), row.names=FALSE)
print(prim)
cat("PHASE3B_DONE\n")
sink()
