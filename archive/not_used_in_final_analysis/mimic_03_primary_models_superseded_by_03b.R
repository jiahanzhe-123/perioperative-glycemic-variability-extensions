#!/usr/bin/env Rscript
# 03_primary_models.R — 重构版主要分析
# 唯一主要研究问题:day-1 landmark 幸存者中,date-anchored 24h SD-based GV 与
# 索引日后30日内死亡(风险随访自 landmark 起)的条件关联。
# 主要模型 = Model B(固定临床协变量 + RCS 灵活 mean glucose + 线性 GV)。
# 主要检验唯一:Model B 30d 中 GV 的 Wald 检验。
# 缺失数据主方案:MICE m=50 maxit=20(Nelson–Aalen H + event + 全部协变量);
# complete-case 为敏感性。标准化常数在目标队列计算一次并固定。
# 效应尺度:per-10 mg/dL(模型内 gv10)与 per-固定队列SD(由 gv10 系数精确换算,P 不变)。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.libPaths(c("~/cardiac_glucose_rebuild_20260728/rlib", .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite)})
ROOT <- normalizePath("~/cardiac_glucose_rebuild_20260728")
dir.create(file.path(ROOT,"results"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(ROOT,"logs/03_primary_models.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base.csv"), stringsAsFactors=FALSE)
# pandas 布尔列读出为 "True"/"False" 字符串,统一转回逻辑值
for (bc in c("landmark_eligible","cabg_flag","open_valve_flag","aortic_surgery_flag",
             "transplant_vad_flag","congenital_cardiac_flag","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
base$event_lm_30 <- as.integer(base$event_lm_30); base$event_lm_365 <- as.integer(base$event_lm_365)
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(base$stay_id), !anyDuplicated(feat$stay_id))
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id", all.x=FALSE, all.y=FALSE)
stopifnot(!anyDuplicated(d$subject_id))   # 每位患者仅一个主要 stay
# 分析暴露一律取 priority 序列版本
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose
d$glucose_count <- d$ps_glucose_count; d$span_hours <- d$ps_span_hours
d$frac_central_lab <- d$ps_frac_central_lab; d$frac_blood_gas <- d$ps_frac_blood_gas
d$frac_poct <- d$ps_frac_poct; d$frac_icu_charted <- d$ps_frac_icu_charted

# ---- 目标队列:landmark 合格 + >=2 次合格测量(priority 序列) ----
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
stopifnot(all(tgt$t_lm_30 >= 0), all(tgt$t_lm_365 >= 0))
cat("目标队列 N =", nrow(tgt), "; 30d 事件 =", sum(tgt$event_lm_30),
    "; 365d 事件 =", sum(tgt$event_lm_365), "\n")

# ---- 固定标准化常数 ----
K <- list(
  gv_mean = mean(tgt$gv), gv_sd = sd(tgt$gv),
  mean_glu_mean = mean(tgt$mean_glu), mean_glu_sd = sd(tgt$mean_glu),
  gv_q25 = unname(quantile(tgt$gv, .25)), gv_q75 = unname(quantile(tgt$gv, .75)),
  mean_glu_knots = unname(quantile(tgt$mean_glu, c(.05,.35,.65,.95))),
  gv_knots_30d = unname(quantile(tgt$gv, c(.10,.50,.90))),
  gv_knots_365d = unname(quantile(tgt$gv, c(.05,.35,.65,.95))),
  cohort = "GV-only, day-1 landmark survivors, >=2 eligible measurements (priority series)",
  seed = SEED)
write_json(K, file.path(ROOT,"results","standardization_constants.json"), pretty=TRUE, auto_unbox=TRUE)
tgt$gv10 <- tgt$gv / 10
kSD <- K$gv_sd / 10   # per-10 -> per-SD 精确换算因子(HR^kSD,SE*kSD,P 不变)

# ---- 固定协变量 ----
tgt$gender <- factor(tgt$gender)
tgt$procedure_cat6 <- factor(tgt$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
EXPAND <- c("albumin_adm_first","hgb_postop_first","wbc_postop_first","platelets_postop_first")
if ("anchor_year_group" %in% names(tgt) && any(!is.na(tgt$anchor_year_group))) {
  tgt$anchor_year_group <- factor(tgt$anchor_year_group); COVS <- c(COVS, "anchor_year_group")
  cat("calendar era 已纳入\n")
} else cat("calendar era 缺失(anchor_year_group 不可用),见限制文档\n")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")

fml_A <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_fml))
fml_B <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml))
fml_C <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml,
           " + log_count + span_hours + frac_central_lab + frac_poct + frac_icu_charted"))

report_row <- function(fit, beta_se, model_id, analysis, cohort, outcome, model, N, events, ph_p=NA){
  cf <- beta_se[1]; se <- beta_se[2]; z <- cf/se
  data.frame(model_id=model_id, analysis=analysis, cohort=cohort, outcome=outcome, model=model,
             N=N, events=events,
             HR_per10=exp(cf), lo_per10=exp(cf-1.96*se), hi_per10=exp(cf+1.96*se), P_per10=2*pnorm(abs(z), lower.tail=FALSE),
             HR_perSD=exp(cf*kSD), lo_perSD=exp((cf-1.96*se)*kSD), hi_perSD=exp((cf+1.96*se)*kSD),
             P_perSD=2*pnorm(abs(z), lower.tail=FALSE), ph_global_p=ph_p, stringsAsFactors=FALSE)
}
beta_gv10 <- function(fit) c(summary(fit)$coefficients["gv10","coef"], summary(fit)$coefficients["gv10","se(coef)"])

# ---- complete-case(敏感性 + 对照) ----
tgt$log_count <- log(tgt$glucose_count)
cc <- tgt[complete.cases(tgt[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","gv","gv10","mean_glu",COVS)]) & tgt$t_lm_30>0,]
cat("complete-case N =", nrow(cc), "; 30d 事件 =", sum(cc$event_lm_30), "; 365d 事件 =", sum(cc$event_lm_365), "\n")
cc_res <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  for (mn in c("A","B","C")) {
    fit <- coxph(get(paste0("fml_",mn))(tm,ev), data=cc, x=TRUE, y=TRUE)
    ph <- tryCatch(cox.zph(fit)$table["GLOBAL","p"], error=function(e) NA)
    cc_res[[length(cc_res)+1]] <- report_row(fit, beta_gv10(fit), sprintf("CC_%s_%sd", mn, hz),
      "complete-case sensitivity", "GV-only landmark",
      paste0("mortality by index day ",hz," among day-1 landmark survivors"),
      paste0("Model ",mn), nrow(cc), sum(cc[[ev]]), ph)
  }
}
cc_tab <- do.call(rbind, cc_res)
write.csv(cc_tab, file.path(ROOT,"results","cc_models_30d_365d.csv"), row.names=FALSE)
print(cc_tab[,c("model_id","N","events","HR_per10","lo_per10","hi_per10","P_per10","HR_perSD","P_perSD","ph_global_p")])

# ---- MICE(主要方案) ----
cat("各变量缺失率(%):\n"); print(round(colMeans(is.na(tgt[,c(COVS,EXPAND)]))*100, 2))
mi_data <- tgt[, c("stay_id","t_lm_30","t_lm_365","event_lm_30","event_lm_365",
                   "gv","gv10","mean_glu","glucose_count","span_hours","log_count",
                   "frac_central_lab","frac_blood_gas","frac_poct","frac_icu_charted",
                   COVS, EXPAND)]
mi_data$diabetes <- as.numeric(mi_data$diabetes)
mi_data$H30  <- nelsonaalen(mi_data, "t_lm_30", "event_lm_30")
mi_data$H365 <- nelsonaalen(mi_data, "t_lm_365", "event_lm_365")
pred <- make.predictorMatrix(mi_data); pred[,] <- 0L
meth <- make.method(mi_data); meth[] <- ""
cont_imp <- intersect(c("bmi","lactate_postop_first","creat_postop_first","sofa_24h", EXPAND), names(mi_data))
cont_imp <- cont_imp[colSums(is.na(mi_data[,cont_imp]))>0]
impute_preds <- setdiff(names(mi_data), "stay_id")
for (v in cont_imp) { meth[v] <- "pmm"; pred[v, impute_preds] <- 1L }
cat("MICE 插补变量:", paste(cont_imp, collapse=", "), "\n")
t0 <- Sys.time()
imp <- mice(mi_data, m=50, maxit=20, method=meth, predictorMatrix=pred, seed=SEED, printFlag=FALSE)
cat("MICE 耗时:", round(difftime(Sys.time(), t0, units="mins"),1), "分钟\n")
saveRDS(imp, file.path(ROOT,"results","mice_m50_object.rds"))
lg <- imp$loggedEvents
if (!is.null(lg) && nrow(lg)) { write.csv(data.frame(lg), file.path(ROOT,"results","mice_logged_events.csv"), row.names=FALSE)
  cat("loggedEvents:", nrow(lg), "\n") } else cat("loggedEvents: 0\n")

rubin <- function(betas, vars){
  m <- length(betas); qb <- mean(unlist(betas)); U <- mean(unlist(vars))
  B <- var(unlist(betas)); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, fmi=(1+1/m)*B/Tm)
}
pool_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  betas <- list(A=list(b=c(),v=c()), B=list(b=c(),v=c()), C=list(b=c(),v=c()))
  conv_warn <- 0L
  for (i in 1:50) {
    dd <- complete(imp, i)
    for (mn in c("A","B","C")) {
      fit <- tryCatch(coxph(get(paste0("fml_",mn))(tm,ev), data=dd),
                      warning=function(w){conv_warn <<- conv_warn+1L; invokeRestart("muffleWarning")},
                      error=function(e) NULL)
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
      n_imputations_fit=length(betas[[mn]]$b), convergence_warnings=conv_warn, stringsAsFactors=FALSE)
  }
}
pool_tab <- do.call(rbind, pool_rows)
write.csv(pool_tab, file.path(ROOT,"results","mice_pooled_models.csv"), row.names=FALSE)
print(pool_tab[,c("model_id","N","events","HR_per10","lo_per10","hi_per10","P_per10","HR_perSD","P_perSD","fmi","n_imputations_fit")])

# ---- 主要/关键次要结果(04_primary_results.csv) ----
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
cat("PHASE3_DONE\n")
sink()
