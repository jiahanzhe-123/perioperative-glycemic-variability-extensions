# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 09_sensitivity.R — open-core、Model A vs B 描述、替代 GV 指标(同样本配对)、
# 极端血糖与 winsorized、expanded covariate、极简模型。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/09_sensitivity.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$gender <- factor(d$gender)
d$procedure_cat6 <- factor(d$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
EXPAND <- c("albumin_adm_first","hgb_postop_first","wbc_postop_first","platelets_postop_first")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
rows <- list()
add <- function(...) { rows[[length(rows)+1]] <<- data.frame(..., stringsAsFactors=FALSE) }
fit_B <- function(dd, hz, xvar){
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  dd <- dd[complete.cases(dd[,c(tm,ev,xvar,"mean_glu",COVS)]) & dd[[tm]]>0,]
  dd$x10 <- dd[[xvar]]/10
  fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ x10 + ", rcs_mean, " + ", covs_fml)), data=dd)
  s <- summary(fit)
  list(HR=s$coefficients["x10","exp(coef)"], lo=s$conf.int["x10","lower .95"],
       hi=s$conf.int["x10","upper .95"], P=s$coefficients["x10","Pr(>|z|)"],
       N=nrow(dd), events=sum(dd[[ev]]))
}

# ---- 1. open-core 队列(次要 #3) ----
oc <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2 & d$final_v2_open_core==1,]
for (hz in c("30","365")) {
  r <- fit_B(oc, hz, "gv")
  add(model_id=paste0("OPENCORE_B_",hz,"d"), analysis="key secondary #3 (open-core)",
      cohort="open-core landmark", outcome=paste0("mortality by index day ",hz," among day-1 landmark survivors"),
      model="Model B", N=r$N, events=r$events, HR_per10=r$HR, lo_per10=r$lo, hi_per10=r$hi, P=r$P)
}

# ---- 2. 替代 GV 指标(同样本配对 SD 对照) ----
cc <- tgt[complete.cases(tgt[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","gv","mean_glu",COVS)]) & tgt$t_lm_30>0,]
alt_metrics <- list(
  list(var="ps_cv", label="CV", min_n=2, scale=100),
  list(var="ps_arv", label="ARV", min_n=3, scale=10),
  list(var="ps_tw_arv", label="time-weighted ARV", min_n=3, scale=10),
  list(var="ps_mad_glucose", label="MAD", min_n=2, scale=10),
  list(var="ps_iqr_glucose", label="IQR", min_n=2, scale=10))
for (am in alt_metrics) {
  dd <- cc[!is.na(cc[[am$var]]) & cc$glucose_count>=am$min_n,]
  for (hz in c("30","365")) {
    r1 <- fit_B(dd, hz, am$var)
    add(model_id=paste0("ALT_",am$label,"_",hz,"d"), analysis="sensitivity (alternative GV definition)",
        cohort=paste0("GV-only landmark, ", am$label, "-computable subset"), metric=am$label,
        N=r1$N, events=r1$events, HR=r1$HR, lo=r1$lo, hi=r1$hi, P=r1$P,
        note=paste0("per ", am$scale, " units; same subset as paired SD row below"))
    dd$gv <- dd$gv  # SD comparator on identical patients
    r2 <- fit_B(dd, hz, "gv")
    add(model_id=paste0("ALTpairedSD_",am$label,"_",hz,"d"), analysis="sensitivity (paired SD comparator)",
        cohort=paste0("GV-only landmark, ", am$label, "-computable subset"), metric="SD (paired)",
        N=r2$N, events=r2$events, HR=r2$HR, lo=r2$lo, hi=r2$hi, P=r2$P,
        note="per 10 mg/dL; identical patients as the alternative-metric row")
  }
}

# ---- 3. 极端血糖与 winsorized ----
for (adj in c("prop_lt70","prop_gt180","minmax","burden")) {
  dd <- cc
  extra <- switch(adj,
    prop_lt70=" + ps_prop_lt70", prop_gt180=" + ps_prop_gt180",
    minmax=" + ps_min_glucose + ps_max_glucose", burden=" + ps_prop_lt70 + ps_prop_gt180")
  for (hz in c("30","365")) {
    tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
    dd$gv10 <- dd$gv/10
    wtxt <- ""
    fit <- tryCatch(
      withCallingHandlers(
        coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml, extra)), data=dd),
        warning=function(w){ wtxt <<- conditionMessage(w); invokeRestart("muffleWarning") }),
      error=function(e) NULL)
    if (is.null(fit)) { add(model_id=paste0("EXT_",adj,"_",hz,"d"), note="model failed"); next }
    s <- summary(fit)
    add(model_id=paste0("EXT_",adj,"_",hz,"d"), analysis="sensitivity (extreme-glucose adjustment; definition dependence only)",
        cohort="GV-only landmark", outcome=paste0("mortality by index day ",hz," among day-1 landmark survivors"),
        model=paste0("Model B + ", adj), N=nrow(dd), events=sum(dd[[ev]]),
        HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
        hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"],
        note=paste0("数学耦合,不作'独立于极端值的生物效应'解释", if(nzchar(wtxt)) paste0("; WARNING: ", wtxt)))
  }
}
# winsorized 序列(队列 0.5/99.5 分位截断后重算 GV)
ser <- read.csv(file.path(ROOT,"data","series_priority.csv"), stringsAsFactors=FALSE)
lo <- quantile(ser$value, .005); hi <- quantile(ser$value, .995)
ser$w <- pmin(pmax(ser$value, lo), hi)
wg <- aggregate(w ~ stay_id, data=ser, FUN=function(x) if(length(x)>=2) sd(x) else NA_real_)
names(wg)[2] <- "gv_wins"
dd <- merge(cc[, !(names(cc) %in% c("gv"))], wg, by="stay_id")
for (hz in c("30","365")) {
  r <- fit_B(dd, hz, "gv_wins")
  add(model_id=paste0("WINS_",hz,"d"), analysis="sensitivity (winsorized sequence, 0.5/99.5 percentile)",
      cohort="GV-only landmark", N=r$N, events=r$events, HR_per10=r$HR, lo_per10=r$lo, hi_per10=r$hi, P=r$P)
}

# ---- 4. expanded covariate model & 极简模型 ----
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  dd <- cc[complete.cases(cc[,c(EXPAND)]),]
  dd$gv10 <- dd$gv/10
  fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml, " + ",
           paste(EXPAND, collapse=" + "))), data=dd)
  s <- summary(fit)
  add(model_id=paste0("EXPAND_",hz,"d"), analysis="sensitivity (expanded covariate model)",
      cohort="GV-only landmark (expanded-complete subset)", N=nrow(dd), events=sum(dd[[ev]]),
      HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
      hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"])
  dd2 <- cc[complete.cases(cc[,c("age_at_admission","gender","charlson_without_diabetes","sofa_24h")]),]
  dd2$gv10 <- dd2$gv/10
  fit2 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + age_at_admission + gender + charlson_without_diabetes + sofa_24h")), data=dd2)
  s2 <- summary(fit2)
  add(model_id=paste0("MINIMAL_",hz,"d"), analysis="sensitivity (minimal no-BMI model)",
      cohort="GV-only landmark", N=nrow(dd2), events=sum(dd2[[ev]]),
      HR_per10=s2$coefficients["gv10","exp(coef)"], lo=s2$conf.int["gv10","lower .95"],
      hi=s2$conf.int["gv10","upper .95"], P=s2$coefficients["gv10","Pr(>|z|)"])
}

allcols <- unique(unlist(lapply(rows, names)))
rows2 <- lapply(rows, function(x){ for (cn in setdiff(allcols, names(x))) x[[cn]] <- NA; x[, allcols, drop=FALSE] })
tab <- do.call(rbind, rows2)
write.csv(tab, file.path(ROOT,"results","06_sensitivity_results.csv"), row.names=FALSE)
print(tab[, intersect(c("model_id","analysis","cohort","N","events","HR_per10","lo_per10","hi_per10","P","HR","lo","hi","note"), names(tab))], max=80)
cat("PHASE9_DONE\n")
sink()
