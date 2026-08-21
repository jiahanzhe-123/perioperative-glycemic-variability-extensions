# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 07_random_intercept_rerun.R — eICU hospital random-intercept logistic(M4 协调帧重跑)。
# 迁移自 stats_fix_phase1/scripts/02_eicu_random_intercept_rerun.R(逻辑不变,路径配置化)。
# Frame: 01_fit_harmonized_models.R ee filters + complete.cases on M4 variables
#   (gv10 + age + sex + procedure_category + diabetes + creatinine + apachescore + ns(mean_glu,3) + log_count + span_hours)
# Expected: N=7,115; events=130; hospitals=67 (must match 10_cross_database_results.csv M2/M3/M4).
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
# rm(list=ls()) 会清除页首 source 的配置,须重新加载:
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(lme4); library(splines)})
ROOT <- PGV("mimic_record_work")
dir.create(file.path(ROOT,"logs"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(ROOT,"logs","02_eicu_random_intercept.log"), split=TRUE)
ee <- read.csv(file.path(PGV("replication_work"),"03_eicu_harmonized","eicu_model_data.csv"), stringsAsFactors=FALSE)
ee <- ee[ee$in_main %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
ee <- ee[ee$in_landmark %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
ee$y <- as.integer(ee$hosp_mortality %in% c(TRUE,"t","True","TRUE","true",1,"1") &
                   !(ee$died_within_24h %in% c(TRUE,"t","True","TRUE","true",1,"1")))
ee$gv10 <- ee$glucose_sd_24h/10; ee$mean_glu <- ee$glucose_mean_24h
ee$diabetes <- as.integer(ee$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
ee$sex <- factor(ee$sex); ee$procedure_category <- factor(ee$procedure_category)
ee$log_count <- log(ee$glucose_n)
ee$span_hours <- ee$measurement_span_minutes/60
ee <- ee[!is.na(ee$gv10) & !is.na(ee$y),]
cat("ee after 08 filters (OLD wrong frame): N =", nrow(ee), "; events =", sum(ee$y), "\n")
m4vars <- c("y","gv10","age","sex","procedure_category","diabetes","creatinine","apachescore","mean_glu","log_count","span_hours")
ee2 <- ee[complete.cases(ee[, m4vars]),]
nh <- length(unique(ee2$hospitalid))
cat("M4 complete-case frame: N =", nrow(ee2), "; events =", sum(ee2$y), "; hospitals =", nh, "\n")
stopifnot(nrow(ee2)==7115, sum(ee2$y)==130, nh==67)
f <- glmer(y ~ gv10 + age + sex + procedure_category + diabetes + creatinine + apachescore + (1|hospitalid),
           data=ee2, family=binomial, control=glmerControl(optimizer="bobyqa"))
s <- summary(f)$coefficients["gv10",]
vc <- as.data.frame(VarCorr(f))
conv <- f@optinfo$conv$lme4$messages
tab <- data.frame(model_id="HARM_eICU_randinteract_rerun", database="eICU-CRD",
  outcome="post-landmark hospital mortality", model="hospital random-intercept logistic (sensitivity; M4 frame rerun)",
  N=nrow(ee2), events=sum(ee2$y), hospitals=nh,
  OR_per10=exp(s["Estimate"]), lo=exp(s["Estimate"]-1.96*s["Std. Error"]),
  hi=exp(s["Estimate"]+1.96*s["Std. Error"]), P=s["Pr(>|z|)"],
  hospital_variance=vc$vcov[1], hospital_sd=sqrt(vc$vcov[1]),
  convergence=if(is.null(conv)) "converged (no messages)" else paste(conv, collapse="; "),
  note="Rerun of the prespecified random-intercept sensitivity on the correct M4 frame (7,115/130/67 hospitals). The previously archived output (N=8,460, events=153) used the pre-complete-case frame and is superseded.",
  stringsAsFactors=FALSE)
write.csv(tab, file.path(ROOT,"results","EICU_RANDOM_INTERCEPT_RERUN.csv"), row.names=FALSE)
print(tab)
cat("EICU_RI_DONE\n"); sink()
