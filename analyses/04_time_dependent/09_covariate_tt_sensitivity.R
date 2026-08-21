# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 09_covariate_tt_sensitivity.R — 365d 协变量 PH 敏感性(age/lactate/SOFA log(time) 项),MICE 池化。
# 迁移自 stats_fix_phase1/scripts/03_mimic_365d_covariate_tt_sensitivity.R(逻辑不变,路径配置化)。
#   D1) 365d Model B + tt(age)+tt(lactate)+tt(SOFA): GV average estimate under covariate PH handling.
#   D2) 365d interval model + same tt terms: checks interval GV directions under covariate PH handling.
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
# rm(list=ls()) 会清除页首 source 的配置,须重新加载:
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
P35  <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs","03_mimic_365d_covariate_tt_sensitivity.log"), split=TRUE)

K <- fromJSON(file.path(P35,"results","standardization_constants.json"))
imp <- readRDS(file.path(P35,"results","mice_m50_object.rds"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
TT <- c("age_at_admission","lactate_postop_first","sofa_24h")
tt_fml <- paste(paste0("tt(", TT, ")"), collapse=" + ")

rubin <- function(betas, vars){
  m <- length(betas); qb <- mean(unlist(betas)); U <- mean(unlist(vars))
  B <- var(unlist(betas)); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, fmi=(1+1/m)*B/Tm)
}
grab_pool <- function(bs, vs){
  rb <- rubin(bs, vs); q <- qt(.975, rb$df)
  list(HR=exp(rb$beta), lo=exp(rb$beta-q*rb$se), hi=exp(rb$beta+q*rb$se),
       P=2*pt(abs(rb$beta/rb$se), df=rb$df, lower.tail=FALSE), fmi=rb$fmi, m=length(bs))
}

# ---------- D1) 365d Model B + covariate tt terms ----------
cat("== D1) 365d Model B + tt(age,lactate,SOFA) ==\n")
tm <- "t_lm_365"; ev <- "event_lm_365"
fmlD1 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml, " + ", tt_fml))
bs <- c(); vs <- c()
for (i in 1:50) {
  dd <- complete(imp, i); dd$gv10 <- dd$gv/10
  fit <- tryCatch(coxph(fmlD1, data=dd, tt=function(x,t,...) x*log(pmax(t,0.5))), error=function(e) NULL)
  if (is.null(fit)) { cat("imputation", i, "D1 fit failed\n"); next }
  si <- summary(fit)$coefficients
  bs <- c(bs, si["gv10","coef"]); vs <- c(vs, si["gv10","se(coef)"]^2)
}
r1 <- grab_pool(bs, vs)
tabD1 <- data.frame(model_id="MICE_B_365d_covariate_tt_sensitivity",
  term="tt(age_at_admission)+tt(lactate_postop_first)+tt(sofa_24h)",
  N=10561, events=745, HR_per10=r1$HR, lo=r1$lo, hi=r1$hi, P=r1$P, fmi=r1$fmi, n_imputations=r1$m,
  note="365-day Model B with prespecified time-varying log(time) terms for PH-violating covariates age, lactate, SOFA; compare uncorrected average HR 0.985 (0.942-1.030) P=0.506")
write.csv(tabD1, file.path(ROOT,"results","PH365_COVARIATE_TT_SENSITIVITY_MICE_POOLED.csv"), row.names=FALSE)
cat(sprintf("D1 GV: HR=%.4f (%.4f-%.4f) P=%.4f fmi=%.4f m=%d\n", r1$HR, r1$lo, r1$hi, r1$P, r1$fmi, r1$m))

# ---------- D2) 365d interval model + covariate tt terms ----------
cat("== D2) 365d intervals + tt(age,lactate,SOFA) ==\n")
INTS <- c("d1-7","d8-30","d31-365")
pool_int <- setNames(vector("list",3), INTS)
for (lv in INTS) pool_int[[lv]] <- list(b=c(), v=c())
fmlD2 <- as.formula(paste0("Surv(tstart,tstop,evsplit) ~ gv10:interval + ", rcs_mean, " + ", covs_fml, " + ", tt_fml))
for (i in 1:50) {
  dd <- complete(imp, i); dd$gv10 <- dd$gv/10
  sp <- survSplit(as.formula(paste0("Surv(",tm,",",ev,") ~ .")),
                  data=dd[, c("stay_id",tm,ev,"gv10","mean_glu",COVS)],
                  cut=c(7,30), zero=-0.5, episode="interval", start="tstart", end="tstop", event="evsplit")
  sp$interval <- factor(sp$interval, labels=INTS)
  sp <- sp[!is.na(sp$evsplit),]
  fit <- tryCatch(coxph(fmlD2, data=sp, tt=function(x,t,...) x*log(pmax(t,0.5))), error=function(e) NULL)
  if (is.null(fit)) { cat("imputation", i, "D2 fit failed\n"); next }
  si <- summary(fit)$coefficients
  for (lv in INTS) {
    term <- paste0("gv10:interval", lv)
    pool_int[[lv]]$b <- c(pool_int[[lv]]$b, si[term,"coef"])
    pool_int[[lv]]$v <- c(pool_int[[lv]]$v, si[term,"se(coef)"]^2)
  }
}
rowsD2 <- list()
for (lv in INTS) {
  r <- grab_pool(pool_int[[lv]]$b, pool_int[[lv]]$v)
  rowsD2[[lv]] <- data.frame(model_id="MICE_B_365d_intervals_covariate_tt", interval=lv,
    HR_per10=r$HR, lo=r$lo, hi=r$hi, P=r$P, fmi=r$fmi, n_imputations=r$m, stringsAsFactors=FALSE)
  cat(sprintf("D2 %s: HR=%.4f (%.4f-%.4f) P=%.4f fmi=%.4f\n", lv, r$HR, r$lo, r$hi, r$P, r$fmi))
}
tabD2 <- do.call(rbind, rowsD2)
write.csv(tabD2, file.path(ROOT,"results","PH365_INTERVAL_COVARIATE_TT_MICE_POOLED.csv"), row.names=FALSE)
cat("done\n")
sink()
