# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 17_ccmask_fix_rerun.R — 修复 complete-case mask 范围错误(版本 v3)
# 根因:RECON_I2_30d 的 CC mask 含模型外变量(asa/bmi/sex/charlson 等),
#       误排 10 人(7 ASA-missing + 3 BMI-missing)及其中 4 例事件。
# 30d 模型公式仅用 gv + RCS(mean) + age(+I3 的 log_count/span),frame 内全部完整
# → 正确 CC = 全 frame 1,353/27,无需 MICE;
# 365d 全协变量模型按冻结协议使用 MICE m=50 为主,CC 为敏感性对照。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite); library(rms)})
ROOT <- PGV("inspire_work")
sink(file.path(ROOT,"logs","17_ccmask_fix.log"), split=TRUE)

r30 <- read.csv(file.path(ROOT,"data","outcome_30d_reconciled.csv"), stringsAsFactors=FALSE)
cm <- read.csv(file.path(ROOT,"data","comorbidity.csv"), stringsAsFactors=FALSE)
b  <- read.csv(file.path(ROOT,"data","inspire_base.csv"), stringsAsFactors=FALSE)
d <- merge(r30[,c("subject_id","op_id","landmark_24h","day30_cutoff","death_time_composite","event_30d_reconciled")],
           b, by=c("subject_id","op_id"))
d <- merge(d, cm, by="subject_id", all.x=TRUE)
d <- d[is.na(d$death_time_composite) | d$death_time_composite > d$landmark_24h,]
d$gv <- d$gv_sd; d$mean_glu <- d$mean_glucose; d$n_gv <- d$n_glucose_0_24h
d$bmi <- d$weight/(d$height/100)^2
d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$sex <- factor(d$sex); d$asa_f <- factor(d$asa); d$emop_f <- as.integer(d$emop %in% c("1",1))
d$surgery_group <- factor(d$surgery_group)
d$event_30d <- as.integer(d$event_30d_reconciled %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$event_365d <- as.integer(!is.na(d$death_time_composite) & d$death_time_composite > d$landmark_24h &
                           d$death_time_composite <= d$opend_time + 365*1440)
d$fu_end <- pmin(ifelse(is.na(d$death_time_composite), Inf, d$death_time_composite), d$discharge_time)
d$t30 <- (pmin(d$fu_end, d$day30_cutoff) - d$landmark_24h)/1440.0
d$t365 <- (pmin(d$fu_end, d$opend_time + 365*1440) - d$landmark_24h)/1440.0
stopifnot(all(d$t30 >= 0))
d$gv10 <- d$gv/10; d$log_count <- log(d$n_gv)
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
kSD <- K$gv_sd/10
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
cat("frame N =", nrow(d), "; 30d 事件 =", sum(d$event_30d), "; 365d 事件 =", sum(d$event_365d), "\n")

# ---- 正确的 CC mask:仅模型实际变量 ----
mask30 <- c("t30","event_30d","gv","gv10","mean_glu","age","log_count","span_hours")
stopifnot(all(colSums(is.na(d[,mask30]))==0))
cc30 <- d[complete.cases(d[,mask30]) & d$t30>0,]
stopifnot(nrow(cc30)==1353 && sum(cc30$event_30d)==27)
rep <- function(fit, model_id, model, N, events, note){
  s <- summary(fit)
  data.frame(model_id=model_id, frame="INSPIRE_OPEND_24H_LANDMARK_V2",
    outcome_version="30d_reconciled_composite_v2 + ccmask_v3", model=model, N=N, events=events,
    HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
    hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"],
    HR_perSD=exp(s$coefficients["gv10","coef"]*kSD),
    ph_global=tryCatch(cox.zph(fit)$table["GLOBAL","p"], error=function(e) NA),
    note=note, stringsAsFactors=FALSE)
}
rows30 <- list(
  rep(coxph(Surv(t30,event_30d) ~ gv10 + age, data=cc30, x=TRUE, y=TRUE),
      "CCMASKFIX_I1_30d","Model I1", nrow(cc30), sum(cc30$event_30d),
      "corrected CC mask (model variables only); N=1353, events=27"),
  rep(coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean, " + age")), data=cc30, x=TRUE, y=TRUE),
      "CCMASKFIX_I2_30d","Model I2 (primary)", nrow(cc30), sum(cc30$event_30d),
      "corrected CC mask (model variables only); no MICE needed (no missing model variables)"),
  rep(coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean, " + age + log_count + span_hours")), data=cc30, x=TRUE, y=TRUE),
      "CCMASKFIX_I3_30d","Model I3", nrow(cc30), sum(cc30$event_30d),
      "corrected CC mask (model variables only)")
)
tab30 <- do.call(rbind, rows30)
write.csv(tab30, file.path(ROOT,"results","27_corrected_30d_results_v3.csv"), row.names=FALSE)
print(tab30[,c("model_id","N","events","HR_per10","lo","hi","P","ph_global")])

# ---- 30 日绝对风险(全 frame 1,353) ----
risk_c <- function(dd){
  fit <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean, " + age")), data=dd)
  bh <- basehaz(fit, centered=FALSE); bh <- bh[bh$time<=29,]; H <- max(bh$hazard)
  lp25 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q25; nd$gv10 <- nd$gv/10; nd}, type="lp", reference="zero")
  lp75 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q75; nd$gv10 <- nd$gv/10; nd}, type="lp", reference="zero")
  r25 <- mean(1 - exp(-H*exp(lp25))); r75 <- mean(1 - exp(-H*exp(lp75)))
  c(rd=r75-r25, rr=r75/r25, r25=r25, r75=r75)
}
pt <- risk_c(cc30)
B <- 1000; vals <- matrix(NA, B, 4); fails <- 0L; n <- nrow(cc30)
for (b_ in 1:B) {
  idx <- sample.int(n, n, replace=TRUE)
  rb <- tryCatch(risk_c(cc30[idx,]), error=function(e) NULL)
  if (is.null(rb) || any(!is.finite(rb))) { fails <- fails+1L; next }
  vals[b_,] <- rb
}
colnames(vals) <- c("rd","rr","r25","r75")
ci <- apply(vals, 2, quantile, probs=c(.025,.975), na.rm=TRUE)
ab30 <- data.frame(model_id="CCMASKFIX_ABSRISK_30d", horizon=30, N=nrow(cc30), events=sum(cc30$event_30d),
  gv_q25=K$gv_q25, gv_q75=K$gv_q75, risk_q25=pt["r25"], risk_q75=pt["r75"],
  rd=pt["rd"], rd_lo=ci[1,"rd"], rd_hi=ci[2,"rd"], rr=pt["rr"], rr_lo=ci[1,"rr"], rr_hi=ci[2,"rr"],
  n_boot=B, n_boot_failed=fails, note="corrected CC mask; full frame 1353/27", stringsAsFactors=FALSE)
write.csv(ab30, file.path(ROOT,"results","27b_corrected_30d_absrisk_v3.csv"), row.names=FALSE)
print(ab30)

# ---- 365 日:MICE m=50(协议主方案)+ CC 对照 ----
CLIN365 <- c("age","sex","diabetes","charlson_without_diabetes","asa_f","emop_f","surgery_group","bmi")
covs365 <- paste(CLIN365, collapse=" + ")
cc365 <- d[complete.cases(d[,c("t365","event_365d","gv","gv10","mean_glu",CLIN365)]) & d$t365>0,]
cat("365d CC N =", nrow(cc365), "; 事件 =", sum(cc365$event_365d), "\n")
f365_cc <- coxph(as.formula(paste0("Surv(t365,event_365d) ~ gv10 + ", rcs_mean, " + ", covs365)), data=cc365, x=TRUE, y=TRUE)
s_cc <- summary(f365_cc)

mi_data <- d[, c("subject_id","t365","event_365d","gv","gv10","mean_glu","n_gv","span_hours","log_count",
                 "age","sex","diabetes","charlson_without_diabetes","asa_f","emop_f","surgery_group","bmi")]
mi_data$sex <- as.numeric(mi_data$sex); mi_data$surgery_group <- as.numeric(mi_data$surgery_group)
mi_data$asa_num <- suppressWarnings(as.numeric(as.character(mi_data$asa_f)))
mi_data <- mi_data[, names(mi_data)!="asa_f"]
mi_data$H365 <- nelsonaalen(mi_data, "t365", "event_365d")
pred <- make.predictorMatrix(mi_data); pred[,] <- 0L
meth <- make.method(mi_data); meth[] <- ""
tgt_imp <- intersect(c("bmi","asa_num"), names(mi_data))
tgt_imp <- tgt_imp[colSums(is.na(mi_data[,tgt_imp]))>0]
impute_preds <- setdiff(names(mi_data), "subject_id")
for (v in tgt_imp) { meth[v] <- "pmm"; pred[v, impute_preds] <- 1L }
cat("MICE 插补变量:", paste(tgt_imp, collapse=", "), "\n")
imp <- mice(mi_data, m=50, maxit=20, method=meth, predictorMatrix=pred, seed=SEED, printFlag=FALSE)
saveRDS(imp, file.path(ROOT,"results","inspire_mice_m50_365d_v3.rds"))
lg <- imp$loggedEvents
cat("loggedEvents:", if(is.null(lg)) 0 else nrow(lg), "\n")

rubin <- function(bs, vs){
  m <- length(bs); qb <- mean(bs); U <- mean(vs); B <- var(bs); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, p=2*pt(abs(qb/sqrt(Tm)), df=df, lower.tail=FALSE), fmi=(1+1/m)*B/Tm)
}
bs <- c(); vs <- c()
for (i in 1:50) {
  dd <- complete(imp, i)
  dd$sex <- factor(dd$sex); dd$surgery_group <- factor(dd$surgery_group, levels=c(1,2), labels=c("OPEN_CABG","OPEN_VALVE"))
  dd$asa_f <- factor(dd$asa_num)
  f <- tryCatch(coxph(as.formula(paste0("Surv(t365,event_365d) ~ gv10 + ", rcs_mean, " + ", covs365)), data=dd), error=function(e) NULL)
  if (is.null(f)) next
  s <- summary(f)
  bs <- c(bs, s$coefficients["gv10","coef"]); vs <- c(vs, s$coefficients["gv10","se(coef)"]^2)
}
rb <- rubin(bs, vs)
tab365 <- data.frame(
  model_id=c("CC_I2_365d","MICE_I2_365d"), method=c("complete-case","MICE m=50 (protocol primary)"),
  N=c(nrow(cc365), nrow(d)), events=c(sum(cc365$event_365d), sum(d$event_365d)),
  HR_per10=c(s_cc$coefficients["gv10","exp(coef)"], exp(rb$beta)),
  lo=c(s_cc$conf.int["gv10","lower .95"], exp(rb$beta-qt(.975,rb$df)*rb$se)),
  hi=c(s_cc$conf.int["gv10","upper .95"], exp(rb$beta+qt(.975,rb$df)*rb$se)),
  P=c(s_cc$coefficients["gv10","Pr(>|z|)"], rb$p),
  HR_perSD=c(exp(s_cc$coefficients["gv10","coef"]*kSD), exp(rb$beta*kSD)),
  fmi=c(NA, rb$fmi), frame="INSPIRE_OPEND_24H_LANDMARK_V2",
  outcome_version="365d_exploratory_registry_linked_v2 + ccmask_v3",
  note=c("CC sensitivity (ASA/BMI complete subset)","MICE primary; impute bmi/asa_num; 365d exploratory"),
  stringsAsFactors=FALSE)
write.csv(tab365, file.path(ROOT,"results","28_corrected_365d_results_v3.csv"), row.names=FALSE)
print(tab365[,c("model_id","method","N","events","HR_per10","lo","hi","P")])
cat("PHASE17_DONE\n")
sink()
