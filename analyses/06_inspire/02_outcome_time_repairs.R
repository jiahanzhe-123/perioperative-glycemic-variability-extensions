# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 15_outcome_time_repairs_models.R — 修复后模型重跑
#  A) 30d 复合结局(30d_reconciled_composite_v2, 27 事件)全依赖模型重跑;
#  C) 正确 48h landmark(INSPIRE_OPEND_48H_LANDMARK_V2)分析 + 24h/48h 对照。
#  365d 按决策 B 仅保留 exploratory HR(不重跑绝对风险/RMST/校准)。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite); library(rms)})
ROOT <- PGV("inspire_work")
dir.create(file.path(ROOT,"results"), showWarnings=FALSE, recursive=TRUE)
sink(file.path(ROOT,"logs","15_repairs_models.log"), split=TRUE)

psql <- function(sql){
  # 由 bash 预导出到 data/;此处按 SQL 片段映射文件名
  if (grepl("outcome_30d_reconciled", sql)) return(read.csv(file.path(ROOT,"data","outcome_30d_reconciled.csv"), stringsAsFactors=FALSE))
  if (grepl("cohort_48h_landmark", sql)) return(read.csv(file.path(ROOT,"data","cohort_48h_landmark.csv"), stringsAsFactors=FALSE))
  stop("no mapping for SQL")
}
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))

# ============ A. 30d 复合结局重跑 ============
r30 <- psql("SELECT * FROM derived.outcome_30d_reconciled")
cm <- read.csv(file.path(ROOT,"data","comorbidity.csv"), stringsAsFactors=FALSE)
b  <- read.csv(file.path(ROOT,"data","inspire_base.csv"), stringsAsFactors=FALSE)
d <- merge(r30[,c("subject_id","op_id","landmark_24h","day30_cutoff","death_time_composite","event_30d_reconciled")],
           b, by=c("subject_id","op_id"))
d <- merge(d, cm, by="subject_id", all.x=TRUE)
# 复合结局下的 landmark 资格:composite 死亡时间必须晚于 landmark(2 例 EHR 死亡早于 landmark → 排除)
d <- d[is.na(d$death_time_composite) | d$death_time_composite > d$landmark_24h,]
d$gv <- d$gv_sd; d$mean_glu <- d$mean_glucose; d$n_gv <- d$n_glucose_0_24h
d$bmi <- d$weight/(d$height/100)^2
d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$sex <- factor(d$sex); d$asa_f <- factor(d$asa); d$emop_f <- as.integer(d$emop %in% c("1",1))
d$surgery_group <- factor(d$surgery_group)
# 复合结局时间:自 landmark 起至 min(death, last_observed, day30_cutoff)
d$fu_end <- pmin(ifelse(is.na(d$death_time_composite), Inf, d$death_time_composite), d$discharge_time)
d$t30 <- (pmin(d$fu_end, d$day30_cutoff) - d$landmark_24h)/1440.0
d$event_30d <- as.integer(d$event_30d_reconciled %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$gv10 <- d$gv/10
d$log_count <- log(d$n_gv)
kSD <- K$gv_sd/10
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
stopifnot(all(d$t30 >= 0), sum(d$event_30d) == 27)
cat("30d 复合队列 N =", nrow(d), "; 事件 =", sum(d$event_30d), "\n")

run30 <- function(dd, label){
  ev <- "event_30d"; tm <- "t30"
  fI1 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + age"))
  fI2 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + age"))
  fI3 <- as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + age + log_count + span_hours"))
  out <- list()
  for (mk in c("I1","I2","I3")) {
    fit <- coxph(get(paste0("f",mk)), data=dd, x=TRUE, y=TRUE)
    s <- summary(fit)
    out[[mk]] <- data.frame(model_id=paste0("RECON_",mk,"_30d"), frame="INSPIRE_OPEND_24H_LANDMARK_V2",
      outcome_version="30d_reconciled_composite_v2", model=paste0("Model ",mk), N=nrow(dd), events=sum(dd[[ev]]),
      HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
      hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"],
      HR_perSD=exp(s$coefficients["gv10","coef"]*kSD),
      ph_global=tryCatch(cox.zph(fit)$table["GLOBAL","p"], error=function(e) NA),
      note=label, stringsAsFactors=FALSE)
  }
  out
}
cc <- d[complete.cases(d[,c("t30","event_30d","gv","gv10","mean_glu","age","sex","diabetes",
  "charlson_without_diabetes","asa_f","emop_f","surgery_group","bmi","log_count","span_hours")]) & d$t30>0,]
res30 <- do.call(rbind, run30(cc, "complete-case; core covariates per frozen EPV rule"))
write.csv(res30, file.path(ROOT,"results","20_corrected_30d_results.csv"), row.names=FALSE)
print(res30[,c("model_id","N","events","HR_per10","lo","hi","P")])

# 绝对风险(30d 复合,Q75/Q25,1000 bootstrap)
risk_c <- function(dd){
  fit <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean, " + age")), data=dd)
  bh <- basehaz(fit, centered=FALSE); bh <- bh[bh$time<=29,]; H <- max(bh$hazard)
  lp25 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q25; nd$gv10 <- nd$gv/10; nd}, type="lp", reference="zero")
  lp75 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q75; nd$gv10 <- nd$gv/10; nd}, type="lp", reference="zero")
  r25 <- mean(1 - exp(-H*exp(lp25))); r75 <- mean(1 - exp(-H*exp(lp75)))
  c(rd=r75-r25, rr=r75/r25, r25=r25, r75=r75)
}
pt <- risk_c(cc)
B <- 1000; vals <- matrix(NA, B, 4); fails <- 0L; n <- nrow(cc)
for (b_ in 1:B) {
  idx <- sample.int(n, n, replace=TRUE)
  rb <- tryCatch(risk_c(cc[idx,]), error=function(e) NULL)
  if (is.null(rb) || any(!is.finite(rb))) { fails <- fails+1L; next }
  vals[b_,] <- rb
}
colnames(vals) <- c("rd","rr","r25","r75")
ci <- apply(vals, 2, quantile, probs=c(.025,.975), na.rm=TRUE)
ab30 <- data.frame(model_id="RECON_ABSRISK_30d", horizon=30, N=nrow(cc), events=sum(cc$event_30d),
  gv_q25=K$gv_q25, gv_q75=K$gv_q75, risk_q25=pt["r25"], risk_q75=pt["r75"],
  rd=pt["rd"], rd_lo=ci[1,"rd"], rd_hi=ci[2,"rd"], rr=pt["rr"], rr_lo=ci[1,"rr"], rr_hi=ci[2,"rr"],
  n_boot=B, n_boot_failed=fails, note="reconciled 30d composite", stringsAsFactors=FALSE)
write.csv(ab30, file.path(ROOT,"results","20b_corrected_30d_absrisk.csv"), row.names=FALSE)
print(ab30)

# ============ C. 正确 48h landmark ============
c48 <- psql("SELECT * FROM derived.cohort_48h_landmark")
d48 <- merge(c48, cm, by="subject_id", all.x=TRUE)
d48$gv <- d48$gv_sd; d48$mean_glu <- d48$mean_glucose
d48$bmi <- d48$weight/(d48$height/100)^2
d48$diabetes <- as.integer(d48$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d48$sex <- factor(d48$sex)
# 48h 队列自有常数(一次冻结,不沿用 24h)
K48 <- list(gv_mean=mean(d48$gv), gv_sd=sd(d48$gv), gv_median=median(d48$gv),
  gv_q25=unname(quantile(d48$gv,.25)), gv_q75=unname(quantile(d48$gv,.75)),
  mean_glu_mean=mean(d48$mean_glu), mean_glu_sd=sd(d48$mean_glu),
  mean_glu_knots=unname(quantile(d48$mean_glu, c(.05,.35,.65,.95))),
  N=nrow(d48), events_30d=sum(d48$event_30d_48h %in% c(TRUE,"t","True","TRUE","true",1,"1")),
  frame="INSPIRE_OPEND_48H_LANDMARK_V2", seed=SEED)
write_json(K48, file.path(ROOT,"results","inspire_48h_constants.json"), pretty=TRUE, auto_unbox=TRUE)
d48$event_30d <- as.integer(d48$event_30d_48h %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d48$fu_end <- pmin(ifelse(is.na(d48$death_time_composite), Inf, d48$death_time_composite), d48$discharge_time)
d48$t30 <- (pmin(d48$fu_end, d48$day30_cutoff) - d48$landmark_48h)/1440.0
stopifnot(all(d48$t30 >= 0))
d48$gv10 <- d48$gv/10; d48$log_count <- log(d48$n_glucose_0_48h)
rcs_mean48 <- paste0("rms::rcs(mean_glu, c(", paste(format(K48$mean_glu_knots, digits=10), collapse=","), "))")
cc48 <- d48[complete.cases(d48[,c("t30","event_30d","gv","gv10","mean_glu","age","sex","diabetes")]) & d48$t30>0,]
cat("48h landmark 队列 N =", nrow(d48), "; 事件 =", sum(d48$event_30d), "; CC =", nrow(cc48), "\n")
f48 <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ gv10 + ", rcs_mean48, " + age")), data=cc48, x=TRUE, y=TRUE)
s48 <- summary(f48)
r48 <- data.frame(model_id="CORR_48H_I2_30d", frame="INSPIRE_OPEND_48H_LANDMARK_V2",
  outcome_version="30d_reconciled_composite_v2", N=nrow(cc48), events=sum(cc48$event_30d),
  HR_per10=s48$coefficients["gv10","exp(coef)"], lo=s48$conf.int["gv10","lower .95"],
  hi=s48$conf.int["gv10","upper .95"], P=s48$coefficients["gv10","Pr(>|z|)"],
  HR_perSD=exp(s48$coefficients["gv10","coef"]*K48$gv_sd/10),
  ph_global=tryCatch(cox.zph(f48)$table["GLOBAL","p"], error=function(e) NA),
  note="correct 48h landmark; exposure strictly < 48h landmark; risk starts at 48h", stringsAsFactors=FALSE)
write.csv(r48, file.path(ROOT,"results","23_corrected_48h_landmark_results.csv"), row.names=FALSE)
print(r48)

# ============ 24h vs 48h 对照 ============
j48 <- merge(d[,c("subject_id","gv","mean_glu","n_gv")], d48[,c("subject_id","gv","mean_glu","n_glucose_0_48h")],
             by="subject_id", suffixes=c("_24","_48"))
cmp <- data.frame(
  item=c("N_24h_cohort","events_24h","N_48h_cohort","events_48h","N_intersection",
         "deaths_between_24h_and_48h_landmarks",
         "cor_gv_24h_vs_48h","cor_mean_24h_vs_48h",
         "mean_n_glucose_24h","mean_n_glucose_48h","mean_span_hours_24h","mean_span_hours_48h",
         "reclassified_into_Q4_from_le_Q3_pct"),
  value=c(nrow(d), sum(d$event_30d), nrow(d48), sum(d48$event_30d), nrow(j48),
          sum((d$death_time_composite > d$landmark_24h) & (d$death_time_composite <= d$landmark_24h + 1440), na.rm=TRUE),
          cor(j48$gv_24, j48$gv_48), cor(j48$mean_glu_24, j48$mean_glu_48),
          mean(d$n_gv), mean(d48$n_glucose_0_48h),
          mean(d$span_hours, na.rm=TRUE), mean(d48$span_hours, na.rm=TRUE),
          NA))
cmp$value[13] <- mean(j48$gv_24 > quantile(d$gv,.75) & j48$gv_48 > quantile(d48$gv,.75), na.rm=TRUE)
write.csv(cmp, file.path(ROOT,"results","22_24h_vs_48h_cohort_comparison.csv"), row.names=FALSE)

# 测量过程审计
mpa <- data.frame(
  item=c("records in 24h window (shared)","records in 48h window (total)","records added by 24-48h extension",
         "patients gaining >=2 measurements only by 48h extension",
         "cor(n_measurements_24h, gv_48h)","cor(n_measurements_48h, gv_48h)"),
  value=c(sum(d$n_gv), sum(d48$n_glucose_0_48h), sum(d48$n_glucose_0_48h) - sum(d$n_gv[d$subject_id %in% d48$subject_id]),
          sum(!(d48$subject_id %in% d$subject_id)),
          cor(d$n_gv[d$subject_id %in% d48$subject_id], d48$gv[match(d$subject_id[d$subject_id %in% d48$subject_id], d48$subject_id)]),
          cor(d48$n_glucose_0_48h, d48$gv)))
write.csv(mpa, file.path(ROOT,"results","24_48h_measurement_process_audit.csv"), row.names=FALSE)
print(cmp); print(mpa)

# 图:24h vs 48h GV 散点 + 重分类
suppressMessages(library(ggplot2))
j48$reclass <- ifelse(j48$gv_48 > quantile(d48$gv,.75), "Q4(48h)",
                ifelse(j48$gv_24 > quantile(d$gv,.75), "Q4(24h only)", "not-Q4"))
p <- ggplot(j48, aes(gv_24, gv_48, color=reclass)) +
  geom_point(alpha=0.4, size=1) + geom_abline(slope=1, intercept=0, linetype=2, color="grey50") +
  scale_color_manual(values=c("Q4(48h)"="#D55E00","Q4(24h only)"="#0072B2","not-Q4"="grey70")) +
  labs(title="GV: 24h vs 48h window (same patients)", x="GV (0-24h)", y="GV (0-48h)") + theme_bw()
ggsave(file.path(ROOT,"figures","24h_vs_48h_gv_scatter.png"), p, width=6.5, height=5.5, dpi=600)
tab <- as.data.frame(table(j48$reclass)); names(tab)[1] <- "reclass"
p2 <- ggplot(tab, aes(reclass, Freq, fill=reclass)) + geom_col() +
  scale_fill_manual(values=c("Q4(48h)"="#D55E00","Q4(24h only)"="#0072B2","not-Q4"="grey70")) +
  labs(title="Quartile reclassification 24h vs 48h", x=NULL, y="patients") + theme_bw()
ggsave(file.path(ROOT,"figures","24h_vs_48h_reclassification.png"), p2, width=6, height=4, dpi=600)
cat("PHASE15_DONE\n")
sink()
