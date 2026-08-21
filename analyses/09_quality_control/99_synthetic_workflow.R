#!/usr/bin/env Rscript
# 99_synthetic_workflow.R — 最小完整流程:血糖清洗 → 特征 → Cox → modified Poisson → 图。
# 仅证明代码可执行;输出不构成研究估计,并写入 results/machine_readable/synthetic/。
rm(list=ls()); options(stringsAsFactors=FALSE)
suppressMessages({library(survival); library(sandwich); library(lmtest); library(ggplot2)})
ROOT <- normalizePath(getwd())
if (!file.exists(file.path(ROOT, "data/synthetic/synthetic_cohort.csv"))) {
  if (file.exists("../data/synthetic/synthetic_cohort.csv")) ROOT <- normalizePath("..")
  else if (file.exists("../../data/synthetic/synthetic_cohort.csv")) ROOT <- normalizePath("../..")
}
OUT <- file.path(ROOT, "results", "machine_readable", "synthetic")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20260726L)

glu <- read.csv(file.path(ROOT, "data/synthetic/synthetic_glucose_long.csv"), stringsAsFactors=FALSE)
co  <- read.csv(file.path(ROOT, "data/synthetic/synthetic_cohort.csv"), stringsAsFactors=FALSE)
stopifnot(all(co$SYNTHETIC %in% c(TRUE,"TRUE","True","true",1,"1")))

# ---- 1. 清洗:完全重复剔除 → 同分钟来源优先级 → 中位数 ----
prio <- c(central_lab=1, blood_gas=2, poct=3, icu_charted=4)
glu <- glu[!duplicated(glu[,c("subject_id","minute","value","source_class")]),]
glu$pr <- prio[glu$source_class]
glu <- glu[order(glu$subject_id, glu$minute, glu$pr),]
best <- ave(glu$pr, glu$subject_id, glu$minute, FUN=min)
gm <- aggregate(value ~ subject_id + minute, data=glu[glu$pr==best,], FUN=median)
feat <- do.call(rbind, lapply(split(gm, gm$subject_id), function(g) {
  g <- g[order(g$minute),]; v <- g$value; n <- length(v)
  data.frame(subject_id=g$subject_id[1], n=n,
    mean_glucose=mean(v),
    gv_sd=if(n>=2) sd(v) else NA_real_,
    span_hours=(max(g$minute)-min(g$minute))/60)
}))
d <- merge(co, feat, by="subject_id")
d <- d[!is.na(d$gv_sd) & d$n>=2,]
stopifnot(nrow(d) > 50)
d$gv10 <- d$gv_sd/10
K <- list(gv_mean=mean(d$gv_sd), gv_sd=sd(d$gv_sd))
d$gv_z <- (d$gv_sd - K$gv_mean)/K$gv_sd
d$mean_z <- (d$mean_glucose - mean(d$mean_glucose))/sd(d$mean_glucose)

# ---- 2. 简化 Cox(I2 形式:GV + mean + age) ----
f1 <- coxph(Surv(t_lm_30, event_lm_30) ~ gv10 + mean_z + age, data=d)
s1 <- summary(f1)
cox_out <- data.frame(model="synthetic_cox_I2form", N=nrow(d), events=sum(d$event_lm_30),
  HR_per10=s1$coefficients["gv10","exp(coef)"],
  lo=s1$conf.int["gv10","lower .95"], hi=s1$conf.int["gv10","upper .95"],
  P=s1$coefficients["gv10","Pr(>|z|)"])

# ---- 3. modified Poisson + 稳健方差 ----
f2 <- glm(event_lm_30 ~ gv10 + mean_z + age, data=d, family=poisson(link="log"))
V <- vcovHC(f2, type="HC1"); ct <- coeftest(f2, vcov.=V)
mp_out <- data.frame(model="synthetic_modified_poisson", N=nrow(d), events=sum(d$event_lm_30),
  RR_per10=exp(ct["gv10","Estimate"]),
  lo=exp(ct["gv10","Estimate"]-1.96*ct["gv10","Std. Error"]),
  hi=exp(ct["gv10","Estimate"]+1.96*ct["gv10","Std. Error"]),
  P=ct["gv10","Pr(>|z|)"])

# ---- 4. 图 ----
p <- ggplot(d, aes(x=gv_sd)) + geom_histogram(bins=40, fill="#0072B2", alpha=0.8) +
  labs(title="Synthetic GV distribution (NOT study data)", x="GV (mg/dL)", y="patients") + theme_bw()
ggsave(file.path(OUT, "synthetic_gv_hist.png"), p, width=6, height=4, dpi=300)

cox2 <- data.frame(model=cox_out$model, N=cox_out$N, events=cox_out$events,
  effect_label="HR_per10", effect=cox_out$HR_per10, lo=cox_out$lo, hi=cox_out$hi, P=cox_out$P)
mp2 <- data.frame(model=mp_out$model, N=mp_out$N, events=mp_out$events,
  effect_label="RR_per10", effect=mp_out$RR_per10, lo=mp_out$lo, hi=mp_out$hi, P=mp_out$P)
write.csv(rbind(cox2, mp2), file.path(OUT, "synthetic_workflow_results.csv"), row.names=FALSE)
write.csv(feat, file.path(OUT, "synthetic_features.csv"), row.names=FALSE)
cat("SYNTHETIC_WORKFLOW_OK N=", nrow(d), " events=", sum(d$event_lm_30), "\n", sep="")
