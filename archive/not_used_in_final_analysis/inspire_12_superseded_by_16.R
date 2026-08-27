#!/usr/bin/env Rscript
# 12_inspire_sensitivity.R — 时间锚点敏感性(完全相同样本)与量化分辨率敏感性(协议 §7/§8)
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.libPaths(c(file.path("private", "r-library"), .libPaths()))
suppressMessages({library(survival); library(jsonlite); library(rms)})
ROOT <- normalizePath(file.path("private", "inspire_project"), mustWork=FALSE)
sink(file.path(ROOT,"logs","12_inspire_sensitivity.log"), split=TRUE)

b <- read.csv(file.path(ROOT,"data","inspire_base.csv"), stringsAsFactors=FALSE)
cm <- read.csv(file.path(ROOT,"data","comorbidity.csv"), stringsAsFactors=FALSE)
ga <- read.csv(file.path(ROOT,"data","glucose_features_anchor.csv"), stringsAsFactors=FALSE)
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
d <- merge(b, cm, by="subject_id", all.x=TRUE)
d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$sex <- factor(d$sex)
d$event_30d <- as.integer(!is.na(d$allcause_death_time) & d$allcause_death_time <= d$opend_time + 30*1440)
d$event_365d <- as.integer(!is.na(d$allcause_death_time) & d$allcause_death_time <= d$opend_time + 365*1440)
d$followup_days <- (pmin(ifelse(is.na(d$allcause_death_time), Inf, d$allcause_death_time), d$discharge_time) - d$opend_time)/1440
d$t30 <- pmin(d$followup_days, 30) - 1; d$t365 <- pmin(d$followup_days, 365) - 1
CLIN_365 <- c("age","sex","diabetes","charlson_without_diabetes")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")

rows <- list()
add <- function(...) { rows[[length(rows)+1]] <<- data.frame(..., stringsAsFactors=FALSE) }

# ---- 1. 时间锚点敏感性(完全相同样本:该锚点窗内有 >=2 测量且 landmark 合格者) ----
anchors <- c("opend_0_24h","opstart_0_24h","orout_0_24h","icuin_0_24h","opend_0_12h","opend_0_48h")
anchor_stats <- list()
for (a in anchors) {
  f <- ga[ga$anchor==a & ga$n>=2,]
  dd <- d[d$subject_id %in% f$subject_id & d$t30>=0 &
          d$landmark_eligible_allcause %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
  dd <- merge(dd, f[,c("subject_id","n","mean_glucose","gv_sd")], by="subject_id", suffixes=c("","_aw"))
  dd$gv10 <- dd$gv_sd_aw/10; dd$mean_glu <- dd$mean_glucose_aw
  anchor_stats[[a]] <- data.frame(anchor=a, n_patients=nrow(dd),
    events_30d=sum(dd$event_30d), events_365d=sum(dd$event_365d),
    mean_of_mean_glucose=mean(dd$mean_glu), mean_gv_sd=mean(dd$gv_sd_aw),
    cor_gv_mean=cor(dd$gv_sd_aw, dd$mean_glu),
    gv_q25=unname(quantile(dd$gv_sd_aw,.25)), gv_q75=unname(quantile(dd$gv_sd_aw,.75)))
  for (hz in c(30,365)) {
    tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
    dd2 <- dd[complete.cases(dd[,c(tm,ev,"gv10","mean_glu",CLIN_365)]),]
    if (nrow(dd2) < 100 || sum(dd2[[ev]]) < 10) {
      add(model_id=paste0("ANCH_",a,"_",hz,"d"), anchor=a, N=nrow(dd2), events=sum(dd2[[ev]]),
        note="insufficient N/events (<100 N or <10 events)")
      next
    }
    fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", paste(CLIN_365, collapse=' + '))), data=dd2)
    s <- summary(fit)
    add(model_id=paste0("ANCH_",a,"_",hz,"d"), anchor=a, N=nrow(dd2), events=sum(dd2[[ev]]),
      HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
      hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"],
      note="Model I2 form (RCS mean + reduced clinical); NOT selected by min P")
  }
}
anch_tab <- do.call(rbind, rows)
write.csv(anch_tab, file.path(ROOT,"results","05_time_anchor_results.csv"), row.names=FALSE)
write.csv(do.call(rbind, anchor_stats), file.path(ROOT,"results","05b_anchor_descriptives.csv"), row.names=FALSE)
print(anch_tab)

# ---- 2. 量化分辨率敏感性(同一 opend 主窗;阈值与指标;全部使用 inspire_base 自带列) ----
tgt <- d[d$age>=18 & !is.na(d$opend_time) & d$t30>=0 &
         d$landmark_eligible_allcause %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
tgt$n_gv <- tgt$n_glucose_0_24h
tgt$mean_glu <- tgt$mean_glucose
res_rows <- list()
run_res <- function(dd, xvar, scale_div, tag, label){
  dd$x <- dd[[xvar]]/scale_div
  for (hz in c(30,365)) {
    tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
    dd2 <- dd[complete.cases(dd[,c(tm,ev,"x","mean_glu",CLIN_365)]),]
    if (nrow(dd2) < 100 || sum(dd2[[ev]]) < 10) {
      res_rows[[length(res_rows)+1]] <<- data.frame(model_id=paste0("RES_",tag,"_",hz,"d"),
        cohort=label, N=nrow(dd2), events=sum(dd2[[ev]]), note="insufficient N/events", stringsAsFactors=FALSE)
      next
    }
    fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ x + ", rcs_mean, " + ", paste(CLIN_365, collapse=' + '))), data=dd2)
    s <- summary(fit)
    res_rows[[length(res_rows)+1]] <<- data.frame(model_id=paste0("RES_",tag,"_",hz,"d"),
      cohort=label, metric=xvar, scale=paste0("per ", scale_div, " units"),
      N=nrow(dd2), events=sum(dd2[[ev]]),
      HR=s$coefficients["x","exp(coef)"], lo=s$conf.int["x","lower .95"],
      hi=s$conf.int["x","upper .95"], P=s$coefficients["x","Pr(>|z|)"], note="", stringsAsFactors=FALSE)
  }
}
# 阈值组
for (thr in list(list(name="ge2", mask=tgt$n_gv>=2),
                 list(name="ge3tp", mask=tgt$n_gv>=3),
                 list(name="ge3tp_ge3val", mask=tgt$n_gv>=3 & tgt$n_distinct_values>=3))) {
  dd <- tgt[thr$mask,]
  for (m in list(list(v="gv_sd", div=10, lbl="SD per 10 mg/dL"),
                 list(v="arv", div=10, lbl="ARV per 10 mg/dL"),
                 list(v="mad_glucose", div=10, lbl="MAD per 10 mg/dL"),
                 list(v="iqr_glucose", div=10, lbl="IQR per 10 mg/dL")))
    run_res(dd, m$v, m$div, paste0(thr$name,"_",m$v), paste0(thr$name, " (", m$lbl, ")"))
}
res_tab <- do.call(rbind, res_rows)
write.csv(res_tab, file.path(ROOT,"results","12_sensitivity_results.csv"), row.names=FALSE)
print(res_tab[res_tab$metric=="gv_sd" | is.na(res_tab$metric),], max=40)
cat("PHASE12_DONE\n")
sink()
