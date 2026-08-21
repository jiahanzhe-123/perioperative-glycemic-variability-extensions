# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 16_repairs_remaining_reruns.R — 修正结局下的分辨率/锚点/联合模块 30d 部分重跑
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
sink(file.path(ROOT,"logs","16_repairs_remaining.log"), split=TRUE)

r30 <- read.csv(file.path(ROOT,"data","outcome_30d_reconciled.csv"), stringsAsFactors=FALSE)
cm <- read.csv(file.path(ROOT,"data","comorbidity.csv"), stringsAsFactors=FALSE)
b  <- read.csv(file.path(ROOT,"data","inspire_base.csv"), stringsAsFactors=FALSE)
ga <- read.csv(file.path(ROOT,"data","glucose_features_anchor.csv"), stringsAsFactors=FALSE)
K  <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
d <- merge(r30[,c("subject_id","op_id","landmark_24h","day30_cutoff","death_time_composite","event_30d_reconciled")],
           b, by=c("subject_id","op_id"))
d <- merge(d, cm, by="subject_id", all.x=TRUE)
d <- d[is.na(d$death_time_composite) | d$death_time_composite > d$landmark_24h,]
d$gv <- d$gv_sd; d$mean_glu <- d$mean_glucose; d$n_gv <- d$n_glucose_0_24h
d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$sex <- factor(d$sex)
d$event_30d <- as.integer(d$event_30d_reconciled %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$fu_end <- pmin(ifelse(is.na(d$death_time_composite), Inf, d$death_time_composite), d$discharge_time)
d$t30 <- (pmin(d$fu_end, d$day30_cutoff) - d$landmark_24h)/1440.0
stopifnot(all(d$t30>=0))
d$gv10 <- d$gv/10
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
rows <- list()
add <- function(...) { rows[[length(rows)+1]] <<- data.frame(..., stringsAsFactors=FALSE) }
fit_rep <- function(dd, xvar, div, tag, label){
  dd$x <- dd[[xvar]]/div
  dd2 <- dd[complete.cases(dd[,c("t30","event_30d","x","mean_glu","age")]) & dd$t30>0,]
  if (nrow(dd2) < 100 || sum(dd2$event_30d) < 8) {
    add(model_id=paste0("RECON_",tag,"_30d"), cohort=label, N=nrow(dd2), events=sum(dd2$event_30d),
        note="not reliably estimable (<100 N or <8 events)")
    return(invisible(NULL))
  }
  fit <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ x + ", rcs_mean, " + age")), data=dd2)
  s <- summary(fit)
  if (!"x" %in% rownames(s$coefficients)) {
    add(model_id=paste0("RECON_",tag,"_30d"), cohort=label, N=nrow(dd2), events=sum(dd2$event_30d),
        note="exposure term dropped (singular/constant); not estimable")
    return(invisible(NULL))
  }
  add(model_id=paste0("RECON_",tag,"_30d"), cohort=label, N=nrow(dd2), events=sum(dd2$event_30d),
      HR=s$coefficients["x","exp(coef)"], lo=s$conf.int["x","lower .95"],
      hi=s$conf.int["x","upper .95"], P=s$coefficients["x","Pr(>|z|)"],
      outcome_version="30d_reconciled_composite_v2", note="")
}

# ---- 1. 分辨率敏感性(修正结局) ----
for (thr in list(list(name="ge2", mask=d$n_gv>=2),
                 list(name="ge3tp", mask=d$n_gv>=3),
                 list(name="ge3tp_ge3val", mask=d$n_gv>=3 & d$n_distinct_values>=3))) {
  dd <- d[thr$mask,]
  for (m in list(list(v="gv_sd", div=10, lbl="SD per 10"),
                 list(v="arv", div=10, lbl="ARV per 10"),
                 list(v="mad_glucose", div=10, lbl="MAD per 10"),
                 list(v="iqr_glucose", div=10, lbl="IQR per 10")))
    fit_rep(dd, m$v, m$div, paste0(thr$name,"_",m$v), paste0(thr$name," (",m$lbl,")"))
}

# ---- 2. 锚点敏感性(修正结局;完全相同样本) ----
for (a in c("opend_0_24h","opstart_0_24h","orout_0_24h","icuin_0_24h","opend_0_12h","opend_0_48h")) {
  f <- ga[ga$anchor==a & ga$n>=2,]
  dd <- d[d$subject_id %in% f$subject_id,]
  dd <- merge(dd, f[,c("subject_id","n","mean_glucose","gv_sd")], by="subject_id", suffixes=c("","_aw"))
  dd$gv_aw <- dd$gv_sd_aw; dd$mean_glu_aw <- dd$mean_glucose_aw
  # 24h landmark 窗锚点(除 opend_0_48h 外);opend_0_48h 使用 48h landmark 版本已在 15 完成,此处仅为对照记录
  dd2 <- dd[complete.cases(dd[,c("t30","event_30d","gv_aw","mean_glu_aw","age")]) & dd$t30>0,]
  if (nrow(dd2) < 100 || sum(dd2$event_30d) < 8) {
    add(model_id=paste0("RECON_ANCH_",a,"_30d"), cohort=paste0("anchor=",a), N=nrow(dd2), events=sum(dd2$event_30d),
        note="not reliably estimable (<100 N or <8 events)")
    next
  }
  dd2$x <- dd2$gv_aw/10
  dd2$mean_glu <- dd2$mean_glu_aw
  fit <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ x + ", rcs_mean, " + age")), data=dd2)
  s <- summary(fit)
  if (!"x" %in% rownames(s$coefficients)) {
    add(model_id=paste0("RECON_ANCH_",a,"_30d"), cohort=paste0("anchor=",a), N=nrow(dd2), events=sum(dd2$event_30d),
        note="exposure term dropped (singular/constant); not estimable")
    next
  }
  add(model_id=paste0("RECON_ANCH_",a,"_30d"), cohort=paste0("anchor=",a), N=nrow(dd2), events=sum(dd2$event_30d),
      HR=s$coefficients["x","exp(coef)"], lo=s$conf.int["x","lower .95"],
      hi=s$conf.int["x","upper .95"], P=s$coefficients["x","Pr(>|z|)"],
      outcome_version="30d_reconciled_composite_v2",
      note=if(a=="opend_0_48h") "24h landmark retained here for comparability; correct 48h version in results/23" else "")
}

# ---- 3. 联合模块 30d 部分(修正结局) ----
j <- d[!is.na(d$hba1c_pct) & !is.na(d$eag_mg_dl) & d$eag_mg_dl > 0,]
j$shr <- j$mean_glu / j$eag_mg_dl
JK <- fromJSON(file.path(ROOT,"results","inspire_joint_constants.json"))
j$shr_zf <- (j$shr - JK$shr_mean)/JK$shr_sd
j$gv_zf <- (j$gv - JK$gv_mean)/JK$gv_sd
j2 <- j[complete.cases(j[,c("t30","event_30d","shr_zf","gv_zf","mean_glu","age","diabetes")]) & j$t30>0,]
if (nrow(j2) >= 100 && sum(j2$event_30d) >= 8) {
  fJ0 <- coxph(Surv(t30,event_30d) ~ age + diabetes, data=j2)
  fJ1 <- coxph(Surv(t30,event_30d) ~ shr_zf + gv_zf + age + diabetes, data=j2)
  V <- vcov(fJ1)[c("shr_zf","gv_zf"),c("shr_zf","gv_zf")]; bb <- coef(fJ1)[c("shr_zf","gv_zf")]
  wstat <- as.numeric(t(bb) %*% solve(V) %*% bb) / 2
  add(model_id="RECON_JOMNI_J1_vs_J0_30d", cohort="joint SHR-GV, reconciled outcome",
      N=nrow(j2), events=sum(j2$event_30d), stat=wstat, df1=2, P=pchisq(wstat*2, df=2, lower.tail=FALSE),
      outcome_version="30d_reconciled_composite_v2", note="2-df Wald (CC)")
  s1 <- summary(fJ1)
  for (trm in c("shr_zf","gv_zf"))
    add(model_id=paste0("RECON_J1_30d_",trm), cohort="joint SHR-GV", N=nrow(j2), events=sum(j2$event_30d),
        term=trm, HR=s1$coefficients[trm,"exp(coef)"], lo=s1$conf.int[trm,"lower .95"],
        hi=s1$conf.int[trm,"upper .95"], P=s1$coefficients[trm,"Pr(>|z|)"], outcome_version="30d_reconciled_composite_v2")
  # D2 vs D1
  kt3 <- paste(format(JK$mean_knots3, digits=10), collapse=",")
  fD1 <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ rcs(mean_glu, c(", kt3, ")) + hba1c_pct + age + diabetes")), data=j2)
  fD2 <- coxph(as.formula(paste0("Surv(t30,event_30d) ~ rcs(mean_glu, c(", kt3, ")) + hba1c_pct + gv_zf + age + diabetes")), data=j2)
  s2 <- summary(fD2)
  add(model_id="RECON_D2_vs_D1_30d", cohort="joint SHR-GV", N=nrow(j2), events=sum(j2$event_30d),
      stat=s2$coefficients["gv_zf","coef"]^2/s2$coefficients["gv_zf","se(coef)"]^2, df1=1,
      P=s2$coefficients["gv_zf","Pr(>|z|)"], gv_HR=s2$coefficients["gv_zf","exp(coef)"],
      outcome_version="30d_reconciled_composite_v2", note="GV added-information beyond RCS(mean)+HbA1c")
} else {
  add(model_id="RECON_JOMNI_J1_vs_J0_30d", cohort="joint SHR-GV, reconciled outcome",
      N=nrow(j2), events=sum(j2$event_30d), note="not reliably estimable (<100 N or <8 events)")
}

tab <- do.call(rbind, lapply(rows, function(x){ allc <- unique(unlist(lapply(rows, names)))
  for (cn in setdiff(allc, names(x))) x[[cn]] <- NA; x[, allc, drop=FALSE] }))
write.csv(tab, file.path(ROOT,"results","20c_corrected_30d_sensitivity_joint.csv"), row.names=FALSE)
print(tab[, intersect(c("model_id","cohort","N","events","HR","lo","hi","P","stat","df1","note"), names(tab))], max=60)
cat("PHASE16_DONE\n")
sink()
