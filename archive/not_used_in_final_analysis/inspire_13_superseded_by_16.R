#!/usr/bin/env Rscript
# 13_inspire_joint_perf.R — INSPIRE SHR–GV 联合模块(预设次要)+ 性能 + QC
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.libPaths(c("~/cardiac_glucose_rebuild_20260728/rlib", .libPaths()))
suppressMessages({library(survival); library(mice); library(jsonlite); library(rms)})
ROOT <- normalizePath("~/inspire_cardiac_20260729")
sink(file.path(ROOT,"logs","13_inspire_joint_perf.log"), split=TRUE)

b <- read.csv(file.path(ROOT,"data","inspire_base.csv"), stringsAsFactors=FALSE)
cm <- read.csv(file.path(ROOT,"data","comorbidity.csv"), stringsAsFactors=FALSE)
hb <- read.csv(file.path(ROOT,"data","hba1c_baseline.csv"), stringsAsFactors=FALSE)
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
d <- merge(b, cm, by="subject_id", all.x=TRUE)
d$gv <- d$gv_sd; d$mean_glu <- d$mean_glucose; d$n_gv <- d$n_glucose_0_24h
d$bmi <- d$weight / (d$height/100)^2
d$diabetes <- as.integer(d$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
d$sex <- factor(d$sex)
d$event_30d <- as.integer(!is.na(d$allcause_death_time) & d$allcause_death_time <= d$opend_time + 30*1440)
d$event_365d <- as.integer(!is.na(d$allcause_death_time) & d$allcause_death_time <= d$opend_time + 365*1440)
d$followup_days <- (pmin(ifelse(is.na(d$allcause_death_time), Inf, d$allcause_death_time), d$discharge_time) - d$opend_time)/1440
d$t30 <- pmin(d$followup_days, 30) - 1; d$t365 <- pmin(d$followup_days, 365) - 1

# ---- SHR 联合队列(inspire_base 自带 hba1c_pct/eag_mg_dl,无需合并) ----
tgt <- d[d$age>=18 & !is.na(d$opend_time) & !is.na(d$gv) & d$n_gv>=2 &
         d$landmark_eligible_allcause %in% c(TRUE,"t","True","TRUE","true",1,"1") & d$t30>=0,]
j <- tgt[!is.na(tgt$hba1c_pct) & !is.na(tgt$eag_mg_dl) & tgt$eag_mg_dl > 0,]
j$shr <- j$mean_glu / j$eag_mg_dl
stopifnot(!anyDuplicated(j$subject_id))
cat("联合队列 N =", nrow(j), "; 30d 事件 =", sum(j$event_30d), "; 365d 事件 =", sum(j$event_365d), "\n")

# ---- 联合队列固定常数(一次冻结) ----
JK <- list(shr_mean=mean(j$shr), shr_sd=sd(j$shr), shr_median=median(j$shr),
  shr_q25=unname(quantile(j$shr,.25)), shr_q75=unname(quantile(j$shr,.75)),
  gv_mean=mean(j$gv), gv_sd=sd(j$gv), gv_median=median(j$gv),
  gv_q25=unname(quantile(j$gv,.25)), gv_q75=unname(quantile(j$gv,.75)),
  mean_knots3=unname(quantile(j$mean_glu, c(.10,.50,.90))),
  N_joint=nrow(j), events_30d=sum(j$event_30d), events_365d=sum(j$event_365d),
  cohort="INSPIRE strict pre-op 1-90d HbA1c, opend landmark, >=2 glucose", seed=SEED)
write_json(JK, file.path(ROOT,"results","inspire_joint_constants.json"), pretty=TRUE, auto_unbox=TRUE)
j$shr_zf <- (j$shr - JK$shr_mean)/JK$shr_sd
j$gv_zf  <- (j$gv  - JK$gv_mean )/JK$gv_sd
j$gv10 <- j$gv/10
kt_mean3 <- paste(format(JK$mean_knots3, digits=10), collapse=",")

# EPV 规则:30d 事件<20 时临床集简化为 age+diabetes(与主分析一致的预设逻辑);365d 用扩展集
CLIN30 <- "age + diabetes"
CLIN365 <- "age + sex + diabetes + charlson_without_diabetes"
clin_of <- function(hz) if(hz==30) CLIN30 else CLIN365

rows_lin <- list(); rows_omni <- list(); rows_dec <- list()
add_lin  <- function(...) rows_lin [[length(rows_lin )+1]] <<- data.frame(..., stringsAsFactors=FALSE)
add_omni <- function(...) rows_omni[[length(rows_omni)+1]] <<- data.frame(..., stringsAsFactors=FALSE)
add_dec  <- function(...) rows_dec [[length(rows_dec )+1]] <<- data.frame(..., stringsAsFactors=FALSE)

cc <- j[complete.cases(j[,c("t30","t365","event_30d","event_365d","shr","gv","mean_glu","hba1c_pct",
  "age","sex","diabetes","charlson_without_diabetes","bmi")]) & j$t30>0,]
cat("联合 CC N =", nrow(cc), "\n")

for (hz in c(30,365)) {
  tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
  covs <- clin_of(hz)
  fJ0 <- as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs))
  fJ1 <- as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs))
  fD1 <- as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + ", covs))
  fD2 <- as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + gv_zf + ", covs))
  J0 <- coxph(fJ0, data=cc); J1 <- coxph(fJ1, data=cc); D1 <- coxph(fD1, data=cc); D2 <- coxph(fD2, data=cc)
  # 联合 Wald(J1 中 shr_zf+gv_zf,2 df,CC)
  V <- vcov(J1)[c("shr_zf","gv_zf"),c("shr_zf","gv_zf")]; bb <- coef(J1)[c("shr_zf","gv_zf")]
  wstat <- as.numeric(t(bb) %*% solve(V) %*% bb) / 2
  add_omni( model_id=paste0("JOMNI_J1_vs_J0_",hz,"d"), analysis=if(hz==30) "lead secondary (joint)" else "key secondary (joint)",
    N=nrow(cc), events=sum(cc[[ev]]), stat=wstat, df1=2, P=pchisq(wstat*2, df=2, lower.tail=FALSE),
    method="2-df Wald (complete-case)", note=if(hz==30) "EPV-limited" else "")
  s1 <- summary(J1)
  for (trm in c("shr_zf","gv_zf"))
    add_lin( model_id=paste0("J1_",hz,"d_",trm), N=nrow(cc), events=sum(cc[[ev]]), term=trm,
      HR=s1$coefficients[trm,"exp(coef)"], lo=s1$conf.int[trm,"lower .95"],
      hi=s1$conf.int[trm,"upper .95"], P=s1$coefficients[trm,"Pr(>|z|)"])
  # D2 vs D1 added-information(Wald gv_zf)
  s2 <- summary(D2)
  add_dec( model_id=paste0("D2_vs_D1_",hz,"d"), N=nrow(cc), events=sum(cc[[ev]]),
    stat=s2$coefficients["gv_zf","coef"]^2 / s2$coefficients["gv_zf","se(coef)"]^2, df1=1,
    P=s2$coefficients["gv_zf","Pr(>|z|)"], gv_HR=s2$coefficients["gv_zf","exp(coef)"],
    gv_lo=s2$conf.int["gv_zf","lower .95"], gv_hi=s2$conf.int["gv_zf","upper .95"],
    note="GV added-information beyond RCS(mean)+HbA1c")
  # 交互(仅 365d 事件充分时)
  if (hz==365) {
    j2 <- cc; j2$shr_c <- j2$shr - JK$shr_median; j2$gv_c <- j2$gv - JK$gv_median
    fI <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_c*gv_c + ", covs)), data=j2)
    sI <- summary(fI)
    add_dec( model_id="INT_SHRxGV_365d", N=nrow(j2), events=sum(j2[[ev]]),
      interaction_HR=exp(sI$coefficients["shr_c:gv_c","coef"]),
      int_P=sI$coefficients["shr_c:gv_c","Pr(>|z|)"], note="exploratory, centered")
  }
}
lin <- do.call(rbind, rows_lin); omni <- do.call(rbind, rows_omni); allc <- unique(unlist(lapply(rows_dec, names)))
dec <- do.call(rbind, lapply(rows_dec, function(x){ for (cn in setdiff(allc, names(x))) x[[cn]] <- NA; x[, allc, drop=FALSE] }))
write.csv(lin, file.path(ROOT,"results","joint_linear_results.csv"), row.names=FALSE)
write.csv(omni, file.path(ROOT,"results","09_joint_omnibus.csv"), row.names=FALSE)
write.csv(dec, file.path(ROOT,"results","10_component_decomposition.csv"), row.names=FALSE)
write.csv(lin, file.path(ROOT,"results","08_hba1c_shr_results.csv"), row.names=FALSE)

# ---- 角点绝对风险(共同支持检查后) ----
ch <- chull(cc$shr, cc$gv); hp <- cc[ch, c("shr","gv")]
in_hull <- function(x, y){
  px <- hp$shr; py <- hp$gv; n <- length(px); out <- logical(length(x))
  for (k in seq_along(x)) {
    cr <- 0
    for (i in 1:n) { j2 <- if(i==n) 1 else i+1
      if (((py[i] > y[k]) != (py[j2] > y[k])) && (x[k] < (px[j2]-px[i])*(y[k]-py[i])/(py[j2]-py[i]) + px[i])) cr <- cr+1 }
    out[k] <- (cr %% 2)==1 }
  out
}
corners <- data.frame(scenario=c("SHR P25/GV P25","SHR P75/GV P25","SHR P25/GV P75","SHR P75/GV P75"),
  shr=c(JK$shr_q25,JK$shr_q75,JK$shr_q25,JK$shr_q75), gv=c(JK$gv_q25,JK$gv_q25,JK$gv_q75,JK$gv_q75))
corners$inside <- in_hull(corners$shr, corners$gv)
write.csv(corners, file.path(ROOT,"results","joint_corner_support.csv"), row.names=FALSE)
risk_rows <- list()
for (hz in c(30,365)) {
  tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
  tmax <- if(hz==30) 29 else 364
  covs <- clin_of(hz)
  fit_risk <- function(dd){
    f <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs)), data=dd)
    bh <- basehaz(f, centered=FALSE); H <- max(bh$hazard[bh$time<=tmax])
    sapply(1:4, function(k){
      nd <- dd; nd$shr <- corners$shr[k]; nd$gv <- corners$gv[k]
      nd$shr_zf <- (nd$shr - JK$shr_mean)/JK$shr_sd; nd$gv_zf <- (nd$gv - JK$gv_mean)/JK$gv_sd
      lp <- predict(f, newdata=nd, type="lp", reference="zero")
      mean(1 - exp(-H*exp(lp)))
    })
  }
  if (hz==365 && sum(cc$discharge_time - cc$opend_time >= 30*1440) < 50) {
    risk_rows[[length(risk_rows)+1]] <- data.frame(horizon=365, note="365d absolute risk not reliably estimable (<50 patients observed >=30d in-hospital)")
    next
  }
  pt <- fit_risk(cc)
  B <- 1000; vals <- matrix(NA, B, 4); fails <- 0L; n <- nrow(cc)
  for (b_ in 1:B) {
    idx <- sample.int(n, n, replace=TRUE)
    rb <- tryCatch(fit_risk(cc[idx,]), error=function(e) NULL)
    if (is.null(rb) || any(!is.finite(rb))) { fails <- fails+1L; next }
    vals[b_,] <- rb
  }
  ci <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
  for (k in 1:4)
    risk_rows[[length(risk_rows)+1]] <- data.frame(model_id=paste0("JRISK_",hz,"d_c",k), horizon=hz,
      scenario=corners$scenario[k], inside_support=corners$inside[k],
      risk=pt[k], risk_lo=ci(vals[,k])[1], risk_hi=ci(vals[,k])[2],
      rd_vs_P25P25=pt[k]-pt[1], rd_lo=ci(vals[,k]-vals[,1])[1], rd_hi=ci(vals[,k]-vals[,1])[2],
      n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
}
allc_r <- unique(unlist(lapply(risk_rows, names)))
risk_tab <- do.call(rbind, lapply(risk_rows, function(x){ for (cn in setdiff(allc_r, names(x))) x[[cn]] <- NA; x[, allc_r, drop=FALSE] }))
write.csv(risk_tab, file.path(ROOT,"results","joint_absolute_risk.csv"), row.names=FALSE)
print(risk_tab)

# ---- ratio vs component 成对性能(1000 bootstrap) ----
cstat <- function(fit, tm, ev, dd){
  lp <- predict(fit, newdata=dd, type="lp")
  1 - survival::concordance(survival::Surv(dd[[tm]], dd[[ev]]) ~ lp)$concordance
}
pairs_rows <- list()
for (hz in c(30,365)) {
  tm <- paste0("t",hz); ev <- paste0("event_",hz,"d")
  covs <- clin_of(hz)
  fmls <- list(
    P0=as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs)),
    P3=as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs)),
    P5=as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + gv_zf + ", covs)))
  fits0 <- lapply(fmls, function(f) coxph(f, data=cc))
  ap <- sapply(fits0, function(f) cstat(f,tm,ev,cc))
  B <- 1000; n <- nrow(cc); opt <- matrix(0,B,3); dcs <- matrix(NA,B,2); fails <- 0L
  for (b_ in 1:B) {
    idx <- sample.int(n, n, replace=TRUE); db <- cc[idx,]
    r <- tryCatch({
      fb <- lapply(fmls, function(f) coxph(f, data=db)); names(fb) <- names(fmls)
      c(sapply(fb, function(f) cstat(f,tm,ev,db)) - sapply(fb, function(f) cstat(f,tm,ev,cc)),
        cstat(fb$P3,tm,ev,cc)-cstat(fb$P0,tm,ev,cc), cstat(fb$P3,tm,ev,cc)-cstat(fb$P5,tm,ev,cc))
    }, error=function(e) NULL)
    if (is.null(r) || any(!is.finite(r))) { fails <- fails+1L; next }
    opt[b_,] <- r[1:3]; dcs[b_,] <- r[4:5]
  }
  ci <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
  for (i in 1:3)
    pairs_rows[[length(pairs_rows)+1]] <- data.frame(model_id=paste0("PERF_",names(fmls)[i],"_",hz,"d"),
      model=c("clinical (P0)","ratio SHR+GV (P3)","component mean+HbA1c+GV (P5)")[i],
      horizon=hz, c_apparent=ap[i], optimism=mean(opt[,i]), c_corrected=ap[i]-mean(opt[,i]),
      AIC=AIC(fits0[[i]]), n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  pairs_rows[[length(pairs_rows)+1]] <- data.frame(model_id=paste0("PERF_pairs_",hz,"d"), horizon=hz,
    dC_P3_vs_P0=mean(dcs[,1]), dC_P3P0_lo=ci(dcs[,1])[1], dC_P3P0_hi=ci(dcs[,1])[2],
    dC_P3_vs_P5=mean(dcs[,2]), dC_P3P5_lo=ci(dcs[,2])[1], dC_P3P5_hi=ci(dcs[,2])[2],
    note="paired bootstrap; ratio vs component", n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
}
allc_p <- unique(unlist(lapply(pairs_rows, names)))
perf <- do.call(rbind, lapply(pairs_rows, function(x){ for (cn in setdiff(allc_p, names(x))) x[[cn]] <- NA; x[, allc_p, drop=FALSE] }))
write.csv(perf, file.path(ROOT,"results","11_performance.csv"), row.names=FALSE)
print(perf)
cat("PHASE13_DONE\n")
sink()
