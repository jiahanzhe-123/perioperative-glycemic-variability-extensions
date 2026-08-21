# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 04_ph_rcs_absrisk.R — PH 诊断与处理、GV 非线性、Royston–Parmar 对照、绝对风险 bootstrap。
# 输入:data/analysis_base.csv + data/features_priority.csv + results/standardization_constants.json
# 输出:results/07_ph_diagnostics.csv, results/rcs_gv_nonlinearity.csv, results/rcs_gv_curve_source.csv,
#       results/12_absolute_risk_results.csv, results/bootstrap_absrisk_iterations.csv,
#       results/flexsurv_crosscheck.csv, results/collinearity_support.csv,
#       results/common_support_by_mean_decile.csv, figures/rcs_gv_*.png|pdf
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite); library(rms); library(ggplot2); library(flexsurv)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/04_ph_rcs_absrisk.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$span_hours <- d$ps_span_hours; d$frac_central_lab <- d$ps_frac_central_lab
d$frac_blood_gas <- d$ps_frac_blood_gas; d$frac_poct <- d$ps_frac_poct; d$frac_icu_charted <- d$ps_frac_icu_charted
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
tgt$gv10 <- tgt$gv/10; tgt$log_count <- log(tgt$glucose_count)
tgt$gender <- factor(tgt$gender)
tgt$procedure_cat6 <- factor(tgt$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
cc <- tgt[complete.cases(tgt[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","gv","gv10","mean_glu",COVS)]) & tgt$t_lm_30>0,]
cat("analysis(complete-case) N =", nrow(cc), "\n")

# ---- 1. PH 诊断与处理 ----
ph_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  fB <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml)), data=cc, x=TRUE, y=TRUE)
  z <- cox.zph(fB, transform="km")
  tab <- as.data.frame(z$table); tab$term <- rownames(tab)
  for (i in seq_len(nrow(tab)))
    ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("ModelB_",hz,"d"), term=tab$term[i],
      chisq=tab$chisq[i], df=tab$df[i], p=tab$p[i],
      interval=NA, interval_events=NA, HR_per10=NA, note="", stringsAsFactors=FALSE)
  pv <- tab$p; names(pv) <- tab$term
  cat(sprintf("[%sd] PH global=%.4g, gv10=%.4g\n", hz, pv["GLOBAL"], pv["gv10"]))
  # GV 违反或常规报告:预设区间平均效应(1–7, 8–30;365d 加 31–365;区间事件<10 不单独报告)
  cuts <- if(hz=="30") c(0, 7, 30) else c(0, 7, 30, 365)
  dd_split <- cc[, c("stay_id", tm, ev, "gv","gv10","mean_glu", COVS)]
  sp <- survSplit(as.formula(paste0("Surv(",tm,",",ev,") ~ .")), data=dd_split, cut=cuts[-length(cuts)],
                  episode="interval", start="tstart", end="tstop", event="evsplit")
  sp$interval <- factor(sp$interval, labels=paste0("d", cuts[-length(cuts)]+1, "-", cuts[-1]))
  sp <- sp[!is.na(sp$evsplit),]
  evn <- tapply(sp$evsplit, sp$interval, sum)
  f_int <- coxph(as.formula(paste0("Surv(tstart,tstop,evsplit) ~ gv10:interval + ", rcs_mean, " + ", covs_fml)), data=sp)
  si <- summary(f_int)$coefficients
  b0 <- coef(f_int); V <- vcov(f_int)
  for (lv in levels(sp$interval)) {
    term <- paste0("gv10:interval", lv)
    insufficient <- is.na(evn[lv]) || evn[lv] < 10
    if (!term %in% rownames(si)) next
    # 区间主效应 = gv10:interval 项本身(参照区间为第一段)
    ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("ModelB_",hz,"d_GV_intervals"), term=term,
      chisq=NA, df=NA, p=if(insufficient) NA else si[term,"Pr(>|z|)"],
      interval=lv, interval_events=as.integer(evn[lv]),
      HR_per10=if(insufficient) NA else exp(si[term,"coef"]),
      note=if(insufficient) "insufficient events" else "", stringsAsFactors=FALSE)
  }
  # 连续协变量违反 → log(time) 交互(记录处理动作)
  viol_cont <- intersect(names(pv)[pv < 0.05],
    c("age_at_admission","bmi","charlson_without_diabetes","lactate_postop_first","creat_postop_first","sofa_24h"))
  if (length(viol_cont)) {
    tt_fml <- paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml, " + ",
                     paste(paste0("tt(", viol_cont, ")"), collapse=" + "))
    f_tt <- tryCatch(coxph(as.formula(tt_fml), data=cc, tt=function(x, t, ...) x*log(t)), error=function(e) NULL)
    ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("ModelB_",hz,"d_tt_corrected"),
      term=paste(viol_cont, collapse="+"), chisq=NA, df=NA, p=NA, interval=NA, interval_events=NA_integer_,
      HR_per10=if(!is.null(f_tt)) exp(summary(f_tt)$coefficients["gv10","coef"]) else NA,
      note=paste0("log(time) interaction added for: ", paste(viol_cont, collapse=","),
                  if(is.null(f_tt)) " (fit failed)" else ""), stringsAsFactors=FALSE)
  }
}
ph_tab <- do.call(rbind, ph_rows)
write.csv(ph_tab, file.path(ROOT,"results","07_ph_diagnostics.csv"), row.names=FALSE)

# ---- 2. GV 非线性(RCS;30d 3 knots,365d 4 knots;曲线 2.5–97.5 百分位,参考=中位数) ----
lrt_p <- function(a,b){ ddf <- attr(logLik(b),"df")-attr(logLik(a),"df")
  if (!is.finite(ddf)||ddf<=0) return(NA_real_); pchisq(2*as.numeric(logLik(b)-logLik(a)), df=ddf, lower.tail=FALSE) }
nl_rows <- list(); curve_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  knots <- if(hz=="30") K$gv_knots_30d else K$gv_knots_365d
  kt <- paste(format(knots, digits=10), collapse=",")
  fS <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(gv, c(", kt, ")) + ", rcs_mean, " + ", covs_fml)), data=cc, x=TRUE, y=TRUE)
  fL <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv + ", rcs_mean, " + ", covs_fml)), data=cc, x=TRUE, y=TRUE)
  f0 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ ", rcs_mean, " + ", covs_fml)), data=cc, x=TRUE, y=TRUE)
  op <- lrt_p(f0,fS); np <- lrt_p(fL,fS)
  nl_rows[[length(nl_rows)+1]] <- data.frame(model_id=paste0("RCS_GV_",hz,"d"), horizon=hz,
    N=nrow(cc), events=sum(cc[[ev]]), knots=paste(round(knots,3),collapse=";"),
    overall_p=op, nonlinear_p=np, stringsAsFactors=FALSE)
  sc <- which(attr(fS$x,"assign")==1); b <- coef(fS)[sc]; V <- vcov(fS)[sc,sc,drop=FALSE]
  grid <- seq(quantile(cc$gv,.025), quantile(cc$gv,.975), length.out=200)
  Bg <- rcs(grid, knots); Br <- rcs(median(cc$gv), knots); D <- sweep(Bg,2,Br,"-")
  lp <- as.vector(D%*%b); se <- sqrt(pmax(0,rowSums((D%*%V)*D)))
  curve_rows[[length(curve_rows)+1]] <- data.frame(horizon=hz, gv=grid, HR=exp(lp),
    lo=exp(lp-1.96*se), hi=exp(lp+1.96*se), reference=median(cc$gv),
    overall_p=op, nonlinear_p=np, N=nrow(cc), events=sum(cc[[ev]]))
}
nl_tab <- do.call(rbind, nl_rows); curves <- do.call(rbind, curve_rows)
write.csv(nl_tab, file.path(ROOT,"results","rcs_gv_nonlinearity.csv"), row.names=FALSE)
write.csv(curves, file.path(ROOT,"results","rcs_gv_curve_source.csv"), row.names=FALSE)
rug <- data.frame(gv=cc$gv)
for (hz in c("30","365")) {
  cv <- curves[curves$horizon==hz,]
  p <- ggplot(cv, aes(x=gv)) +
    geom_rug(data=rug, aes(x=gv), inherit.aes=FALSE, alpha=0.03, color="#0E7C7B") +
    geom_ribbon(aes(ymin=lo,ymax=hi), fill="#0E7C7B", alpha=0.18) +
    geom_line(aes(y=HR), color="#0E7C7B", linewidth=0.9) +
    geom_hline(yintercept=1, linetype="dashed", color="grey55") +
    geom_vline(xintercept=median(cc$gv), linetype="dotted", color="grey40") +
    annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=1.5, size=3.2, color="grey30",
      label=sprintf("overall P=%.3f\nnonlinear P=%.3f\nN=%s, events=%s",
        unique(cv$overall_p), unique(cv$nonlinear_p), format(unique(cv$N),big.mark=","), unique(cv$events))) +
    labs(title=sprintf("GV spline, mortality by index day %s (day-1 landmark survivors)", hz),
         x="GV (SD of date-anchored 24-h glucose, mg/dL)", y="Adjusted HR (95% CI)") +
    theme_bw(base_size=10.5)
  ggsave(file.path(ROOT,"figures",paste0("rcs_gv_",hz,"d.png")), p, width=7, height=5, dpi=600)
  ggsave(file.path(ROOT,"figures",paste0("rcs_gv_",hz,"d.pdf")), p, width=7, height=5)
}

# ---- 3. 共线性与支持域 ----
cc$gv_zf <- (cc$gv - K$gv_mean)/K$gv_sd; cc$mean_zf <- (cc$mean_glu - K$mean_glu_mean)/K$mean_glu_sd
pear <- cor(cc$gv, cc$mean_glu); spear <- cor(cc$gv, cc$mean_glu, method="spearman")
numv <- c("age_at_admission","bmi","charlson_without_diabetes","sofa_24h","lactate_postop_first","creat_postop_first")
l1 <- lm(as.formula(paste("gv_zf ~ mean_zf +", paste(numv, collapse="+"))), data=cc)
vif_gv <- 1/(1-summary(l1)$r.squared)
l2 <- lm(as.formula(paste("mean_zf ~ gv_zf +", paste(numv, collapse="+"))), data=cc)
vif_mean <- 1/(1-summary(l2)$r.squared)
X <- scale(cc[,c("gv","mean_glu",numv)])
sv <- svd(X)$d; cidx <- max(sv)/sv
col_tab <- data.frame(item=c("pearson_gv_mean","spearman_gv_mean","vif_gv","vif_mean","condition_index_max","condition_index_min"),
                      value=c(pear,spear,vif_gv,vif_mean,max(cidx),min(cidx)))
write.csv(col_tab, file.path(ROOT,"results","collinearity_support.csv"), row.names=FALSE)
cc$mean_dec <- cut(cc$mean_glu, breaks=quantile(cc$mean_glu, seq(0,1,.1)), include.lowest=TRUE)
sup <- aggregate(gv ~ mean_dec, data=cc,
  FUN=function(x) paste0("n=",length(x),"; q25=",round(quantile(x,.25),1),"; q75=",round(quantile(x,.75),1)))
write.csv(sup, file.path(ROOT,"results","common_support_by_mean_decile.csv"), row.names=FALSE)

# ---- 4. Royston–Parmar 对照 ----
fs_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  fs <- tryCatch(flexsurvspline(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_fml)),
                 data=cc, k=3, scale="hazard"), error=function(e) e)
  if (inherits(fs,"flexsurvreg")) {
    cf <- fs$coefficients["gv10"]; se <- sqrt(diag(vcov(fs)))["gv10"]
    fs_rows[[length(fs_rows)+1]] <- data.frame(model_id=paste0("RP_",hz,"d"), horizon=hz,
      HR_per10=exp(cf), lo=exp(cf-1.96*se), hi=exp(cf+1.96*se), k=3,
      note="Royston-Parmar flexible baseline; covariates as Model A; cross-check vs Cox", stringsAsFactors=FALSE)
  } else fs_rows[[length(fs_rows)+1]] <- data.frame(model_id=paste0("RP_",hz,"d"), horizon=hz,
      HR_per10=NA, lo=NA, hi=NA, k=3, note=paste("flexsurvspline failed:", conditionMessage(fs)), stringsAsFactors=FALSE)
}
write.csv(do.call(rbind,fs_rows), file.path(ROOT,"results","flexsurv_crosscheck.csv"), row.names=FALSE)

# ---- 5. 绝对风险(Q75 vs Q25 标准化;1000 次患者级 bootstrap,每次完整重拟合) ----
risk_contrast <- function(dd, hz){
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  tmax <- if(hz=="30") 29 else 364
  fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv + ", rcs_mean, " + ", covs_fml)), data=dd)
  bh <- basehaz(fit, centered=FALSE)
  bh <- bh[bh$time <= tmax,]
  tt <- c(0, bh$time); H <- c(0, bh$hazard); dt <- diff(c(tt, tmax))
  lp25 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q25; nd}, type="lp", reference="zero")
  lp75 <- predict(fit, newdata={nd <- dd; nd$gv <- K$gv_q75; nd}, type="lp", reference="zero")
  risk_of <- function(lpv, Ht) 1 - exp(-Ht * exp(lpv))
  r25 <- risk_of(lp25, max(H)); r75 <- risk_of(lp75, max(H))
  rmst_of <- function(lpv){
    S <- exp(-outer(H, exp(lpv)))
    colSums(S[-nrow(S), , drop=FALSE] * dt)
  }
  m25 <- mean(r25); m75 <- mean(r75)
  list(risk25=m25, risk75=m75, rd=m75-m25, rr=m75/m25,
       rmst_diff=if(hz=="365") mean(rmst_of(lp75)) - mean(rmst_of(lp25)) else NA_real_)
}
boot_out <- list(); iter_rows <- list()
for (hz in c("30","365")) {
  ev <- paste0("event_lm_",hz)
  pt <- risk_contrast(cc, hz)
  B <- 1000; vals <- matrix(NA, nrow=B, ncol=5); fails <- 0L; n <- nrow(cc)
  t0 <- Sys.time()
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE)
    dd <- cc[idx,]
    rb <- tryCatch({ r <- risk_contrast(dd, hz); c(r$rd, r$rr, r$risk25, r$risk75, ifelse(is.na(r$rmst_diff), NA, r$rmst_diff)) },
                   error=function(e) NULL)
    if (is.null(rb) || any(!is.finite(rb[c(1,2,3,4)]))) { fails <- fails+1L; next }
    vals[b,] <- rb
  }
  cat(sprintf("[absrisk %sd] %d 次,%.1f 分钟,失败 %d\n", hz, B, as.numeric(difftime(Sys.time(),t0,units="mins")), fails))
  colnames(vals) <- c("rd","rr","risk25","risk75","rmst_diff")
  iter_rows[[length(iter_rows)+1]] <- data.frame(horizon=hz, iteration=1:B, vals)
  ci <- apply(vals, 2, quantile, probs=c(.025,.975), na.rm=TRUE)
  boot_out[[length(boot_out)+1]] <- data.frame(
    model_id=paste0("ABSRISK_",hz,"d"), horizon=hz, N=nrow(cc), events=sum(cc[[ev]]),
    gv_q25=K$gv_q25, gv_q75=K$gv_q75,
    risk_q25=pt$risk25, risk_q75=pt$risk75,
    rd=pt$rd, rd_lo=ci[1,"rd"], rd_hi=ci[2,"rd"],
    rr=pt$rr, rr_lo=ci[1,"rr"], rr_hi=ci[2,"rr"],
    rmst_diff=pt$rmst_diff, rmst_lo=if(hz=="365") ci[1,"rmst_diff"] else NA, rmst_hi=if(hz=="365") ci[2,"rmst_diff"] else NA,
    n_boot=B, n_boot_failed=fails, method="percentile bootstrap, patient-level, full model refit",
    stringsAsFactors=FALSE)
}
write.csv(do.call(rbind,boot_out), file.path(ROOT,"results","12_absolute_risk_results.csv"), row.names=FALSE)
write.csv(do.call(rbind,iter_rows), file.path(ROOT,"results","bootstrap_absrisk_iterations.csv"), row.names=FALSE)
cat("PHASE4_DONE\n")
sink()
