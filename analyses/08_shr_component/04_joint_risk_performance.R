# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 22_joint_risk_perf.R — 联合绝对风险(4 角点 + 二维风险面 + 1000 bootstrap)、
# 模型性能成对比较(P0–P5)、bootstrap 稳定性、影响点诊断。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite); library(rms); library(ggplot2)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/22_joint_risk_perf.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$eAG <- 28.7*d$hba1c_pct - 46.7; d$shr <- d$mean_glu/d$eAG
d$gender <- factor(d$gender)
d$procedure_cat6 <- factor(d$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
JK <- fromJSON(file.path(ROOT,"results","shr_gv_joint_constants.json"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
j  <- tgt[tgt$hba1c_tier=="1-90 days pre-index" & tgt$eAG>0,]
cc <- j[complete.cases(j[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","shr","gv","mean_glu",COVS)]) & j$t_lm_30>0,]
cat("absolute-risk CC N =", nrow(cc), "\n")
prep <- function(dd){
  dd$shr_zf <- (dd$shr - JK$shr_mean)/JK$shr_sd
  dd$gv_zf  <- (dd$gv  - JK$gv_mean )/JK$gv_sd
  dd
}
cc <- prep(cc)
kt_mean3 <- paste(format(JK$mean_knots3, digits=10), collapse=",")

# ---- 1. 绝对风险:四角情景 + bootstrap ----
corners <- data.frame(
  scenario=c("SHR P25 / GV P25","SHR P75 / GV P25","SHR P25 / GV P75","SHR P75 / GV P75"),
  shr=c(JK$shr_q25, JK$shr_q75, JK$shr_q25, JK$shr_q75),
  gv=c(JK$gv_q25, JK$gv_q25, JK$gv_q75, JK$gv_q75))
fit_risk <- function(dd, hz){
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  tmax <- if(hz=="30") 29 else 364
  fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml)), data=dd)
  bh <- basehaz(fit, centered=FALSE); bh <- bh[bh$time<=tmax,]; H <- max(bh$hazard)
  risks <- sapply(1:nrow(corners), function(k){
    nd <- dd; nd$shr <- corners$shr[k]; nd$gv <- corners$gv[k]
    nd$shr_zf <- (nd$shr - JK$shr_mean)/JK$shr_sd; nd$gv_zf <- (nd$gv - JK$gv_mean)/JK$gv_sd
    lp <- predict(fit, newdata=nd, type="lp", reference="zero")
    mean(1 - exp(-H*exp(lp)))
  })
  rmst365 <- NA
  if (hz=="365") {
    tt <- c(0, bh$time); Hv <- c(0, bh$hazard); dt <- diff(c(tt, tmax))
    rmst365 <- sapply(1:nrow(corners), function(k){
      nd <- dd; nd$shr <- corners$shr[k]; nd$gv <- corners$gv[k]
      nd$shr_zf <- (nd$shr - JK$shr_mean)/JK$shr_sd; nd$gv_zf <- (nd$gv - JK$gv_mean)/JK$gv_sd
      lp <- predict(fit, newdata=nd, type="lp", reference="zero")
      S <- exp(-outer(Hv, exp(lp)))
      mean(colSums(S[-nrow(S),,drop=FALSE]*dt))
    })
  }
  list(risks=risks, rmst=rmst365)
}
risk_rows <- list(); iter_rows <- list()
for (hz in c("30","365")) {
  ev <- paste0("event_lm_",hz)
  pt <- fit_risk(cc, hz)
  B <- 1000; vals <- matrix(NA, nrow=B, ncol=8); fails <- 0L; n <- nrow(cc)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE); db <- cc[idx,]
    r <- tryCatch(fit_risk(db, hz), error=function(e) NULL)
    if (is.null(r)) { fails <- fails+1L; next }
    vals[b,] <- c(r$risks, if(hz=="365") r$rmst else rep(NA,4))
  }
  ci <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
  for (k in 1:4) {
    rd <- pt$risks[k] - pt$risks[1]; rr <- pt$risks[k]/pt$risks[1]
    rd_b <- vals[,k] - vals[,1]; rr_b <- vals[,k]/vals[,1]
    risk_rows[[length(risk_rows)+1]] <- data.frame(
      model_id=paste0("ABSRISK_joint_",hz,"d_corner",k), horizon=hz, scenario=corners$scenario[k],
      N=nrow(cc), events=sum(cc[[ev]]),
      risk=pt$risks[k], risk_lo=ci(vals[,k])[1], risk_hi=ci(vals[,k])[2],
      rd_vs_P25P25=rd, rd_lo=ci(rd_b)[1], rd_hi=ci(rd_b)[2],
      rr_vs_P25P25=rr, rr_lo=ci(rr_b)[1], rr_hi=ci(rr_b)[2],
      rmst=if(hz=="365") pt$rmst[k] else NA,
      rmst_diff_vs_P25P25=if(hz=="365") pt$rmst[k]-pt$rmst[1] else NA,
      n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  }
  it <- as.data.frame(vals); names(it) <- paste0("corner",1:4, rep(c("_risk","_rmst")[1:2], each=4))
  it$horizon <- hz; iter_rows[[length(iter_rows)+1]] <- it
  cat(sprintf("[risk %sd] 失败 %d\n", hz, fails))
}
risk_tab <- do.call(rbind, risk_rows)
write.csv(risk_tab, file.path(ROOT,"results","shr_gv_absolute_risk_surface.csv"), row.names=FALSE)
write.csv(do.call(rbind, iter_rows), file.path(ROOT,"results","shr_gv_abrisk_bootstrap_iterations.csv"), row.names=FALSE)
print(risk_tab[,c("horizon","scenario","risk","rd_vs_P25P25","rd_lo","rd_hi","rr_vs_P25P25","rr_lo","rr_hi")])

# ---- 2. 二维标准化风险面(30d) ----
fit30 <- coxph(as.formula(paste0("Surv(t_lm_30, event_lm_30) ~ shr_zf + gv_zf + ", covs_fml)), data=cc)
bh <- basehaz(fit30, centered=FALSE); H29 <- max(bh$hazard[bh$time<=29])
gx <- seq(quantile(cc$shr,.05), quantile(cc$shr,.95), length.out=60)
gy <- seq(quantile(cc$gv,.05),  quantile(cc$gv,.95),  length.out=60)
grid <- expand.grid(shr=gx, gv=gy)
# 共同支持遮罩(凸包)
ch <- chull(cc$shr, cc$gv); hull_pts <- cc[ch, c("shr","gv")]
in_hull <- function(x, y){
  # 射线法
  px <- hull_pts$shr; py <- hull_pts$gv; n <- length(px)
  inside <- logical(length(x))
  for (k in seq_along(x)) {
    cr <- 0
    for (i in 1:n) {
      j2 <- if(i==n) 1 else i+1
      if (((py[i] > y[k]) != (py[j2] > y[k])) &&
          (x[k] < (px[j2]-px[i])*(y[k]-py[i])/(py[j2]-py[i]) + px[i])) cr <- cr+1
    }
    inside[k] <- (cr %% 2)==1
  }
  inside
}
grid$inside <- in_hull(grid$shr, grid$gv)
risk_at <- function(s_, g_){
  nd <- cc; nd$shr <- s_; nd$gv <- g_
  nd$shr_zf <- (s_ - JK$shr_mean)/JK$shr_sd; nd$gv_zf <- (g_ - JK$gv_mean)/JK$gv_sd
  lp <- predict(fit30, newdata=nd, type="lp", reference="zero")
  mean(1 - exp(-H29*exp(lp)))
}
grid$risk <- mapply(risk_at, grid$shr, grid$gv)
grid$risk_masked <- ifelse(grid$inside, grid$risk, NA)
write.csv(grid, file.path(ROOT,"results","shr_gv_risk_surface_grid.csv"), row.names=FALSE)
p <- ggplot(grid, aes(shr, gv)) +
  geom_raster(aes(fill=risk_masked)) +
  geom_point(data=cc, alpha=0.04, size=0.5, color="grey30") +
  geom_path(data=hull_pts, aes(shr, gv), color="red", linewidth=0.6) +
  scale_fill_viridis_c(name="standardized\n30-day risk", na.value="grey90") +
  geom_point(data=corners, aes(shr, gv), shape=23, fill="white", size=3) +
  labs(title="SHR–GV standardized absolute-risk surface (masked outside common support)",
       subtitle="day-1 landmark cohort, strictly pre-index 1-90d HbA1c; risk by index day 30",
       x="SHR", y="GV (SD, mg/dL)") + theme_bw()
ggsave(file.path(ROOT,"figures","joint_risk_surface.png"), p, width=7.5, height=6, dpi=600)
ggsave(file.path(ROOT,"figures","joint_risk_surface.pdf"), p, width=7.5, height=6)
# 四角情景图
rt30 <- risk_tab[risk_tab$horizon=="30",]
p2 <- ggplot(rt30, aes(x=scenario, y=risk)) +
  geom_point(size=2.5, color="#0072B2") +
  geom_errorbar(aes(ymin=risk_lo, ymax=risk_hi), width=0.25, color="#0072B2") +
  geom_point(aes(y=risk[1]), color="grey50", shape=4) +
  coord_flip() +
  labs(title="Prespecified joint risk scenarios (30-day, standardized)",
       y="standardized risk by index day 30", x=NULL) + theme_bw()
ggsave(file.path(ROOT,"figures","joint_risk_corners.png"), p2, width=7, height=4, dpi=600)
ggsave(file.path(ROOT,"figures","joint_risk_corners.pdf"), p2, width=7, height=4)

# ---- 3. 模型性能(P0–P5,成对 bootstrap) ----
mk_fml <- function(kind, tm, ev){
  switch(kind,
    P0=as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_fml)),
    P1=as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + ", covs_fml)),
    P2=as.formula(paste0("Surv(",tm,",",ev,") ~ gv_zf + ", covs_fml)),
    P3=as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml)),
    P4=as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + ", covs_fml)),
    P5=as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + gv_zf + ", covs_fml)))
}
cstat <- function(fit, tm, ev, dd){
  lp <- predict(fit, newdata=dd, type="lp")
  1 - survival::concordance(survival::Surv(dd[[tm]], dd[[ev]]) ~ lp)$concordance
}
brier_ipc <- function(fit, tm, ev, dd, tmax){
  bh <- basehaz(fit, centered=FALSE); H <- max(bh$hazard[bh$time<=tmax])
  lp <- predict(fit, newdata=dd, type="lp", reference="zero")
  p <- 1 - exp(-H*exp(lp))
  y <- dd[[ev]]; y[dd[[tm]]<tmax & dd[[ev]]==0] <- NA
  mean((y-p)^2, na.rm=TRUE)
}
slope_of <- function(fit, tm, ev, dd){
  lp <- predict(fit, newdata=dd, type="lp")
  coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ lp")), data=dd)$coefficients["lp"]
}
perf_rows <- list(); pairs_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz); tmax <- if(hz=="30") 29 else 364
  kinds <- paste0("P",0:5)
  fits0 <- lapply(kinds, function(k) coxph(mk_fml(k,tm,ev), data=cc))
  names(fits0) <- kinds
  ap <- sapply(fits0, function(f) cstat(f,tm,ev,cc))
  br <- sapply(fits0, function(f) brier_ipc(f,tm,ev,cc,tmax))
  sl <- sapply(fits0, function(f) slope_of(f,tm,ev,cc))
  aics <- sapply(fits0, AIC)
  B <- 1000; n <- nrow(cc)
  opt <- matrix(0, B, 6); dcs <- matrix(NA, B, 3)
  set.seed(SEED + as.integer(hz))
  fails <- 0L
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE); db <- cc[idx,]
    r <- tryCatch({
      fb <- lapply(kinds, function(k) coxph(mk_fml(k,tm,ev), data=db))
      names(fb) <- kinds
      c(sapply(fb, function(f) cstat(f,tm,ev,db)) - sapply(fb, function(f) cstat(f,tm,ev,cc)),
        cstat(fb$P3,tm,ev,cc)-cstat(fb$P0,tm,ev,cc),
        cstat(fb$P5,tm,ev,cc)-cstat(fb$P4,tm,ev,cc),
        cstat(fb$P3,tm,ev,cc)-cstat(fb$P5,tm,ev,cc))
    }, error=function(e) NULL)
    if (is.null(r) || any(!is.finite(r))) { fails <- fails+1L; next }
    opt[b,] <- r[1:6]; dcs[b,] <- r[7:9]
  }
  ci <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
  for (i in 1:6) {
    perf_rows[[length(perf_rows)+1]] <- data.frame(model_id=paste0("PERF_P",i-1,"_",hz,"d"), horizon=hz,
      model=c("clinical","clinical+SHR","clinical+GV","clinical+SHR+GV (J1)",
              "clinical+flex mean+HbA1c","clinical+flex mean+HbA1c+GV (D2)")[i],
      N=nrow(cc), events=sum(cc[[ev]]),
      c_apparent=ap[i], optimism=mean(opt[,i]), c_corrected=ap[i]-mean(opt[,i]),
      brier=br[i], calibration_slope_apparent=sl[i], AIC=aics[i],
      n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  }
  pairs_rows[[length(pairs_rows)+1]] <- data.frame(model_id=paste0("PERF_pairs_",hz,"d"), horizon=hz,
    dC_P3_vs_P0=mean(dcs[,1]), dC_P3P0_lo=ci(dcs[,1])[1], dC_P3P0_hi=ci(dcs[,1])[2],
    dC_P5_vs_P4=mean(dcs[,2]), dC_P5P4_lo=ci(dcs[,2])[1], dC_P5P4_hi=ci(dcs[,2])[2],
    dC_P3_vs_P5=mean(dcs[,3]), dC_P3P5_lo=ci(dcs[,3])[1], dC_P3P5_hi=ci(dcs[,3])[2],
    note="paired bootstrap, identical resamples; ratio formulation=P3, component formulation=P5",
    n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  cat(sprintf("[perf %sd] 失败 %d\n", hz, fails))
}
perf_tab <- do.call(rbind, perf_rows); pairs_tab <- do.call(rbind, pairs_rows)
write.csv(perf_tab, file.path(ROOT,"results","shr_gv_model_performance.csv"), row.names=FALSE)
write.csv(pairs_tab, file.path(ROOT,"results","shr_gv_model_performance_pairs.csv"), row.names=FALSE)
print(perf_tab[,c("model_id","model","c_corrected","brier","calibration_slope_apparent","AIC")])
print(pairs_tab)

# ---- 4. bootstrap 稳定性(J1 系数) ----
stab_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  B <- 1000; n <- nrow(cc); coefs <- matrix(NA, B, 2); fails <- 0L
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE); db <- cc[idx,]
    f <- tryCatch(coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml)), data=db),
                  error=function(e) NULL)
    if (is.null(f)) { fails <- fails+1L; next }
    coefs[b,] <- coef(f)[c("shr_zf","gv_zf")]
  }
  for (k in 1:2) {
    nm <- c("shr_zf","gv_zf")[k]
    x <- coefs[,k]
    stab_rows[[length(stab_rows)+1]] <- data.frame(model_id=paste0("STAB_J1_",hz,"d_",nm), horizon=hz, term=nm,
      boot_mean=mean(x,na.rm=TRUE), boot_sd=sd(x,na.rm=TRUE),
      p2.5=quantile(x,.025,na.rm=TRUE), p97.5=quantile(x,.975,na.rm=TRUE),
      frac_sign_reversal=mean(sign(x)!=sign(mean(x,na.rm=TRUE)), na.rm=TRUE),
      frac_abs_gt_1=mean(abs(x)>1, na.rm=TRUE), n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  }
}
stab <- do.call(rbind, stab_rows)
write.csv(stab, file.path(ROOT,"results","shr_gv_bootstrap_stability.csv"), row.names=FALSE)

# ---- 5. 影响点诊断(J1 CC) ----
fJ1_30 <- coxph(as.formula(paste0("Surv(t_lm_30, event_lm_30) ~ shr_zf + gv_zf + ", covs_fml)), data=cc, x=TRUE, y=TRUE)
dfb <- residuals(fJ1_30, type="dfbeta")
dev <- residuals(fJ1_30, type="deviance")
mar <- residuals(fJ1_30, type="martingale")
infl <- data.frame(stay_id=cc$stay_id,
  max_abs_dfbeta_shr=apply(abs(dfb[, grep("shr", colnames(dfb)), drop=FALSE]), 1, max),
  max_abs_dfbeta_gv=apply(abs(dfb[, grep("gv_zf", colnames(dfb)), drop=FALSE]), 1, max),
  deviance=dev, martingale=mar,
  shr=cc$shr, gv=cc$gv, hba1c=cc$hba1c_pct, event_30d=cc$event_lm_30)
infl <- infl[order(-infl$deviance),]
write.csv(head(infl, 50), file.path(ROOT,"results","shr_gv_influence_diagnostics.csv"), row.names=FALSE)
p3 <- ggplot(infl, aes(x=seq_len(nrow(infl)), y=deviance)) + geom_point(alpha=0.4, size=0.8) +
  labs(title="Deviance residuals (J1, 30d)", x="patient (sorted)", y="deviance residual") + theme_bw()
ggsave(file.path(ROOT,"figures","joint_influence_deviance.png"), p3, width=7, height=4, dpi=600)
cat("PHASE22_DONE\n")
sink()
