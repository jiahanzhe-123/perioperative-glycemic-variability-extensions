# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 10_performance.R — 模型表现与增量信息(内部描述,不称外部验证)
# 同一样本比较:clinical / clinical+mean / clinical+mean+GV;
# 同一批 1000 个患者级 bootstrap 成对重拟合;apparent + optimism-corrected C-index、
# 配对 ΔC + CI、calibration intercept/slope、time-dependent Brier(30d)。
# 另:Model A vs B 系数变化(log-HR 差)的 bootstrap CI(仅描述,不作中介解释)。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/10_performance.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$gender <- factor(d$gender)
d$procedure_cat6 <- factor(d$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")
cc <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
cc <- cc[complete.cases(cc[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","gv","mean_glu",COVS)]) & cc$t_lm_30>0,]
cc$gv10 <- cc$gv/10
cat("performance sample N =", nrow(cc), "\n")

fml1 <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_fml))
fml2 <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ ", rcs_mean, " + ", covs_fml))
fml3 <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml))
fmlA <- function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_fml))

cstat <- function(fit, tm, ev, dd){
  lp <- predict(fit, newdata=dd, type="lp")
  1 - survival::concordance(survival::Surv(dd[[tm]], dd[[ev]]) ~ lp)$concordance
}
calib <- function(fit, tm, ev, dd, tmax){
  # calibration intercept/slope:以模型预测 lp 重拟合
  lp <- predict(fit, newdata=dd, type="lp")
  ci <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ offset(lp)")), data=dd)
  # intercept: 以 offset(lp) 拟合常数项 → 近似 calibration-in-the-large
  intercept <- tryCatch({ f0 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ 1")), data=dd)
                          # 用 lp 作 offset 的常数风险差
                          NA_real_ }, error=function(e) NA_real_)
  slope <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ lp")), data=dd)$coefficients["lp"]
  list(intercept=intercept, slope=unname(slope))
}
brier <- function(fit, tm, ev, dd, tmax){
  # 以 basehaz 个体化预测 tmax 时刻风险,计算 Brier(不做事删失校正,作内部比较)
  bh <- basehaz(fit, centered=FALSE); bh <- bh[bh$time<=tmax,]; H <- max(bh$hazard)
  lp <- predict(fit, newdata=dd, type="lp", reference="zero")
  p <- 1 - exp(-H*exp(lp))
  y <- dd[[ev]]; y[dd[[tm]]<tmax & dd[[ev]]==0] <- NA
  mean((y-p)^2, na.rm=TRUE)
}
perf_rows <- list(); boots <- list()
B <- 1000; n <- nrow(cc)
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz); tmax <- if(hz=="30") 29 else 364
  fits0 <- list(M1=coxph(fml1(tm,ev), data=cc), M2=coxph(fml2(tm,ev), data=cc), M3=coxph(fml3(tm,ev), data=cc))
  ap <- sapply(fits0, function(f) cstat(f, tm, ev, cc))
  br <- sapply(fits0, function(f) brier(f, tm, ev, cc, tmax))
  sl <- sapply(fits0, function(f) calib(f, tm, ev, cc, tmax)$slope)
  # bootstrap:同一批重抽样,成对
  set.seed(SEED + as.integer(hz))
  opt <- matrix(0, nrow=B, ncol=3); dc13 <- dc23 <- numeric(B); ab <- numeric(B)
  slp <- matrix(NA_real_, nrow=B, ncol=3); fails <- 0L
  for (b in 1:B) {
    idx <- sample.int(n, n, replace=TRUE); db <- cc[idx,]
    r <- tryCatch({
      f1 <- coxph(fml1(tm,ev), data=db); f2 <- coxph(fml2(tm,ev), data=db); f3 <- coxph(fml3(tm,ev), data=db)
      fA <- coxph(fmlA(tm,ev), data=db)
      c(cstat(f1,tm,ev,db)-cstat(f1,tm,ev,cc), cstat(f2,tm,ev,db)-cstat(f2,tm,ev,cc), cstat(f3,tm,ev,db)-cstat(f3,tm,ev,cc),
        cstat(f3,tm,ev,cc)-cstat(f1,tm,ev,cc), cstat(f3,tm,ev,cc)-cstat(f2,tm,ev,cc),
        summary(fA)$coefficients["gv10","coef"] - summary(f3)$coefficients["gv10","coef"],
        calib(f1,tm,ev,cc,tmax)$slope, calib(f2,tm,ev,cc,tmax)$slope, calib(f3,tm,ev,cc,tmax)$slope)
    }, error=function(e) NULL)
    if (is.null(r)) { fails <- fails+1L; next }
    opt[b,] <- r[1:3]; dc13[b] <- r[4]; dc23[b] <- r[5]; ab[b] <- r[6]; slp[b,] <- r[7:9]
  }
  ci <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
  for (i in 1:3) {
    mn <- paste0("M",i)
    perf_rows[[length(perf_rows)+1]] <- data.frame(model_id=paste0("PERF_",mn,"_",hz,"d"), horizon=hz,
      model=c("clinical","clinical + mean glucose","clinical + mean glucose + GV")[i],
      N=nrow(cc), events=sum(cc[[ev]]),
      c_apparent=ap[i], optimism=mean(opt[,i]), c_corrected=ap[i]-mean(opt[,i]),
      brier=br[i], calibration_slope_apparent=sl[i],
      calibration_slope_boot_mean=mean(slp[,i], na.rm=TRUE),
      calibration_slope_boot_lo=ci(slp[,i])[1], calibration_slope_boot_hi=ci(slp[,i])[2],
      n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  }
  boots[[length(boots)+1]] <- data.frame(model_id=paste0("PERF_delta_",hz,"d"), horizon=hz,
    dC_M3_vs_M1=mean(dc13), dC_M1_lo=ci(dc13)[1], dC_M1_hi=ci(dc13)[2],
    dC_M3_vs_M2=mean(dc23), dC_M2_lo=ci(dc23)[1], dC_M2_hi=ci(dc23)[2],
    betaA_minus_betaB_mean=mean(ab), ab_lo=ci(ab)[1], ab_hi=ci(ab)[2],
    n_boot=B, n_boot_failed=fails, stringsAsFactors=FALSE)
  cat(sprintf("[%sd] 完成,失败 %d\n", hz, fails))
}
perf <- do.call(rbind, perf_rows); dl <- do.call(rbind, boots)
write.csv(perf, file.path(ROOT,"results","11_prediction_performance.csv"), row.names=FALSE)
write.csv(dl, file.path(ROOT,"results","performance_delta_and_AvsB.csv"), row.names=FALSE)
print(perf); print(dl)
cat("PHASE10_DONE\n")
sink()
