# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 06_shr_hba1c.R — SHR 独立次要模块 + HbA1c 可用性 IPW(选择偏倚敏感性)
# 主要 SHR 分析仅用索引日前 1–90 日 HbA1c;索引日值不称术前。
# 不把 SHR/mean/HbA1c 三个作同模独立线性暴露;以"mean + HbA1c 分别进入"判别分子/分母驱动。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite); library(rms); library(ggplot2); library(Hmisc)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/06_shr_hba1c.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$eAG <- 28.7 * d$hba1c_pct - 46.7
d$shr <- d$mean_glu / d$eAG
d$gender <- factor(d$gender)
d$procedure_cat6 <- factor(d$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
cat("HbA1c 层级(landmark 目标队列内):\n"); print(table(tgt$hba1c_tier))

empty_row <- function() data.frame(model_id=character(), analysis=character(), cohort=character(),
  outcome=character(), model=character(), N=integer(), events=integer(), term=character(),
  HR=double(), lo=double(), hi=double(), P=double(), overall_p=double(), nonlinear_p=double(),
  knots=character(), scale=character(), note=character(), stringsAsFactors=FALSE)
rows <- list()
add_row <- function(...) { rows[[length(rows)+1]] <<- merge(empty_row(), data.frame(..., stringsAsFactors=FALSE), all=TRUE) }

lrt <- function(a,b){ ddf<-attr(logLik(b),"df")-attr(logLik(a),"df")
  if(!is.finite(ddf)||ddf<=0) return(NA_real_); pchisq(2*as.numeric(logLik(b)-logLik(a)),df=ddf,lower.tail=FALSE) }

run_shr_block <- function(dd, tag, cohort_label){
  if (nrow(dd) < 200) {
    add_row(model_id=paste0("SHR_",tag), analysis="secondary (SHR module)", cohort=cohort_label,
            note=paste0("N<200, not estimable (N=", nrow(dd), ")")); return(invisible(NULL))
  }
  shr_sd <- sd(dd$shr); dd$shr_zf <- (dd$shr - mean(dd$shr))/shr_sd
  dd$mean_zf <- (dd$mean_glu - mean(dd$mean_glu))/sd(dd$mean_glu)
  dd$a1c_zf <- (dd$hba1c_pct - mean(dd$hba1c_pct))/sd(dd$hba1c_pct)
  knots <- quantile(dd$shr, c(.05,.35,.65,.95)); kt <- paste(format(knots, digits=10), collapse=",")
  for (hz in c("30","365")) {
    tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
    oc <- paste0("mortality by index day ",hz," among day-1 landmark survivors")
    # (a) SHR 线性 + 固定临床协变量
    f1 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + ", covs_fml)), data=dd)
    s1 <- summary(f1)
    add_row(model_id=paste0("SHR_",tag,"_",hz,"d_linear"), analysis="secondary (SHR module)", cohort=cohort_label,
      outcome=oc, model="SHR + clinical covariates", N=nrow(dd), events=sum(dd[[ev]]), term="shr_zf",
      HR=s1$coefficients["shr_zf","exp(coef)"], lo=s1$conf.int["shr_zf","lower .95"],
      hi=s1$conf.int["shr_zf","upper .95"], P=s1$coefficients["shr_zf","Pr(>|z|)"],
      scale=paste0("per fixed-cohort SD (SD=", round(shr_sd,4), ")"))
    # (b) SHR 样条
    fS <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(shr, c(", kt, ")) + ", covs_fml)), data=dd, x=TRUE, y=TRUE)
    fL <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr + ", covs_fml)), data=dd, x=TRUE, y=TRUE)
    f0 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_fml)), data=dd, x=TRUE, y=TRUE)
    add_row(model_id=paste0("SHR_",tag,"_",hz,"d_spline"), analysis="secondary (SHR module)", cohort=cohort_label,
      outcome=oc, model="SHR RCS (4 knots)", N=nrow(dd), events=sum(dd[[ev]]),
      overall_p=lrt(f0,fS), nonlinear_p=lrt(fL,fS), knots=paste(round(knots,3),collapse=";"))
    # (c) mean glucose 与 HbA1c 分别进入(分子/分母驱动判别)
    f2 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ mean_zf + a1c_zf + ", covs_fml)), data=dd)
    s2 <- summary(f2)
    for (t2 in c("mean_zf","a1c_zf"))
      add_row(model_id=paste0("SHRdecomp_",tag,"_",hz,"d_",t2), analysis="secondary (SHR module)", cohort=cohort_label,
        outcome=oc, model="mean glucose + HbA1c entered separately", N=nrow(dd), events=sum(dd[[ev]]), term=t2,
        HR=s2$coefficients[t2,"exp(coef)"], lo=s2$conf.int[t2,"lower .95"],
        hi=s2$conf.int[t2,"upper .95"], P=s2$coefficients[t2,"Pr(>|z|)"], scale="per within-cohort SD")
    # (d) 交互(仅探索)
    f3 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ mean_zf*a1c_zf + ", covs_fml)), data=dd)
    s3 <- summary(f3)
    add_row(model_id=paste0("SHRinteract_",tag,"_",hz,"d"), analysis="exploratory", cohort=cohort_label,
      outcome=oc, model="mean glucose x HbA1c interaction", N=nrow(dd), events=sum(dd[[ev]]), term="mean_zf:a1c_zf",
      HR=s3$coefficients["mean_zf:a1c_zf","exp(coef)"], lo=s3$conf.int["mean_zf:a1c_zf","lower .95"],
      hi=s3$conf.int["mean_zf:a1c_zf","upper .95"], P=s3$coefficients["mean_zf:a1c_zf","Pr(>|z|)"])
  }
  invisible(NULL)
}
mk_set <- function(tiers){
  dd <- tgt[tgt$hba1c_tier %in% tiers & !is.na(tgt$eAG) & tgt$eAG > 0,]
  dd[complete.cases(dd[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","shr","gv","mean_glu",COVS)]) & dd$t_lm_30>0,]
}
run_shr_block(mk_set("1-90 days pre-index"), "pre90", "SHR cohort (HbA1c 1-90d pre-index)")
run_shr_block(mk_set(c("1-90 days pre-index","91-365 days pre-index")), "pre365", "SHR cohort (strictly pre-index 1-365d)")
run_shr_block(mk_set(c("1-90 days pre-index","91-365 days pre-index","index calendar day")), "hier", "SHR cohort (original hierarchical rule)")
run_shr_block(mk_set("91-365 days pre-index"), "pre91365", "SHR cohort (91-365d only)")
shr_tab <- do.call(rbind, rows)
write.csv(shr_tab, file.path(ROOT,"results","shr_module_results.csv"), row.names=FALSE)
print(shr_tab[, c("model_id","cohort","model","N","events","term","HR","lo","hi","P","overall_p","nonlinear_p","note")], max=60)

# ---- IPW(HbA1c availability,稳定化权重,99% 截断为主方案) ----
ipw_d <- tgt[!is.na(tgt$hba1c_tier),]
ipw_d$A <- as.integer(!is.na(ipw_d$hba1c_pct))
ipw_d <- ipw_d[complete.cases(ipw_d[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","gv","mean_glu",COVS)]) & ipw_d$t_lm_30>0,]
ipw_d$gv10 <- ipw_d$gv/10
ps_fit <- glm(A ~ age_at_admission + gender + bmi + diabetes + charlson_without_diabetes +
              procedure_cat6 + lactate_postop_first + creat_postop_first + sofa_24h,
              data=ipw_d, family=binomial)
ipw_d$p_hat <- predict(ps_fit, type="response")
ipw_d$w_stab <- mean(ipw_d$A==1) / ipw_d$p_hat
w1 <- ipw_d[ipw_d$A==1,]

smd_of <- function(dd, w){
  out <- list()
  for (v in c("age_at_admission","bmi","diabetes","charlson_without_diabetes","lactate_postop_first","creat_postop_first","sofa_24h")) {
    x1 <- dd[dd$A==1,v]; x0 <- dd[dd$A==0,v]; w1v <- w[dd$A==1]
    m1 <- weighted.mean(x1, w1v); m0 <- mean(x0)
    s1 <- sqrt(wtd.var(x1, w1v)); s0 <- sd(x0)
    out[[v]] <- (m1-m0)/sqrt((s1^2+s0^2)/2)
  }
  out
}
ipw_rows <- list()
for (tr in c(1, 0.995, 0.99, 0.975, 0.95)) {
  q <- if(tr==1) Inf else unname(quantile(w1$w_stab, tr))
  w_use <- pmin(w1$w_stab, q)
  ess <- sum(w_use)^2 / sum(w_use^2)
  for (hz in c("30","365")) {
    tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
    f <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", covs_fml)), data=w1, weights=w_use, robust=TRUE)
    s <- summary(f)
    ipw_rows[[length(ipw_rows)+1]] <- data.frame(model_id=sprintf("IPW_%s_%sd", if(tr==1) "untrunc" else paste0("tr",tr*100), hz),
      analysis="selection-sensitivity (IPW)", truncation=tr, N=nrow(w1), events=sum(w1[[ev]]),
      weighted_events=sum(w_use*w1[[ev]]), ess=ess, ess_ratio=ess/nrow(w1), max_weight=max(w_use),
      HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
      hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"], stringsAsFactors=FALSE)
  }
}
ipw_tab <- do.call(rbind, ipw_rows)
smd99 <- smd_of(ipw_d, pmin(ipw_d$w_stab, quantile(w1$w_stab, .99)))
smd_unw <- smd_of(ipw_d, rep(1, nrow(ipw_d)))
smd_df <- data.frame(variable=names(smd99), smd_unweighted=unlist(smd_unw), smd_weighted_99=unlist(smd99))
w_sum <- data.frame(item=c("p_hat_min","p_hat_p1","p_hat_median","p_hat_p99","p_hat_max",
                           "w_min","w_median","w_p95","w_p99","w_max","ess_untruncated","ess_ratio_untruncated"),
                    value=c(min(ipw_d$p_hat), quantile(ipw_d$p_hat,.01), median(ipw_d$p_hat), quantile(ipw_d$p_hat,.99), max(ipw_d$p_hat),
                            min(w1$w_stab), median(w1$w_stab), unname(quantile(w1$w_stab,.95)), unname(quantile(w1$w_stab,.99)), max(w1$w_stab),
                            sum(w1$w_stab)^2/sum(w1$w_stab^2), (sum(w1$w_stab)^2/sum(w1$w_stab^2))/nrow(w1)))
write.csv(ipw_tab, file.path(ROOT,"results","09_ipw_diagnostics.csv"), row.names=FALSE)
write.csv(smd_df, file.path(ROOT,"results","ipw_balance_smd.csv"), row.names=FALSE)
write.csv(w_sum, file.path(ROOT,"results","ipw_weight_summary.csv"), row.names=FALSE)
p <- ggplot(ipw_d, aes(x=p_hat, fill=factor(A))) + geom_histogram(alpha=0.5, position="identity", bins=50) +
  labs(x="P(HbA1c available)", y="count", fill="A", title="Propensity overlap (HbA1c availability)") + theme_bw()
ggsave(file.path(ROOT,"figures","ipw_positivity.png"), p, width=7, height=4.5, dpi=600)
ggsave(file.path(ROOT,"figures","ipw_positivity.pdf"), p, width=7, height=4.5)
max_smd99 <- max(abs(smd_df$smd_weighted_99), na.rm=TRUE)
insufficient <- max_smd99 > 0.1 || w_sum$value[w_sum$item=="ess_ratio_untruncated"] < 0.5
cat(sprintf("IPW: 99%% 截断后 max|SMD|=%.3f; ESS(untrunc)/N=%.3f; 判定: %s\n",
    max_smd99, w_sum$value[w_sum$item=="ess_ratio_untruncated"],
    if(insufficient) "不足以校正 HbA1c 选择机制(结果仅作诊断)" else "平衡可接受(仍仅作选择偏倚敏感性)"))
cat("PHASE6_DONE\n")
sink()
