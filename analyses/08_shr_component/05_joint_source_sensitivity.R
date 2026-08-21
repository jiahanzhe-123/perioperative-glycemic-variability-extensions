# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 23_joint_source.R — SHR–GV 联合模块测量来源敏感性(固定联合常数,不随来源重算)
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
sink(file.path(ROOT,"logs/23_joint_source.log"), split=TRUE)

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
jcc <- j[complete.cases(j[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","shr","gv","mean_glu",COVS)]) & j$t_lm_30>0,]

rows <- list()
add <- function(...) { rows[[length(rows)+1]] <<- data.frame(...) }
run_j1 <- function(dd, gv_vec, shr_vec, tag, label){
  dd$gv_zf  <- (gv_vec  - JK$gv_mean )/JK$gv_sd
  dd$shr_zf <- (shr_vec - JK$shr_mean)/JK$shr_sd
  for (hz in c("30","365")) {
    tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
    dd2 <- dd[complete.cases(dd[,c(tm,ev,"gv_zf","shr_zf",COVS)]) & dd[[tm]]>0,]
    if (nrow(dd2) < 100 || sum(dd2[[ev]]) < 10) {
      add(model_id=paste0("SRCJ_",tag,"_",hz,"d"), cohort=label, N=nrow(dd2), events=sum(dd2[[ev]]),
        note="insufficient N/events (<100 N or <10 events)")
      next
    }
    f <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml)), data=dd2)
    s <- summary(f)
    # 联合 Wald(2 df,CC)
    V <- vcov(f)[c("shr_zf","gv_zf"),c("shr_zf","gv_zf")]
    b <- coef(f)[c("shr_zf","gv_zf")]
    wstat <- as.numeric(t(b) %*% solve(V) %*% b) / 2
    add(model_id=paste0("SRCJ_",tag,"_",hz,"d"), cohort=label, N=nrow(dd2), events=sum(dd2[[ev]]),
      joint_stat=wstat, joint_df=2, joint_P=pchisq(wstat*2, df=2, lower.tail=FALSE),
      HR_shr=s$coefficients["shr_zf","exp(coef)"], shr_lo=s$conf.int["shr_zf","lower .95"],
      shr_hi=s$conf.int["shr_zf","upper .95"], P_shr=s$coefficients["shr_zf","Pr(>|z|)"],
      HR_gv=s$coefficients["gv_zf","exp(coef)"], gv_lo=s$conf.int["gv_zf","lower .95"],
      gv_hi=s$conf.int["gv_zf","upper .95"], P_gv=s$coefficients["gv_zf","Pr(>|z|)"])
  }
}
# 各来源 J1
for (src in c("poct","centrallab","bloodgas","common")) {
  f <- read.csv(file.path(ROOT,"data",paste0("features_source_",src,".csv")), stringsAsFactors=FALSE)
  dd <- jcc[jcc$stay_id %in% f$stay_id[f$n_src>=2 & !is.na(f$gv)],]
  gv_vec <- f$gv[match(dd$stay_id, f$stay_id)]
  shr_vec <- f$mean_glu[match(dd$stay_id, f$stay_id)] / dd$eAG
  lbl <- switch(src, poct="POCT-only", centrallab="central-laboratory-only",
                bloodgas="blood-gas-only", common="common-source (POCT+central lab)")
  run_j1(dd, gv_vec, shr_vec, src, lbl)
  cat(src, ": N =", nrow(dd), "\n")
}
# all-source same-minute median(用 allmedian 特征)
f <- read.csv(file.path(ROOT,"data","features_allmedian.csv"), stringsAsFactors=FALSE)
dd <- jcc[jcc$stay_id %in% f$stay_id[!is.na(f$gv_sd)],]
gv_vec <- f$gv_sd[match(dd$stay_id, f$stay_id)]
shr_vec <- f$mean_glucose[match(dd$stay_id, f$stay_id)] / dd$eAG
run_j1(dd, gv_vec, shr_vec, "allmedian", "all-source same-minute median")

# measurement-process adjusted(在 priority J1 上加测量过程协变量)
jcc$log_count <- log(jcc$glucose_count)
jcc$gv_zf <- (jcc$gv - JK$gv_mean)/JK$gv_sd; jcc$shr_zf <- (jcc$shr - JK$shr_mean)/JK$shr_sd
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  f2 <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml,
           " + log_count + ps_span_hours + ps_frac_central_lab + ps_frac_poct + ps_frac_icu_charted")), data=jcc)
  s <- summary(f2)
  V <- vcov(f2)[c("shr_zf","gv_zf"),c("shr_zf","gv_zf")]; b <- coef(f2)[c("shr_zf","gv_zf")]
  wstat <- as.numeric(t(b) %*% solve(V) %*% b) / 2
  add(model_id=paste0("SRCJ_process_adj_",hz,"d"), cohort="priority mixed-source + measurement-process adjustment",
    N=nrow(jcc), events=sum(jcc[[ev]]), joint_stat=wstat, joint_df=2, joint_P=pchisq(wstat*2, df=2, lower.tail=FALSE),
    HR_shr=s$coefficients["shr_zf","exp(coef)"], shr_lo=s$conf.int["shr_zf","lower .95"],
    shr_hi=s$conf.int["shr_zf","upper .95"], P_shr=s$coefficients["shr_zf","Pr(>|z|)"],
    HR_gv=s$coefficients["gv_zf","exp(coef)"], gv_lo=s$conf.int["gv_zf","lower .95"],
    gv_hi=s$conf.int["gv_zf","upper .95"], P_gv=s$coefficients["gv_zf","Pr(>|z|)"])
}

# 同患者 POCT vs lab(严格 HbA1c ∩ 双来源充足)
sp <- read.csv(file.path(ROOT,"data","samepatient_poct_lab.csv"), stringsAsFactors=FALSE)
dd <- jcc[jcc$stay_id %in% sp$stay_id,]
dd <- dd[complete.cases(dd[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365",COVS)]),]
sp2 <- sp[sp$stay_id %in% dd$stay_id,]
cat("同患者联合子集 N =", nrow(dd), "; 30d 事件 =", sum(dd$event_lm_30), "; 365d 事件 =", sum(dd$event_lm_365), "\n")
run_j1(dd, sp2$gv_poct, sp2$mean_glu_poct/dd$eAG, "samepatient_poct", "same patients, POCT-derived SHR+GV")
run_j1(dd, sp2$gv_lab,  sp2$mean_glu_lab /dd$eAG, "samepatient_lab",  "same patients, lab-derived SHR+GV")

tab <- do.call(rbind, rows)
write.csv(tab, file.path(ROOT,"results","shr_gv_source_sensitivity.csv"), row.names=FALSE)
print(tab[, intersect(c("model_id","cohort","N","events","joint_P","HR_shr","P_shr","HR_gv","P_gv","note"), names(tab))])
cat("PHASE23_DONE\n")
sink()
