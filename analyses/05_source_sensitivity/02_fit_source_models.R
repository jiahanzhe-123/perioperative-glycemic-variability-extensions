# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 07b_source_models.R — 测量来源敏感性模型(同一 Model B 框架;同患者 POCT vs lab 比较)
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
sink(file.path(ROOT,"logs/07b_source_models.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$mean_glu <- d$ps_mean_glucose; d$glucose_count <- d$ps_glucose_count
d$gender <- factor(d$gender)
d$procedure_cat6 <- factor(d$procedure_cat6,
  levels=c("isolated CABG","isolated open valve","combined CABG + open valve",
           "open aortic surgery (+/- other)","transplant/VAD","congenital/other open cardiac"))
K <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
rcs_mean <- paste0("rms::rcs(mean_glu, c(", paste(format(K$mean_glu_knots, digits=10), collapse=","), "))")

run_one <- function(gv_vec, dd, tag, cohort_label){
  dd$gv10 <- gv_vec/10
  rows <- list()
  for (hz in c("30","365")) {
    tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
    dd2 <- dd[complete.cases(dd[,c(tm,ev,"gv10","mean_glu",COVS)]) & dd[[tm]]>0,]
    if (nrow(dd2) < 100 || sum(dd2[[ev]]) < 10) {
      rows[[length(rows)+1]] <- data.frame(model_id=paste0("SRC_",tag,"_",hz,"d"), cohort=cohort_label,
        outcome=paste0("mortality by index day ",hz," among day-1 landmark survivors"),
        N=nrow(dd2), events=sum(dd2[[ev]]), note="insufficient N/events (<100 N or <10 events)", stringsAsFactors=FALSE)
      next
    }
    fit <- coxph(as.formula(paste0("Surv(",tm,",",ev,") ~ gv10 + ", rcs_mean, " + ", covs_fml)), data=dd2)
    s <- summary(fit)
    rows[[length(rows)+1]] <- data.frame(model_id=paste0("SRC_",tag,"_",hz,"d"), cohort=cohort_label,
      outcome=paste0("mortality by index day ",hz," among day-1 landmark survivors"),
      model="Model B (same as primary)", N=nrow(dd2), events=sum(dd2[[ev]]),
      HR_per10=s$coefficients["gv10","exp(coef)"], lo=s$conf.int["gv10","lower .95"],
      hi=s$conf.int["gv10","upper .95"], P=s$coefficients["gv10","Pr(>|z|)"],
      gv_sd_within=sd(gv_vec[dd$stay_id %in% dd2$stay_id], na.rm=TRUE), note="", stringsAsFactors=FALSE)
  }
  do.call(rbind, rows)
}

tgt <- d[d$landmark_eligible==TRUE,]
out <- list()
# 1. 各来源限制性队列(同一 Model B)
for (src in c("poct","centrallab","bloodgas","common")) {
  f <- read.csv(file.path(ROOT,"data",paste0("features_source_",src,".csv")), stringsAsFactors=FALSE)
  dd <- tgt[tgt$stay_id %in% f$stay_id[f$n_src>=2 & !is.na(f$gv)],]
  gv_vec <- f$gv[match(dd$stay_id, f$stay_id)]
  mg <- f$mean_glu[match(dd$stay_id, f$stay_id)]
  dd$mean_glu <- mg
  lbl <- switch(src, poct="POCT-only series", centrallab="central-laboratory-only series",
                bloodgas="blood-gas-only series", common="common-source (POCT + central lab) series")
  out[[src]] <- run_one(gv_vec, dd, src, lbl)
  cat(src, ": N =", nrow(dd), "\n")
}
# 2. 同患者 POCT vs central lab(完全相同患者、协变量、事件)
sp <- read.csv(file.path(ROOT,"data","samepatient_poct_lab.csv"), stringsAsFactors=FALSE)
dd <- tgt[tgt$stay_id %in% sp$stay_id,]
dd <- dd[complete.cases(dd[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365",COVS)]),]
cat("同患者子集 N =", nrow(dd), "; 30d 事件 =", sum(dd$event_lm_30), "; 365d 事件 =", sum(dd$event_lm_365), "\n")
sp2 <- sp[sp$stay_id %in% dd$stay_id,]
for (which in c("poct","lab")) {
  gv_vec <- sp2[[paste0("gv_",which)]]; mg <- sp2[[paste0("mean_glu_",which)]]
  dd$mean_glu <- mg
  out[[paste0("samepatient_",which)]] <- run_one(gv_vec, dd, paste0("samepatient_",which),
    paste0("same patients (N=", nrow(dd), "), ", toupper(which), "-derived GV"))
}
tab <- do.call(rbind, out)
write.csv(tab, file.path(ROOT,"results","measurement_source_results.csv"), row.names=FALSE)
print(tab[, intersect(c("model_id","cohort","N","events","HR_per10","lo_per10","hi_per10","P","note"), names(tab))])
cat("PHASE7B_DONE\n")
sink()
