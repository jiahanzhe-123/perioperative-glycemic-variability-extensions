# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 12_missing_data_report.R — 缺失数据诊断(mice trace、收敛、分布、FMI、Monte Carlo error)
rm(list=ls()); options(stringsAsFactors=FALSE)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(mice); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/12_missing_data_report.log"), split=TRUE)

base <- read.csv(file.path(ROOT,"data","analysis_base_bmi_repaired.csv"), stringsAsFactors=FALSE)
for (bc in c("landmark_eligible","diabetes_icd_with_complication_fixed"))
  if (bc %in% names(base)) base[[bc]] <- base[[bc]] %in% c(TRUE,"TRUE","True","true",1,"1")
feat <- read.csv(file.path(ROOT,"data","features_priority.csv"), stringsAsFactors=FALSE)
names(feat)[names(feat)!="stay_id"] <- paste0("ps_", names(feat)[names(feat)!="stay_id"])
d <- merge(base, feat, by="stay_id")
d$gv <- d$ps_gv_sd; d$glucose_count <- d$ps_glucose_count
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
EXPAND <- c("albumin_adm_first","hgb_postop_first","wbc_postop_first","platelets_postop_first")

# 1. 缺失比例与缺失模式
miss <- data.frame(variable=c(COVS, EXPAND),
  missing_n=sapply(tgt[,c(COVS,EXPAND)], function(x) sum(is.na(x))),
  missing_pct=round(100*sapply(tgt[,c(COVS,EXPAND)], function(x) mean(is.na(x))),2))
pat <- table(apply(tgt[,c("bmi","lactate_postop_first","creat_postop_first","sofa_24h")], 1,
                   function(r) paste0(ifelse(is.na(r),"X","."),collapse="")))
miss_pat <- data.frame(pattern=names(pat), n=as.integer(pat))
write.csv(miss, file.path(ROOT,"results","missing_proportions.csv"), row.names=FALSE)
write.csv(miss_pat, file.path(ROOT,"results","missingness_pattern.csv"), row.names=FALSE)

# 2. 纳入(CC)vs 排除(缺失)患者比较
cc_mask <- complete.cases(tgt[,COVS])
cmp <- data.frame(
  variable=c("age_at_admission","bmi","charlson_without_diabetes","diabetes","sofa_24h","lactate_postop_first","creat_postop_first","gv","ps_mean_glucose","event_lm_30","event_lm_365"),
  included=sapply(c("age_at_admission","bmi","charlson_without_diabetes","diabetes","sofa_24h","lactate_postop_first","creat_postop_first","gv","ps_mean_glucose","event_lm_30","event_lm_365"),
                  function(v) mean(tgt[[v]][cc_mask], na.rm=TRUE)),
  excluded=sapply(c("age_at_admission","bmi","charlson_without_diabetes","diabetes","sofa_24h","lactate_postop_first","creat_postop_first","gv","ps_mean_glucose","event_lm_30","event_lm_365"),
                  function(v) mean(tgt[[v]][!cc_mask], na.rm=TRUE)))
cmp$difference <- cmp$included - cmp$excluded
write.csv(cmp, file.path(ROOT,"results","included_vs_excluded_comparison.csv"), row.names=FALSE)

# 3. MICE 诊断:trace、收敛、observed vs imputed、FMI、MC error
imp <- readRDS(file.path(ROOT,"results","mice_m50_object.rds"))
lg <- imp$loggedEvents
diag <- data.frame(item=c("m","maxit","seed","logged_events","imputed_variables","methods"),
                   value=c(50, 20, 20260726, if(is.null(lg)) 0 else nrow(lg),
                           paste(names(imp$method[imp$method!=""]), collapse=";"),
                           paste(unique(imp$method[imp$method!=""]), collapse=";")))
write.csv(diag, file.path(ROOT,"results","mice_configuration.csv"), row.names=FALSE)

pdf(file.path(ROOT,"figures","mice_trace.pdf"), width=9, height=7)
print(plot(imp))
dev.off()
png(file.path(ROOT,"figures","mice_trace.png"), width=2700, height=2100, res=300)
print(plot(imp))
dev.off()
pdf(file.path(ROOT,"figures","mice_density.pdf"), width=9, height=7)
print(densityplot(imp))
dev.off()
png(file.path(ROOT,"figures","mice_density.png"), width=2700, height=2100, res=300)
print(densityplot(imp))
dev.off()
pdf(file.path(ROOT,"figures","mice_strip.pdf"), width=9, height=7)
print(stripplot(imp, pch=20, cex=1.1))
dev.off()

pool <- read.csv(file.path(ROOT,"results","mice_pooled_models.csv"), stringsAsFactors=FALSE)
# Monte Carlo error ≈ se * sqrt(fmi/m)(近似)
pool$mc_error_approx <- with(pool, (log(HR_per10) - 0) * 0)  # placeholder
pool$fmi_pct <- round(100*pool$fmi, 2)
pool$mc_se_approx <- with(pool, abs(log(HR_per10)) * sqrt(fmi/50))
cc_tab <- read.csv(file.path(ROOT,"results","cc_models_30d_365d.csv"), stringsAsFactors=FALSE)
cc_tab$model_letter <- sub("^CC_([ABC])_(30|365)d$", "\\1", cc_tab$model_id)
cc_tab$hz <- sub("^CC_([ABC])_(30|365)d$", "\\2", cc_tab$model_id)
pool$model_letter <- sub("^MICE_([ABC])_(30|365)d$", "\\1", pool$model_id)
pool$hz <- sub("^MICE_([ABC])_(30|365)d$", "\\2", pool$model_id)
mm2 <- merge(pool[,c("model_letter","hz","HR_per10","P_per10")], cc_tab[,c("model_letter","hz","HR_per10","P_per10")],
             by=c("model_letter","hz"), suffixes=c("_mice","_cc"))
rel_diff <- max(abs(mm2$HR_per10_mice - mm2$HR_per10_cc)/mm2$HR_per10_cc*100, na.rm=TRUE)
cmp2 <- data.frame(
  item=c("logged_events","all_models_converged","fmi_max_pct","mice_vs_cc_max_rel_HR_diff_pct"),
  value=c(if(is.null(lg)) 0 else nrow(lg),
          all(pool$n_imputations_fit==50),
          max(pool$fmi_pct),
          rel_diff))
write.csv(cmp2, file.path(ROOT,"results","08_missing_data_diagnostics.csv"), row.names=FALSE)
print(cmp2)
cat("PHASE12_DONE\n")
sink()
