# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 07_frame_registry.R
# 为 6 个主分析 frame 建立注册表:每个 frame 输出 stay_id 清单与 MD5 hash。
# 任何补充表声明 frame_hash 即可逐位验证"与主模型同一批患者"。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(readr); library(dplyr); library(tools)})
ROOT <- PGV("mimic_work")
OUT  <- file.path(ROOT,"reviewer_issue_verification_20260728/04_frame_registry")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

dat <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/01_unified_glucose_sequence/final_dataset_for_v3_pipeline.csv"),
                check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("","NA","NULL"))
num <- c("shr","gv","age_at_admission","bmi","charlson_comorbidity_index","lactate_postop_first",
         "creat_postop_first","wbc_postop_first","albumin_adm_first","hgb_postop_first",
         "platelets_postop_first","sofa_24h","survival_time_days","postop_30d_death_flag","one_year_death_flag","diabetes")
for (nm in num) dat[[nm]] <- suppressWarnings(as.numeric(as.character(dat[[nm]])))
dat$gender <- factor(dat$gender)
pg <- tolower(trimws(dat$procedure_group_main))
dat$procedure_group_main_model <- factor(ifelse(pg=="cabg","cabg",ifelse(pg=="open_valve","open_valve",
  ifelse(pg=="transplant_vad","transplant_vad","aortic_congenital_other"))),
  levels=c("cabg","open_valve","aortic_congenital_other","transplant_vad"))
dat$event_30d <- as.integer(dat$postop_30d_death_flag==1)
dat$event_1y <- as.integer(dat$one_year_death_flag==1)
dat$t30 <- pmin(dat$survival_time_days,30); dat$t365 <- pmin(dat$survival_time_days,365)
reduced <- "age_at_admission + gender + charlson_comorbidity_index + procedure_group_main_model + lactate_postop_first + creat_postop_first + sofa_24h"
full <- paste(reduced, "+ bmi + wbc_postop_first + albumin_adm_first + hgb_postop_first + platelets_postop_first")

build_primary_frame <- function(dat, cohort_flag, cohort, horizon){
  d <- dat[dat[[cohort_flag]]==1,]
  ev <- if(horizon==30) "event_30d" else "event_1y"; tm <- if(horizon==30) "t30" else "t365"
  covs <- if(cohort=="GV-only") reduced else if(horizon==30) reduced else full
  fe <- if(cohort=="GV-only" || (cohort=="SHR-GV" && horizon==30)) c("bmi","diabetes") else NULL
  expo_vars <- if(cohort=="GV-only") "gv" else c("gv","shr")
  keep <- complete.cases(d[,c(tm,ev,expo_vars,strsplit(covs," \\+ ")[[1]],fe)]) & d[[tm]]>0
  list(d=d[keep,], ev=ev, tm=tm, covs=covs, fe=fe)
}

reg <- list()
for (cf in c("final_v2_A","final_v2_open_core","final_v2_C")) {
  cohort <- c(final_v2_A="SHR-GV",final_v2_open_core="open-core",final_v2_C="GV-only")[cf]
  for (h in c(30,365)) {
    fr <- build_primary_frame(dat, cf, cohort, h)
    ids <- sort(fr$d$stay_id)
    frame_id <- paste0("FRAME_", gsub("-","",cohort), "_", h, "d")
    lst <- file.path(OUT, paste0(frame_id, "_stay_ids.txt"))
    writeLines(as.character(ids), lst)
    h5 <- unname(md5sum(lst))
    reg[[length(reg)+1]] <- data.frame(
      frame_id=frame_id, cohort=cohort, cohort_flag=cf, horizon_days=h,
      N=length(ids), events=sum(fr$d[[fr$ev]]), md5_stay_ids=h5,
      covariates=fr$covs,
      frame_extra=if(is.null(fr$fe)) "(none)" else paste(fr$fe, collapse="+"),
      exposures_required=paste(if(cohort=="GV-only") "gv" else c("gv","shr"), collapse="+"),
      rule="complete-case on outcome+exposure+covariates+frame_extra, survival_time>0; cohort-specific z within frame",
      stringsAsFactors=FALSE)
  }
}
reg <- dplyr::bind_rows(reg)
write_csv(reg, file.path(OUT,"frame_registry.csv"))
print(reg %>% select(frame_id,N,events,md5_stay_ids))
cat("FRAME_REGISTRY_DONE\n")
