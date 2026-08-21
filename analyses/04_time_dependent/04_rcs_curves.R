# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 09_rcs_curves_framefix.R
# 用主模型 frame 重新生成全部 RCS 曲线(Figure S2/S3/S13 及 figure3 替换件)。
# 数据:final blood-only 数据集;frame 与 01_refit_primary_models.R 完全一致(逐模型断言)。
# 曲线在原始暴露尺度(SHR 比值 / GV mg/dL)上拟合与展示,与既有发布图一致;
# LRT 的 overall/nonlinear P 对标 z 尺度 RCS 表(尺度不变,断言 |Δp|<1e-8)。
# 输出:每条曲线的源数据 CSV、合并源数据、Figure S13 源数据(含 TV 面板)、
#       PNG 600dpi + PDF、图注 MD。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(survival); library(readr); library(dplyr); library(rms); library(ggplot2); library(patchwork)})
ROOT <- PGV("mimic_work")
RDIR <- file.path(ROOT,"reviewer_issue_verification_20260728")
OUT  <- file.path(RDIR,"07_rcs_curves")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE); dir.create(file.path(OUT,"source_data"), showWarnings=FALSE)
sink(file.path(RDIR,"logs/09_rcs_curves.log"), split=TRUE)

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
safe_z <- function(x) as.numeric(scale(suppressWarnings(as.numeric(x))))

build_primary_frame <- function(dat, cohort_flag, cohort, horizon){
  d <- dat[dat[[cohort_flag]]==1,]
  ev <- if(horizon==30) "event_30d" else "event_1y"; tm <- if(horizon==30) "t30" else "t365"
  covs <- if(cohort=="GV-only") reduced else if(horizon==30) reduced else full
  fe <- if(cohort=="GV-only" || (cohort=="SHR-GV" && horizon==30)) c("bmi","diabetes") else NULL
  expo_vars <- if(cohort=="GV-only") "gv" else c("gv","shr")
  keep <- complete.cases(d[,c(tm,ev,expo_vars,strsplit(covs," \\+ ")[[1]],fe)]) & d[[tm]]>0
  list(d=d[keep,], ev=ev, tm=tm, covs=covs, fe=fe)
}
prim <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/02_primary_models_refit/final_primary_models.csv"), stringsAsFactors=FALSE)
rcs_ref <- read.csv(file.path(RDIR,"03_rcs_framefix/rcs_primary_frame_framefix.csv"), stringsAsFactors=FALSE)
lrt_p <- function(a,b){
  la <- logLik(a); lb <- logLik(b); ddf <- attr(lb,"df")-attr(la,"df")
  if (!is.finite(ddf) || ddf<=0) return(NA_real_)
  pchisq(2*(as.numeric(lb)-as.numeric(la)), df=ddf, lower.tail=FALSE)
}

# ---- 单条曲线:raw 尺度 spline Cox + 网格预测(对照参考值=中位数) ----
fit_curve <- function(fr, cohort, horizon, exposure){
  d <- fr$d
  other <- if(cohort=="GV-only") character(0) else setdiff(c("shr","gv"), exposure)
  for (o in other) d[[paste0(o,"_z")]] <- safe_z(d[[o]])
  x <- d[[exposure]]
  knots <- as.numeric(quantile(x, probs=c(.05,.35,.65,.95), names=FALSE, type=7))
  ref  <- median(x)
  mu <- mean(x); sg <- sd(x)
  knot_text <- paste(format(knots, scientific=FALSE, digits=12), collapse=",")
  terms_base <- c(if(length(other)) paste0(other,"_z") else NULL, strsplit(fr$covs," \\+ ")[[1]])
  spline_fml <- as.formula(paste("Surv(",fr$tm,",",fr$ev,") ~ rcs(",exposure,", c(",knot_text,")) +", paste(terms_base, collapse=" + ")))
  linear_fml <- as.formula(paste("Surv(",fr$tm,",",fr$ev,") ~ ",exposure," + ", paste(terms_base, collapse=" + ")))
  noexp_fml  <- as.formula(paste("Surv(",fr$tm,",",fr$ev,") ~ ", paste(terms_base, collapse=" + ")))
  fs <- coxph(spline_fml, data=d, x=TRUE, y=TRUE, model=TRUE, ties="efron")
  fl <- coxph(linear_fml, data=d, x=TRUE, y=TRUE, model=TRUE, ties="efron")
  f0 <- coxph(noexp_fml,  data=d, x=TRUE, y=TRUE, model=TRUE, ties="efron")
  overall_p <- lrt_p(f0, fs); nonlinear_p <- lrt_p(fl, fs)
  # 与 z 尺度 RCS 表核对(尺度不变性)
  refrow <- rcs_ref[rcs_ref$cohort==cohort & rcs_ref$horizon_days==horizon &
                    rcs_ref$exposure==(if(exposure=="shr") "SHR" else "GV"),]
  stopifnot(nrow(refrow)==1, nrow(d)==refrow$N, sum(d[[fr$ev]])==refrow$events)
  if (abs(overall_p-refrow$overall_p)>1e-8 || abs(nonlinear_p-refrow$nonlinear_p)>1e-8)
    stop(paste("SCALE-INVARIANCE CHECK FAILED:", cohort, horizon, exposure))
  # 网格预测
  sc <- which(attr(fs$x,"assign")==1)   # spline 基列
  b <- coef(fs)[sc]; V <- vcov(fs)[sc,sc,drop=FALSE]
  grid <- seq(quantile(x,.01,names=FALSE), quantile(x,.99,names=FALSE), length.out=200)
  Bg <- rcs(grid, knots); Br <- rcs(ref, knots)
  D  <- sweep(Bg, 2, Br, "-")
  lp <- as.vector(D %*% b); se <- sqrt(pmax(0, rowSums((D %*% V) * D)))
  curve <- data.frame(cohort=cohort, horizon_days=horizon,
    endpoint=if(horizon==30) "30-day mortality" else "1-year mortality",
    exposure=if(exposure=="shr") "SHR" else "GV", exposure_var=exposure,
    x_raw=grid, x_z=(grid-mu)/sg, HR=exp(lp), lo=exp(lp-1.96*se), hi=exp(lp+1.96*se),
    reference_raw=ref, reference_z=(ref-mu)/sg, N=nrow(d), events=sum(d[[fr$ev]]),
    overall_p=overall_p, nonlinear_p=nonlinear_p,
    knot_probabilities="0.05;0.35;0.65;0.95",
    knot_values_raw=paste(format(knots,digits=10),collapse=";"),
    knot_values_z=paste(format((knots-mu)/sg,digits=10),collapse=";"),
    exposure_mean=mu, exposure_sd=sg,
    covariate_specification=if(fr$covs==reduced) "M3_reduced" else "M3_full",
    shr_jointly_adjusted=cohort!="GV-only",
    standardization_sample=paste0("z within analytic sample (N=",nrow(d),"); knots within same sample"),
    missing_data_rule="primary complete-case frame (see frame registry)", stringsAsFactors=FALSE)
  list(curve=curve, rug=data.frame(x=x, cohort=cohort, horizon_days=horizon,
                                   exposure=if(exposure=="shr") "SHR" else "GV"),
       meta=curve[1, setdiff(names(curve), c("x_raw","x_z","HR","lo","hi"))])
}

jobs <- list()
for (cf in c("final_v2_A","final_v2_open_core","final_v2_C")) {
  cohort <- c(final_v2_A="SHR-GV",final_v2_open_core="open-core",final_v2_C="GV-only")[cf]
  exposures <- if(cohort=="GV-only") "gv" else c("shr","gv")
  for (h in c(30,365)) {
    fr <- build_primary_frame(dat, cf, cohort, h)
    refp <- prim[prim$cohort==cohort & prim$horizon_days==h,]
    stopifnot(nrow(fr$d)==refp$N, sum(fr$d[[fr$ev]])==refp$events)  # frame 断言
    for (ex in exposures) jobs[[length(jobs)+1]] <- list(fr=fr, cohort=cohort, horizon=h, exposure=ex)
  }
}
fits <- lapply(jobs, function(j) fit_curve(j$fr, j$cohort, j$horizon, j$exposure))
curves <- dplyr::bind_rows(lapply(fits, `[[`, "curve"))
rugs   <- dplyr::bind_rows(lapply(fits, `[[`, "rug"))
metas  <- dplyr::bind_rows(lapply(fits, `[[`, "meta"))
write_csv(curves, file.path(OUT,"source_data","rcs_curve_source_data_all.csv"))
write_csv(metas,  file.path(OUT,"source_data","rcs_curve_metadata.csv"))
cat("曲线条数:", nrow(metas), "; 源数据行:", nrow(curves), "\n")

# ---- 作图 ----
COL <- c(GV="#0E7C7B", SHR="#B36500")
mk_panel <- function(cv, rg, ttl){
  ex <- unique(cv$exposure); col <- COL[[ex]]
  ggplot(cv, aes(x=x_raw)) +
    geom_rug(data=rg, aes(x=x), inherit.aes=FALSE, alpha=0.04, length=unit(0.05,"npc"), color=col) +
    geom_ribbon(aes(ymin=lo, ymax=hi), fill=col, alpha=0.18) +
    geom_line(aes(y=HR), color=col, linewidth=0.9) +
    geom_hline(yintercept=1, linetype="dashed", color="grey55", linewidth=0.5) +
    geom_vline(xintercept=unique(cv$reference_raw), linetype="dotted", color="grey40", linewidth=0.5) +
    geom_point(aes(x=reference_raw, y=1), color=col, size=2.2) +
    annotate("text", x=-Inf, y=Inf, hjust=-0.08, vjust=1.6, size=3.1, color="grey35",
             label=sprintf("overall P=%s\nnonlinear P=%s",
                           formatC(unique(cv$overall_p), format="f", digits=3),
                           formatC(unique(cv$nonlinear_p), format="f", digits=3))) +
    labs(title=ttl, x=if(ex=="GV") "GV (SD of 24-h glucose, mg/dL)" else "SHR (ratio)",
         y="Adjusted HR (95% CI)") +
    theme_bw(base_size=10.5) +
    theme(plot.title=element_text(size=10, face="bold"), panel.grid.minor=element_blank())
}
panel_title <- function(letter, cohort, horizon, N, events)
  sprintf("%s  %s, %s (N=%s, events=%s)", letter, cohort,
          if(horizon==30) "30-day" else "1-year", format(N, big.mark=","), events)
get_cv <- function(cohort, horizon, exposure) curves %>% filter(.data$cohort==.env$cohort, horizon_days==.env$horizon, .data$exposure==.env$exposure)
get_rg <- function(cohort, horizon, exposure) rugs   %>% filter(.data$cohort==.env$cohort, horizon_days==.env$horizon, .data$exposure==.env$exposure)

save_fig <- function(p, name, w, h){
  ggsave(file.path(OUT, paste0(name,".png")), p, width=w, height=h, dpi=600)
  ggsave(file.path(OUT, paste0(name,".pdf")), p, width=w, height=h)
}

# --- Figure S3 替换件:SHR spline(SHR-GV 30d/1y + open-core 30d/1y) ---
pS3 <- (mk_panel(get_cv("SHR-GV",30,"SHR"), get_rg("SHR-GV",30,"SHR"),
                 panel_title("A","SHR–GV cohort",30,metas$N[metas$cohort=="SHR-GV"&metas$horizon_days==30&metas$exposure=="SHR"],metas$events[metas$cohort=="SHR-GV"&metas$horizon_days==30&metas$exposure=="SHR"])) +
        mk_panel(get_cv("SHR-GV",365,"SHR"), get_rg("SHR-GV",365,"SHR"),
                 panel_title("B","SHR–GV cohort",365,4381,279))) /
       (mk_panel(get_cv("open-core",30,"SHR"), get_rg("open-core",30,"SHR"),
                 panel_title("C","Open-core cohort",30,4757,99)) +
        mk_panel(get_cv("open-core",365,"SHR"), get_rg("open-core",365,"SHR"),
                 panel_title("D","Open-core cohort",365,4162,254)))
save_fig(pS3, "figure_S3_shr_rcs_spline_framefix", 10, 8)

# --- figure3 替换件:GV spline(三队列 × 两结局,6 面板) ---
gv_meta <- function(co,h) metas[metas$cohort==co & metas$horizon_days==h & metas$exposure=="GV",]
pG <- (mk_panel(get_cv("SHR-GV",30,"GV"), get_rg("SHR-GV",30,"GV"), panel_title("A","SHR–GV cohort",30,gv_meta("SHR-GV",30)$N,gv_meta("SHR-GV",30)$events)) +
       mk_panel(get_cv("SHR-GV",365,"GV"), get_rg("SHR-GV",365,"GV"), panel_title("B","SHR–GV cohort",365,4381,279)) +
       mk_panel(get_cv("open-core",30,"GV"), get_rg("open-core",30,"GV"), panel_title("C","Open-core cohort",30,4757,99)) +
       mk_panel(get_cv("open-core",365,"GV"), get_rg("open-core",365,"GV"), panel_title("D","Open-core cohort",365,4162,254)) +
       mk_panel(get_cv("GV-only",30,"GV"), get_rg("GV-only",30,"GV"), panel_title("E","GV-only cohort",30,9786,256)) +
       mk_panel(get_cv("GV-only",365,"GV"), get_rg("GV-only",365,"GV"), panel_title("F","GV-only cohort",365,9786,650))) +
      plot_layout(ncol=2)
save_fig(pG, "figure_gv_rcs_spline_framefix", 10, 11)

# --- Figure S13 替换件:SHR-GV SHR spline 30d/1y + TV 面板(SHR 1y,SHR-GV 与 open-core) ---
tv <- read.csv(file.path(ROOT,"final_statistical_freeze_20260727/03_dependent_analyses_refit/v3_pipeline/tables/final_v3_timevarying.csv"), stringsAsFactors=FALSE)
tv <- tv[tv$horizon_days==365 & tv$exposure=="shr_z" & tv$cohort %in% c("A","open_core"),]
tv$cohort <- c("A"="SHR-GV","open_core"="open-core")[tv$cohort]
tv$interval <- factor(tv$time_interval, levels=c("0-1 days","2-7 days","8-30 days","31-365 days"))
tv$sparse <- grepl("sparse", tv$interval_event_status)
mk_tv_panel <- function(df, ttl){
  col <- COL[["SHR"]]
  ggplot(df, aes(y=interval)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey55", linewidth=0.5) +
    geom_errorbarh(aes(xmin=lower95, xmax=upper95), height=0.16, color=col, linewidth=0.7) +
    geom_point(aes(x=HR, shape=sparse), color=col, size=2.6) +
    scale_shape_manual(values=c(`FALSE`=15, `TRUE`=0), guide="none") +
    geom_text(aes(x=upper95, label=sprintf("%.2f (%.2f\u2013%.2f)",HR,lower95,upper95)),
              hjust=-0.15, size=3.0, color="grey25") +
    annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=1.4, size=3.1, color="grey35",
             label=sprintf("time-interaction P=%s", formatC(unique(df$time_interaction_p), format="f", digits=3))) +
    scale_x_log10() + coord_cartesian(xlim=c(0.5,4.2)) +
    labs(title=ttl, x="HR per SD (log scale)", y=NULL) +
    theme_bw(base_size=10.5) +
    theme(plot.title=element_text(size=10, face="bold"), panel.grid.minor=element_blank())
}
tv_a <- tv[tv$cohort=="SHR-GV",]; tv_o <- tv[tv$cohort=="open-core",]
pS13 <- (mk_panel(get_cv("SHR-GV",30,"SHR"), get_rg("SHR-GV",30,"SHR"),
                  panel_title("A","SHR–GV cohort",30,metas$N[metas$cohort=="SHR-GV"&metas$horizon_days==30&metas$exposure=="SHR"],114)) +
         mk_panel(get_cv("SHR-GV",365,"SHR"), get_rg("SHR-GV",365,"SHR"),
                  panel_title("B","SHR–GV cohort",365,4381,279)) +
         mk_tv_panel(tv_a, "C  SHR–GV cohort, 1-year (time-varying SHR)") +
         mk_tv_panel(tv_o, "D  Open-core cohort, 1-year (time-varying SHR)")) +
        plot_layout(ncol=2) +
        plot_annotation(caption="Descriptive time-varying association; not a treatment time window. Date-anchored exposure window; see Methods.",
                        theme=theme(plot.caption=element_text(size=8.5, color="grey45", hjust=0)))
save_fig(pS13, "figure_S13_primary_spline_timevarying_framefix", 11, 8.5)
write_csv(tv, file.path(OUT,"source_data","figure_S13_timevarying_source_data.csv"))
cat("RCS_CURVES_DONE\n")
sink()
