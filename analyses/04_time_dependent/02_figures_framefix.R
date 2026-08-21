# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 08_figures_framefix.R — 用 framefix 数据重生成补充森林图(Figure S15 组件)
# 输入:01_frame_repair 各模块 framefix CSV;输出 PNG 600dpi + PDF。
rm(list=ls()); options(stringsAsFactors=FALSE); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(readr); library(dplyr); library(ggplot2)})
ROOT <- PGV("mimic_work")
RDIR <- file.path(ROOT,"reviewer_issue_verification_20260728")
FIG  <- file.path(RDIR,"figures"); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

mk_forest <- function(df, x, lo, hi, y, label, file, title, xlim=c(0.7,1.8)) {
  p <- ggplot(df, aes(x=.data[[x]], y=.data[[y]])) +
    geom_vline(xintercept=1, linetype=2, color="grey50") +
    geom_errorbarh(aes(xmin=.data[[lo]], xmax=.data[[hi]]), height=0.18) +
    geom_point(size=2.2, color="#0072B2") +
    geom_text(aes(label=.data[[label]]), nudge_y=0.3, size=2.8) +
    scale_x_log10(limits=xlim) +
    labs(x="Adjusted HR per 1-SD increase (log scale)", y=NULL, title=title) +
    theme_bw(base_size=10)
  ggsave(file.path(FIG, paste0(file,".png")), p, width=9, height=max(4, 0.35*nrow(df)+2), dpi=600)
  ggsave(file.path(FIG, paste0(file,".pdf")), p, width=9, height=max(4, 0.35*nrow(df)+2))
  p
}

# ---- Part 2: 替代 GV(clinical vs +mean;含配对 SD 对照行) ----
t2 <- read_csv(file.path(RDIR,"01_frame_repair/alt_gv/alt_gv_models_framefix.csv"), show_col_types=FALSE) %>%
  filter(!is.na(HR_clinical))
d2 <- bind_rows(
  t2 %>% transmute(y=paste0(metric," | ",cohort," | ",horizon_days,"d | clinical"), HR=HR_clinical, lo=lo1, hi=hi1,
                   lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_clinical,lo1,hi1)),
  t2 %>% transmute(y=paste0(metric," | ",cohort," | ",horizon_days,"d | +mean"), HR=HR_plus_mean, lo=lo2, hi=hi2,
                   lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_plus_mean,lo2,hi2)))
d2$y <- factor(d2$y, levels=rev(unique(d2$y)))
mk_forest(d2, "HR","lo","hi","y","lab","part2_alt_gv_forest_framefix","Alternative GV metrics (primary frame): clinical vs + mean glucose")

# ---- Part 3: 限制性分析 ----
t3 <- read_csv(file.path(RDIR,"01_frame_repair/measurement/measurement_restriction_analyses_framefix.csv"), show_col_types=FALSE) %>%
  filter(is.na(note) | note=="", !is.na(HR_gv)) %>%
  mutate(y=factor(paste0(restriction," | ",cohort," | ",horizon_days,"d (N=",N,",ev=",events,")")),
         lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_gv,lo,hi))
d3 <- t3 %>% mutate(y=factor(y, levels=rev(levels(y))))
mk_forest(d3, "HR_gv","lo","hi","y","lab","part3_restrictions_forest_framefix","Measurement-process restriction analyses (primary frame)")

# ---- Part 3: 序列调整 ----
adj_ord <- c("none","count","span","count_span","count_span_source")
t3b <- read_csv(file.path(RDIR,"01_frame_repair/measurement/measurement_sequential_adjustment_framefix.csv"), show_col_types=FALSE) %>%
  mutate(adjustment=factor(adjustment, levels=adj_ord)) %>%
  arrange(cohort, horizon_days, adjustment) %>%
  mutate(y=paste0(adjustment," | ",cohort," | ",horizon_days,"d"),
         lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_gv,lo,hi))
t3b$y <- factor(t3b$y, levels=rev(unique(t3b$y)))
mk_forest(t3b, "HR_gv","lo","hi","y","lab","part3_sequential_forest_framefix","Sequential adjustment for measurement process (primary frame)", xlim=c(0.9,1.55))

# ---- Part 5: 极端血糖负荷调整 ----
ext_ord <- c("none","hypo","hyper","minmax","burden")
t5 <- read_csv(file.path(RDIR,"01_frame_repair/extremes/extremes_extreme_adjustments_framefix.csv"), show_col_types=FALSE) %>%
  filter(is.na(note) | note=="", !is.na(HR_gv)) %>%
  mutate(adjustment=factor(adjustment, levels=ext_ord)) %>%
  arrange(cohort, horizon_days, adjustment) %>%
  mutate(y=paste0(adjustment," | ",cohort," | ",horizon_days,"d"),
         lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_gv,lo,hi))
t5$y <- factor(t5$y, levels=rev(unique(t5$y)))
mk_forest(t5, "HR_gv","lo","hi","y","lab","part5_extremes_forest_framefix","Extreme-glucose burden adjustments (primary frame)", xlim=c(0.9,1.6))

# ---- Part 4: 糖尿病分层 ----
t4 <- read_csv(file.path(RDIR,"01_frame_repair/diabetes_ix/diabetes_ix_interaction_models_framefix.csv"), show_col_types=FALSE) %>%
  filter(db=="MIMIC-IV", exposure=="gv_z")
d4 <- bind_rows(
  t4 %>% transmute(y=paste0(cohort," | ",horizon_days,"d | non-diabetic"), HR=HR_nondiabetic, lo=lo_ndm, hi=hi_ndm,
                   lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_nondiabetic,lo_ndm,hi_ndm)),
  t4 %>% transmute(y=paste0(cohort," | ",horizon_days,"d | diabetic"), HR=HR_diabetic, lo=lo_dm, hi=hi_dm,
                   lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_diabetic,lo_dm,hi_dm)))
d4$y <- factor(d4$y, levels=rev(unique(d4$y)))
mk_forest(d4, "HR","lo","hi","y","lab","part4_diabetes_strata_forest_framefix","GV association by diabetes status (primary frame; interaction P in table)")
cat("FIGURES_FRAMEFIX_DONE\n")
