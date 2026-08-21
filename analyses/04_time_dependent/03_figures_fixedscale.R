# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 08b_figures_fixedscale.R — 用统一固定尺度后的数据重画森林图 part A/B/E
# part C(限制性,子集尺度)与 part D(糖尿病交互,本轮未改)不变。
rm(list=ls()); options(stringsAsFactors=FALSE); set.seed(20260726L)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(readr); library(dplyr); library(ggplot2)})
ROOT <- PGV("mimic_work")
FS   <- file.path(ROOT,"reviewer_issue_verification_20260728/08_fixed_scale")
FIG  <- file.path(PGV("manuscript_work"),"UPDATED_FIGURES")

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

# ---- Part A: 替代 GV(SD 全 frame 行 = 固定尺度;其余行 = 各自样本内尺度) ----
t2 <- read_csv(file.path(FS,"s26_alt_gv_fixedscale_machine.csv"), show_col_types=FALSE) %>%
  filter(!is.na(HR_clinical)) %>%
  mutate(fixed = (metric=="SD (mg/dL)" & min_n==2),
         HRc  = ifelse(fixed, HR_clinical_fixed, HR_clinical),
         loc  = ifelse(fixed, lo1_fixed, lo1),
         hic  = ifelse(fixed, hi1_fixed, hi1),
         HRm  = ifelse(fixed, HR_plus_mean_fixed, HR_plus_mean),
         lom  = ifelse(fixed, lo2_fixed, lo2),
         him  = ifelse(fixed, hi2_fixed, hi2))
d2 <- bind_rows(
  t2 %>% transmute(y=paste0(metric," | ",cohort," | ",horizon_days,"d | clinical"), HR=HRc, lo=loc, hi=hic,
                   lab=sprintf("%.2f (%.2f\u2013%.2f)",HRc,loc,hic)),
  t2 %>% transmute(y=paste0(metric," | ",cohort," | ",horizon_days,"d | +mean"), HR=HRm, lo=lom, hi=him,
                   lab=sprintf("%.2f (%.2f\u2013%.2f)",HRm,lom,him)))
d2$y <- factor(d2$y, levels=rev(unique(d2$y)))
mk_forest(d2, "HR","lo","hi","y","lab","part2_alt_gv_forest_fixedscale",
          "Alternative GV metrics: clinical vs + mean glucose (SD rows: fixed primary SD)")

# ---- Part B: 测量过程序列调整(全部固定尺度) ----
t3 <- read_csv(file.path(FS,"s27a_sequential_fixedscale_machine.csv"), show_col_types=FALSE) %>%
  mutate(adjustment=factor(adjustment, levels=c("none","count","span","count_span","count_span_source"))) %>%
  arrange(cohort, horizon_days, adjustment) %>%
  mutate(y=paste0(adjustment," | ",cohort," | ",horizon_days,"d"),
         lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_gv_fixed,lo_fixed,hi_fixed))
t3$y <- factor(t3$y, levels=rev(unique(t3$y)))
mk_forest(t3, "HR_gv_fixed","lo_fixed","hi_fixed","y","lab","part3_sequential_forest_fixedscale",
          "Sequential adjustment for measurement process (fixed primary SD)", xlim=c(0.9,1.55))

# ---- Part E: 极端血糖负荷(全部固定尺度) ----
t5 <- read_csv(file.path(FS,"s28_extremes_fixedscale_machine.csv"), show_col_types=FALSE) %>%
  filter(is.na(note) | note=="", !is.na(HR_gv)) %>%
  mutate(adjustment=factor(adjustment, levels=c("none","hypo","hyper","minmax","burden"))) %>%
  arrange(cohort, horizon_days, adjustment) %>%
  mutate(y=paste0(adjustment," | ",cohort," | ",horizon_days,"d"),
         lab=sprintf("%.2f (%.2f\u2013%.2f)",HR_gv_fixed,lo_fixed,hi_fixed))
t5$y <- factor(t5$y, levels=rev(unique(t5$y)))
mk_forest(t5, "HR_gv_fixed","lo_fixed","hi_fixed","y","lab","part5_extremes_forest_fixedscale",
          "Extreme-glucose burden adjustments (fixed primary SD)", xlim=c(0.9,1.6))
cat("FIGURES_FIXEDSCALE_DONE\n")
