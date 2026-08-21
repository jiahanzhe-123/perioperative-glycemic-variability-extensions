# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# Tables + figures for harmonized replication
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(readr); library(dplyr); library(ggplot2); library(tidyr)})
ROOT <- PGV("replication_work")
mimic <- read_csv(file.path(ROOT,"02_mimic_harmonized/mimic_model_data.csv"), show_col_types=FALSE)
eicu  <- read_csv(file.path(ROOT,"03_eicu_harmonized/eicu_model_data.csv"), show_col_types=FALSE)
res   <- read_csv(file.path(ROOT,"05_models/all_model_results.csv"), show_col_types=FALSE)
att   <- read_csv(file.path(ROOT,"05_models/attenuation.csv"), show_col_types=FALSE)
med <- function(x) sprintf("%.1f (%.1f\u2013%.1f)", median(x,na.rm=TRUE), quantile(x,.25,na.rm=TRUE), quantile(x,.75,na.rm=TRUE))
pct <- function(x,n) sprintf("%d (%.1f%%)", sum(x,na.rm=TRUE), mean(x,na.rm=TRUE)*100)

# ---------- Main Table A ----------
mm <- mimic %>% filter(in_main, in_landmark, glucose_n>=2)
ee <- eicu %>% filter(in_main, in_landmark, glucose_n>=2)
tabA <- data.frame(
  characteristic = c("N (landmark, \u22652 glucose)","Hospitals, n","Age, years, median (IQR)","Male sex, n (%)",
    "Procedure: CABG, n (%)","Procedure: open valve, n (%)","Procedure: combined CABG+valve, n (%)",
    "Diabetes, n (%)","Severity score, median (IQR)","Creatinine (0\u201324h first), mg/dL",
    "Glucose measurements, n, median (IQR)","Mean glucose, mg/dL","GV (SD), mg/dL","POCT fraction",
    "Post-landmark hospital mortality, n (%)","ICU mortality, n (%)"),
  MIMIC_IV = c(nrow(mm), 1L, med(mm$age), pct(mm$sex=="M",nrow(mm)),
    pct(mm$procedure_category=="cabg",nrow(mm)), pct(mm$procedure_category=="valve",nrow(mm)), pct(mm$procedure_category=="combined",nrow(mm)),
    pct(mm$diabetes==1,nrow(mm)), paste0("SOFA ",med(mm$sofa)," / Charlson ",med(mm$charlson)), med(mm$creatinine),
    med(mm$glucose_n), med(mm$glucose_mean_24h), med(mm$glucose_sd_24h), sprintf("%.2f",mean(mm$poct_fraction,na.rm=TRUE)),
    pct(mm$post_landmark_hosp_mortality,nrow(mm)), pct(mm$icu_mortality,nrow(mm))),
  eICU_CRD = c(nrow(ee), length(unique(ee$hospitalid)), med(ee$age), pct(ee$sex=="Male",nrow(ee)),
    pct(ee$procedure_category=="cabg",nrow(ee)), pct(ee$procedure_category=="valve",nrow(ee)), pct(ee$procedure_category=="combined",nrow(ee)),
    pct(ee$diabetes,nrow(ee)), paste0("APACHE ",med(ee$apachescore)), med(ee$creatinine),
    med(ee$glucose_n), med(ee$glucose_mean_24h), med(ee$glucose_sd_24h), sprintf("%.2f",mean(ee$poct_fraction,na.rm=TRUE)),
    pct(ee$post_landmark_hosp_mortality,nrow(ee)), pct(ee$icu_mortality,nrow(ee))))
write_csv(tabA, file.path(ROOT,"06_tables/main_table_A_cohort_characteristics.csv"))

# ---------- Main Table B ----------
prim <- res %>% filter(outcome=="post_landmark_hosp_mortality", exposure=="gv_z", cohort=="main", model %in% c("M0","M1","M2","M3"))
fmt <- function(r) sprintf("%.3f (%.3f\u2013%.3f)", r$or, r$ci_lo, r$ci_hi)
att_prim <- att %>% filter(cohort=="main", outcome=="post_landmark_hosp_mortality", exposure=="gv_z", (is.na(extra) | extra==""))
tabB <- prim %>% mutate(estimate=mapply(function(o,l,h) sprintf("%.3f (%.3f\u2013%.3f)",o,l,h), or, ci_lo, ci_hi),
                        p_fmt=ifelse(p<0.001,"<0.001",sprintf("%.3f",p))) %>%
  select(db, cohort, outcome, model, n, events, estimate, p_fmt)
att_row <- data.frame(db=c("MIMIC-IV","eICU-CRD"), cohort="main", outcome="post_landmark_hosp_mortality",
                      model="attenuation_pct", n=NA, events=NA,
                      estimate=sprintf("%.1f%%", att_prim$attenuation_pct), p_fmt="")
tabB <- bind_rows(tabB, att_row)
ctx <- data.frame(db="MIMIC-IV (contextual, Cox)",
  cohort=c("SHR-GV","SHR-GV","open-core","open-core","GV-only","GV-only"),
  outcome=c("30-day","1-year","30-day","1-year","30-day","1-year"),
  model="frozen primary Cox", n=c(4869,4381,4757,4162,9787,9787), events=c(114,279,99,254,256,650),
  estimate=c("HR 1.13 (1.02\u20131.26)","HR 1.17 (1.07\u20131.28)","HR 1.34 (1.17\u20131.54)","HR 1.21 (1.09\u20131.34)","HR 1.09 (1.02\u20131.17)","HR 1.10 (1.04\u20131.16)"),
  p_fmt=c("0.022","<0.001","<0.001","<0.001","0.016","<0.001"))
tabB_full <- bind_rows(tabB, ctx)
write_csv(tabB_full, file.path(ROOT,"06_tables/main_table_B_cross_database_estimates.csv"))
write_csv(res, file.path(ROOT,"06_tables/s14_all_models_full.csv"))
write_csv(att, file.path(ROOT,"06_tables/s14b_attenuation_full.csv"))

# ---------- Supplementary tables ----------
# S11 variable mapping
s11 <- read.csv(textConnection(
"concept,MIMIC-IV,eICU-CRD,comparability
Age,age_at_admission (continuous),age ('> 89' censored; set to 90 + flag),direct
Sex,gender,gender,direct
Procedure CABG,ICD procedure codes (frozen codebook),APACHE operative dx + text cross-validation,high but different pathway
Procedure open valve,ICD procedure codes,APACHE operative dx + text cross-validation,high but different pathway
Procedure combined,cabg_flag & open_valve_flag,CABG+OPEN_VALVE categories,direct concept
Diabetes,diagnoses_icd derived,pasthistory paths + apachepredvar,concept-direct
Creatinine 0-24h,labevents first (itemid 50912/52024/52546),lab first ('creatinine'),direct
Severity,Charlson + first-day SOFA,APACHE IV/IVa score,substitute - NOT equivalent
BMI,omr/derived (94.5% complete),admission height/weight (97.4% complete),available; excluded by frozen retention rule
Glucose,chartevents+labevents mg/dL items,lab bedside+chemistry mg/dL,direct; different measurement mix
HbA1c/SHR,47.3% available,essentially absent (9 rows database-wide),NOT comparable
Hospital mortality,hospital_expire_flag + deathtime,hospitaldischargestatus,direct
ICU mortality,deathtime aligned to icustays.outtime,unitdischargestatus,approximate
30-day/1-year mortality,available,NOT available,NOT comparable in eICU
Exposure window,ICU intime + 0-1440 min,ICU admission offset 0-1440 min,direct (harmonized)"
), colClasses="character")
write_csv(s11, file.path(ROOT,"06_tables/s11_variable_mapping.csv"))

# S12 cohort flow
s12 <- data.frame(
 step=c("候选手术 stay","成人且首个 stay","CABG/open valve 主队列","24h landmark 内","\u22652 血糖时间点","主分析 complete-case (Model 2)"),
 MIMIC_IV=c("12,992 (冻结全集)","12,992","9,707","9,690","9,611","9,604"),
 eICU_CRD=c("13,228 (任意证据)","11,840","9,050","9,005","8,460","7,115"))
write_csv(s12, file.path(ROOT,"06_tables/s12_cohort_flow.csv"))

# S13 eICU baseline (已在 Table A 覆盖大部分; 补充 sensC/broad)
write_csv(data.frame(
  cohort=c("main (DEFINITE)","sensC broad (DEFINITE+PROBABLE)"),
  n=c(sum(ee$in_main), sum(eicu %>% filter(in_sensc, in_landmark, glucose_n>=2) %>% nrow())),
  deaths=c(sum(ee$post_landmark_hosp_mortality), eicu %>% filter(in_sensc, in_landmark, glucose_n>=2) %>% summarise(s=sum(post_landmark_hosp_mortality)) %>% pull(s)),
  mort_pct=round(c(mean(ee$post_landmark_hosp_mortality)*100, eicu %>% filter(in_sensc, in_landmark, glucose_n>=2) %>% summarise(m=mean(post_landmark_hosp_mortality)*100) %>% pull(m)),2)),
 file.path(ROOT,"06_tables/s13_eicu_cohorts.csv"))

# S16 hospital distribution
write_csv(read_csv(file.path(ROOT,"04_quality_control/qc_eicu_hospital_distribution.csv"), show_col_types=FALSE),
          file.path(ROOT,"06_tables/s16_hospital_distribution.csv"))

# S17 glucose QC
s17 <- data.frame(
  metric=c("raw glucose rows","exact duplicates removed","minute time points","multi-record minutes (median collapsed)",
           "implausible range","unit check","POCT fraction (main)","central lab fraction","blood gas fraction","ICU charted fraction"),
  MIMIC_IV=c("264,509","28","183,540","79,601 (43.4%)","excluded 20\u20131500 mg/dL outside","mg/dL itemid-restricted",
             sprintf("%.3f",mean(mm$poct_fraction,na.rm=TRUE)), sprintf("%.3f",mean(mm$central_lab_fraction,na.rm=TRUE)),
             sprintf("%.3f",mean(mm$blood_gas_fraction,na.rm=TRUE)), sprintf("%.3f",mean(mm$icu_charted_fraction,na.rm=TRUE))),
  eICU_CRD=c("160,418+","115 (含不合理值)","162,388","1,326 (0.82%)","20\u20131500 mg/dL","mg/dL 全部",
             sprintf("%.3f",mean(ee$poct_fraction,na.rm=TRUE)), sprintf("%.3f",mean(ee$central_lab_fraction,na.rm=TRUE)),"0","0"))
write_csv(s17, file.path(ROOT,"06_tables/s17_glucose_qc.csv"))

# S18 landmark exclusions
s18 <- data.frame(
  item=c("队列 n","24h 内死亡","24h 内出院","landmark 内 n","landmark 前死亡率 %","post-landmark 死亡率 %"),
  MIMIC_IV=c(9707,10,7,9690,"58.8","1.82"),
  eICU_CRD=c(9050,20,25,9005,"44.4","1.73"))
write_csv(s18, file.path(ROOT,"06_tables/s18_landmark.csv"))

# S19 full coefficients
write_csv(read_csv(file.path(ROOT,"05_models/primary_full_coefficients.csv"), show_col_types=FALSE),
          file.path(ROOT,"06_tables/s19_full_coefficients.csv"))
# S15 sensitivities = res 中敏感性行
s15 <- res %>% filter(!(cohort=="main" & outcome=="post_landmark_hosp_mortality" & exposure=="gv_z" & (is.na(extra) | extra=="") & model %in% c("M0","M1","M2","M3")))
write_csv(s15, file.path(ROOT,"06_tables/s15_sensitivity_analyses.csv"))
# standardization params
write_csv(read_csv(file.path(ROOT,"05_models/standardization_params.csv"), show_col_types=FALSE),
          file.path(ROOT,"06_tables/s_standardization_params.csv"))

# ---------- 主图: cross-database forest plot ----------
fp <- res %>% filter(outcome=="post_landmark_hosp_mortality", exposure=="gv_z", cohort=="main", model %in% c("M2","M3"), (is.na(extra) | extra=="")) %>%
  mutate(panel=paste0(db, " — harmonized hospital mortality (logistic OR per 1-SD GV)"),
         label=sprintf("OR %.2f (%.2f\u2013%.2f)", or, ci_lo, ci_hi),
         model=factor(model, levels=c("M3","M2"), labels=c("Model 3: + mean glucose","Model 2: severity-adjusted")))
ctx_df <- data.frame(
  model=factor(c("30-day mortality (Cox)","1-year mortality (Cox)","30-day mortality (Cox)","1-year mortality (Cox)")),
  or=c(1.13,1.17,1.34,1.21), ci_lo=c(1.02,1.07,1.17,1.09), ci_hi=c(1.26,1.28,1.54,1.34),
  panel="MIMIC-IV frozen primary (contextual Cox HR, date-anchored)",
  label=c("HR 1.13 (1.02\u20131.26)","HR 1.17 (1.07\u20131.28)","HR 1.34 (1.17\u20131.54)","HR 1.21 (1.09\u20131.34)"),
  cohort=c("SHR\u2013GV","SHR\u2013GV","open-core","open-core"))
plot_df <- bind_rows(fp %>% select(panel,model,or,ci_lo,ci_hi,label) %>% mutate(cohort="harmonized"), ctx_df)
ctx_lab <- plot_df %>% filter(cohort!="harmonized") %>% group_by(panel, model) %>%
  mutate(hj=ifelse(or==min(or), 1.05, -0.05),
         lab_x=ifelse(or==min(or), ci_lo, ci_hi)) %>% ungroup()
p <- ggplot(plot_df, aes(x=or, y=model, color=grepl("contextual",panel))) +
  geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.2) +
  geom_point(size=2.5) +
  geom_text(data=plot_df %>% filter(cohort=="harmonized"), aes(label=label), nudge_y=0.28, size=3) +
  geom_text(data=ctx_lab, aes(label=label, hjust=hj, x=lab_x), size=2.8) +
  facet_wrap(~panel, ncol=1, scales="free_y") +
  scale_x_log10(limits=c(0.8,1.8)) +
  scale_color_manual(values=c("FALSE"="black","TRUE"="grey55"), guide="none") +
  labs(x="Effect estimate per 1-SD increase in GV (log scale; OR and HR shown separately, no pooled estimate)",
       y=NULL,
       caption="ORs from database-specific logistic models (eICU: hospital-clustered robust SE; MIMIC: robust SE).\nContextual HRs from the frozen MIMIC-IV primary Cox analyses (different outcome horizons and windows; not pooled).") +
  theme_bw(base_size=11)
ggsave(file.path(ROOT,"07_figures/cross_database_forest_plot.png"), p, width=9, height=7, dpi=600)
ggsave(file.path(ROOT,"07_figures/cross_database_forest_plot.pdf"), p, width=9, height=7)
write_csv(plot_df, file.path(ROOT,"07_figures/figure_source_data/forest_plot_data.csv"))

# ---------- 补充图1: eICU cohort flow ----------
flow <- data.frame(step=factor(1:6, labels=c("Any cardiac-surgery evidence\n(n=13,228)","Broad DEFINITE+PROBABLE\nfirst stay (n=11,840)",
  "CABG/open valve DEFINITE\n(n=9,050)","24h landmark\n(n=9,005)","\u22652 glucose time points\n(n=8,460)","Model-2 complete cases\n(n=7,115)")),
  n=c(13228,11840,9050,9005,8460,7115))
pf <- ggplot(flow, aes(x=step, y=n)) + geom_col(fill="steelblue") + geom_text(aes(label=scales::comma(n)), vjust=-0.3, size=3.5) +
  scale_y_continuous(labels=scales::comma, limits=c(0,14000)) + labs(x=NULL, y="Patients", title="eICU-CRD harmonized replication cohort flow") + theme_bw(base_size=11)
ggsave(file.path(ROOT,"07_figures/eicu_cohort_flow.png"), pf, width=9, height=5, dpi=600)
ggsave(file.path(ROOT,"07_figures/eicu_cohort_flow.pdf"), pf, width=9, height=5)
write_csv(flow, file.path(ROOT,"07_figures/figure_source_data/eicu_flow_data.csv"))

# ---------- 补充图2: glucose measurement count 分布 ----------
cnt <- bind_rows(mm %>% select(glucose_n) %>% mutate(db="MIMIC-IV"), ee %>% select(glucose_n) %>% mutate(db="eICU-CRD"))
pc <- ggplot(cnt, aes(x=glucose_n, fill=db)) + geom_histogram(binwidth=2, position="identity", alpha=0.5) +
  labs(x="Glucose time points in ICU 0\u201324h (after per-minute median collapse)", y="Patients", fill=NULL) + theme_bw(base_size=11) +
  facet_wrap(~db, ncol=1)
ggsave(file.path(ROOT,"07_figures/sfig_glucose_count_distribution.png"), pc, width=8, height=6, dpi=600)
ggsave(file.path(ROOT,"07_figures/sfig_glucose_count_distribution.pdf"), pc, width=8, height=6)

# ---------- 补充图3: GV ~ 测量次数 ----------
d2 <- bind_rows(mm %>% select(glucose_n, glucose_sd_24h) %>% mutate(db="MIMIC-IV"),
                ee %>% select(glucose_n, glucose_sd_24h) %>% mutate(db="eICU-CRD"))
pg <- ggplot(d2, aes(x=glucose_n, y=glucose_sd_24h, color=db)) + geom_point(alpha=0.08, size=0.6) +
  geom_smooth(method="loess", se=FALSE) + facet_wrap(~db) +
  labs(x="Glucose time points (0\u201324h)", y="GV (SD, mg/dL)", color=NULL) + theme_bw(base_size=11)
ggsave(file.path(ROOT,"07_figures/sfig_gv_vs_count.png"), pg, width=8, height=5, dpi=600)
ggsave(file.path(ROOT,"07_figures/sfig_gv_vs_count.pdf"), pg, width=8, height=5)

# ---------- 补充图4: eICU 医院分布 ----------
hd <- read_csv(file.path(ROOT,"04_quality_control/qc_eicu_hospital_distribution.csv"), show_col_types=FALSE)
ph <- ggplot(hd %>% arrange(desc(cases)) %>% mutate(rk=1:n()), aes(x=rk, y=cases)) + geom_col(fill="steelblue") +
  labs(x=paste0("Hospitals (ranked, n=",nrow(hd),")"), y="Cases (main cohort)") + theme_bw(base_size=11)
ggsave(file.path(ROOT,"07_figures/sfig_hospital_distribution.png"), ph, width=8, height=5, dpi=600)
ggsave(file.path(ROOT,"07_figures/sfig_hospital_distribution.pdf"), ph, width=8, height=5)

# ---------- 补充图5: 敏感性森林图 ----------
s15p <- res %>% filter(model %in% c("M2","M3"), exposure=="gv_z", (is.na(extra) | extra==""), grepl("random",model)==FALSE) %>%
  mutate(grp=paste0(db," | ",cohort," | ",outcome), estimate=sprintf("%.2f (%.2f\u2013%.2f)",or,ci_lo,ci_hi))
ps <- ggplot(s15p, aes(x=or, y=paste(grp,model), color=db)) + geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.2) + geom_point(size=2) + scale_x_log10() +
  labs(x="OR per 1-SD GV (log scale)", y=NULL, color=NULL, title="Prespecified sensitivity analyses (Model 2 / Model 3)") + theme_bw(base_size=10)
ggsave(file.path(ROOT,"07_figures/sfig_sensitivity_forest.png"), ps, width=9, height=10, dpi=600)
ggsave(file.path(ROOT,"07_figures/sfig_sensitivity_forest.pdf"), ps, width=9, height=10)

cat("TABLES_FIGURES_DONE\n")
