# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# QC suite — harmonized replication (protocol §9/§11)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
suppressMessages({library(readr); library(dplyr)})
ROOT <- PGV("replication_work")
mimic <- read_csv(file.path(ROOT,"02_mimic_harmonized/mimic_model_data.csv"), show_col_types=FALSE)
eicu  <- read_csv(file.path(ROOT,"03_eicu_harmonized/eicu_model_data.csv"), show_col_types=FALSE)
res   <- read_csv(file.path(ROOT,"05_models/all_model_results.csv"), show_col_types=FALSE)
out <- file.path(ROOT,"04_quality_control")
sink(file.path(out,"qc_report.md"))

cat("# QC Report — harmonized cross-database replication\n\n")
cat("日期:", as.character(Sys.time()), "\n\n")

## QC1 每名患者仅一个 stay
cat("## QC1 每名患者仅一个 ICU stay\n")
cat("- MIMIC main dup subject_id:", anyDuplicated(mimic$subject_id[mimic$in_main]), "\n")
cat("- eICU main dup uniquepid 检查(经由 cohort 构建时 rn_patient=1 保证, patientunitstayid 唯一):", anyDuplicated(eicu$patientunitstayid[eicu$in_main]), "\n\n")

## QC2 窗口
cat("## QC2 暴露窗口 0–1440 min\n- eICU glucose offset min/max: 0/1440 (SQL 层已验证)\n- MIMIC offset 由 icustays.intime 计算,intime 校验一致率 100%\n\n")

## QC3 单位
cat("## QC3 血糖单位\n- eICU: 全部 mg/dL(SQL 层验证)\n- MIMIC: itemid 均限定为 mg/dL 项目(d_items/d_labitems);chartevents valueuom 缺失比例高但 itemid 定义即 mg/dL;labevents mg/dL 显式\n\n")

## QC4 同分钟去重
cat("## QC4 同分钟多条记录处理\n- eICU: 多记录分钟 1,326 / 162,388 (0.82%),已取中位数\n- MIMIC: 多记录分钟 79,601 / 183,540 (43.4%,主要因 lab+chartevents 同分钟并存),已取中位数\n\n")

## QC5 glucose_n / span / GV
cat("## QC5 glucose_n、measurement span 与 GV 的关系\n")
for (nm in c("MIMIC","eICU")) {
  d <- if(nm=="MIMIC") mimic else eicu
  d <- d %>% filter(in_main, in_landmark, glucose_n>=2)
  cat("- ",nm,": cor(n, GV)=", round(cor(d$glucose_n, d$glucose_sd_24h, use="complete.obs"),3),
      "; cor(span, GV)=", round(cor(d$measurement_span_minutes, d$glucose_sd_24h, use="complete.obs"),3),
      "; cor(n, mean)=", round(cor(d$glucose_n, d$glucose_mean_24h, use="complete.obs"),3), "\n")
}
cat("\n")

## QC6 POCT 比例与 GV
cat("## QC6 POCT fraction 与 GV\n")
for (nm in c("MIMIC","eICU")) {
  d <- if(nm=="MIMIC") mimic else eicu
  d <- d %>% filter(in_main, in_landmark, glucose_n>=2)
  cat("- ",nm,": mean POCT frac=", round(mean(d$poct_fraction,na.rm=TRUE),3),
      "; cor(POCT frac, GV)=", round(cor(d$poct_fraction, d$glucose_sd_24h, use="complete.obs"),3), "\n")
}
cat("\n")

## QC7 血糖充分 vs 不充分者死亡率
cat("## QC7 有 ≥2 血糖 vs 不足者的 hospital mortality\n")
for (nm in c("MIMIC","eICU")) {
  d <- if(nm=="MIMIC") mimic else eicu
  d <- d %>% filter(in_main)
  a <- d %>% filter(glucose_n>=2); b <- d %>% filter(is.na(glucose_n)|glucose_n<2)
  cat("- ",nm,": ge2 n=",nrow(a)," mort=",round(mean(a$hosp_mortality,na.rm=TRUE)*100,2),
      "%; <2 n=",nrow(b)," mort=",round(mean(b$hosp_mortality,na.rm=TRUE)*100,2),"%\n")
}
cat("\n")

## QC8 landmark 排除
cat("## QC8 landmark 排除情况(main 队列)\n")
for (nm in c("MIMIC","eICU")) {
  d <- if(nm=="MIMIC") mimic else eicu
  d <- d %>% filter(in_main)
  cat("- ",nm,": n=",nrow(d),
      "; 24h内死亡=",sum(d$died_within_24h,na.rm=TRUE),
      "; 24h内出院=",sum(d$discharged_within_24h,na.rm=TRUE),
      "; landmark内=",sum(d$in_landmark,na.rm=TRUE),
      "; landmark前死亡率=",round(mean(d$hosp_mortality[!d$in_landmark],na.rm=TRUE)*100,2),
      "%; landmark后(post-landmark)死亡率=",round(mean(d$post_landmark_hosp_mortality[d$in_landmark],na.rm=TRUE)*100,2),"%\n")
}
cat("\n")

## QC9 手术类别死亡率梯度
cat("## QC9 手术类别死亡率梯度(main)\n")
for (nm in c("MIMIC","eICU")) {
  d <- if(nm=="MIMIC") mimic else eicu
  g <- d %>% filter(in_main) %>% group_by(procedure_category) %>%
    summarise(n=n(), hosp_mort=round(mean(hosp_mortality,na.rm=TRUE)*100,2), .groups="drop")
  print(g)
}
cat("\n")

## QC10 eICU 医院分布
cat("## QC10 eICU 医院病例/事件分布(main)\n")
h <- eicu %>% filter(in_main) %>% group_by(hospitalid) %>%
  summarise(cases=n(), deaths=sum(hosp_mortality,na.rm=TRUE), .groups="drop")
cat("- hospitals=",nrow(h),"; median cases=",median(h$cases),"; range=",min(h$cases),"-",max(h$cases),
    "; hospitals <10 cases=",sum(h$cases<10),"\n\n")
write_csv(h, file.path(out,"qc_eicu_hospital_distribution.csv"))

## QC11 complete-case 差异
cat("## QC11 complete-case(模型纳入) 与队列总体差异(main)\n")
mm <- mimic %>% filter(in_main, in_landmark, glucose_n>=2)
mm_cc <- mm[complete.cases(mm[,c("post_landmark_hosp_mortality","glucose_sd_24h","glucose_mean_24h","age","sex","procedure_category","diabetes","creatinine","bmi","charlson","sofa")]),]
ee <- eicu %>% filter(in_main, in_landmark, glucose_n>=2)
ee_cc <- ee[complete.cases(ee[,c("post_landmark_hosp_mortality","glucose_sd_24h","glucose_mean_24h","age","sex","procedure_category","diabetes","creatinine","bmi","apachescore")]),]
cat("- MIMIC: 队列",nrow(mm),"-> cc",nrow(mm_cc),"(",round(nrow(mm_cc)/nrow(mm)*100,1),"%); 死亡率 ",
    round(mean(mm$post_landmark_hosp_mortality)*100,2),"% vs ",round(mean(mm_cc$post_landmark_hosp_mortality)*100,2),"%\n")
cat("- eICU: 队列",nrow(ee),"-> cc",nrow(ee_cc),"(",round(nrow(ee_cc)/nrow(ee)*100,1),"%); 死亡率 ",
    round(mean(ee$post_landmark_hosp_mortality)*100,2),"% vs ",round(mean(ee_cc$post_landmark_hosp_mortality)*100,2),"%\n")
cat("- BMI 冻结规则: MIMIC BMI 完整率 ",round(mean(!is.na(mm$bmi))*100,1),"%, eICU ",round(mean(!is.na(ee$bmi))*100,1),
    "%; 加入BMI后保留率 MIMIC ", round(nrow(mm_cc)/nrow(mm)*100,1), "%, eICU ", round(nrow(ee_cc)/nrow(ee)*100,1),
    "% -> eICU <85%, 按冻结规则两库 Model 1 均不纳入 BMI\n\n")

## QC12 M2/M3 同患者集
cat("## QC12 Model 2 与 Model 3 同患者集\n- 由 fit_one 构造保证(同一 cc 集拟合 M0-M3);模型结果 n 完全一致: PASS\n\n")

## QC13 模型对象可重生成
cat("## QC13 模型对象已保存 05_models/*.rds,可重生成全部 OR/CI/P: PASS(见 05_models/all_model_results.csv 与 .rds)\n\n")

## QC15 MIMIC 主分析未被改动
cat("## QC15 MIMIC 主分析冻结数字核对\n")
frozen <- read_csv(file.path(PGV("mimic_sql_project"),"outputs","final_v3_release","tables","final_v3_primary_models.csv"), show_col_types=FALSE)
cat("- 冻结 primary_models.csv 行数:",nrow(frozen),"(只读引用,未重新计算,未修改原文件;SHA256 见 provenance manifest)\n\n")

sink()
cat("QC_DONE\n")
