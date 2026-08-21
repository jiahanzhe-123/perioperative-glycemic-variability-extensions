# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 08_harmonized.R — harmonized MIMIC–eICU 短期院内死亡比较(非外部验证;不合并不池化)
# 主结局:24h landmark 后院内死亡。主模型:modified Poisson + 稳健方差(RR + 标准化风险差)。
# 主暴露尺度:per 10 mg/dL GV;数据库内 1-SD 仅作补充(两库 SD 不同,不直接比较)。
# 偏差登记:eICU 原始库本轮不可达,eICU 侧只能用既有全来源抽取;
#          MIMIC 侧主分析用 common-source(POCT+central lab)序列,eICU 侧无法对应限制。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(sandwich); library(lmtest); library(jsonlite)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/08_harmonized.log"), split=TRUE)

# ---- MIMIC:harmonized 队列 + common-source GV ----
mm <- read.csv(file.path(PGV("method_audit_work"),"00_audit","mimic_model_data_charlson_wo_diabetes.csv"), stringsAsFactors=FALSE)
mm <- mm[mm$in_main %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
mm <- mm[mm$in_landmark %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
com <- read.csv(file.path(ROOT,"data","features_source_common.csv"), stringsAsFactors=FALSE)
mm <- merge(mm, com[,c("stay_id","gv","mean_glu","n_src")], by="stay_id", all.x=TRUE)
mm$gv10_common <- mm$gv/10
mm <- mm[!is.na(mm$gv10_common) & !is.na(mm$post_landmark_hosp_mortality),]
mm$y <- as.integer(mm$post_landmark_hosp_mortality %in% c(TRUE,"t","True","TRUE","true",1,"1"))
mm$diabetes <- as.integer(mm$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
mm$sex <- factor(mm$sex); mm$procedure_category <- factor(mm$procedure_category)
mm$log_count <- log(mm$glucose_n)
cat("MIMIC harmonized(common-source): N =", nrow(mm), "; 事件 =", sum(mm$y), "\n")

# ---- eICU:既有 harmonized 抽取(全来源,偏差登记) ----
ee <- read.csv(file.path(PGV("replication_work"),"03_eicu_harmonized","eicu_model_data.csv"), stringsAsFactors=FALSE)
ee <- ee[ee$in_main %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
ee <- ee[ee$in_landmark %in% c(TRUE,"t","True","TRUE","true",1,"1"),]
ee$y <- as.integer(ee$hosp_mortality %in% c(TRUE,"t","True","TRUE","true",1,"1") &
                   !(ee$died_within_24h %in% c(TRUE,"t","True","TRUE","true",1,"1")))
ee$gv10 <- ee$glucose_sd_24h/10; ee$mean_glu <- ee$glucose_mean_24h
ee$diabetes <- as.integer(ee$diabetes %in% c(TRUE,"t","True","TRUE","true",1,"1"))
ee$sex <- factor(ee$sex); ee$procedure_category <- factor(ee$procedure_category)
ee$log_count <- log(ee$glucose_n)
ee <- ee[!is.na(ee$gv10) & !is.na(ee$y),]
cat("eICU harmonized(全来源抽取): N =", nrow(ee), "; 事件 =", sum(ee$y), "\n")

# ---- 模型层级(modified Poisson) ----
fit_mp <- function(dd, fml, cluster=NULL){
  f <- glm(as.formula(fml), data=dd, family=poisson(link="log"))
  V <- if(is.null(cluster)) vcovHC(f, type="HC1") else vcovCL(f, cluster=dd[[cluster]])
  list(fit=f, V=V)
}
grab <- function(obj, term){
  ct <- coeftest(obj$fit, vcov.=obj$V)
  list(RR=unname(exp(ct[term,"Estimate"])), lo=exp(ct[term,"Estimate"]-1.96*ct[term,"Std. Error"]),
       hi=exp(ct[term,"Estimate"]+1.96*ct[term,"Std. Error"]), P=ct[term,"Pr(>|z|)"])
}
risk_diff <- function(dd, fml, cluster=NULL){
  # 标准化风险差:GV 设为该库 Q75 vs Q25(边际平均)
  q <- quantile(dd$gv10*10, c(.25,.75))
  f <- fit_mp(dd, fml, cluster)
  d25 <- dd; d25$gv10 <- q[1]/10; d75 <- dd; d75$gv10 <- q[2]/10
  if ("gv10" %in% names(dd)) {
    p25 <- predict(f$fit, newdata=d25, type="response"); p75 <- predict(f$fit, newdata=d75, type="response")
    return(list(rd=mean(p75)-mean(p25), q25=q[1], q75=q[2]))
  }
  list(rd=NA, q25=q[1], q75=q[2])
}
rows <- list()
run_db <- function(dd, dbname, sev_fml, cluster, gvterm="gv10"){
  if (!"gv10" %in% names(dd)) dd$gv10 <- dd[[gvterm]]
  if (!"span_hours" %in% names(dd))
    dd$span_hours <- if("measurement_span_minutes" %in% names(dd)) dd$measurement_span_minutes/60 else NA_real_
  m1 <- "y ~ gv10 + age + sex + procedure_category + diabetes + creatinine"
  m2 <- paste0(m1, " + ", sev_fml)
  m3 <- paste0(m2, " + splines::ns(mean_glu, df=3)")
  m4 <- paste0(m3, " + log_count + span_hours")
  for (mi in 1:4) {
    fml <- get(paste0("m",mi))
    vars <- all.vars(as.formula(fml))
    dd2 <- dd[complete.cases(dd[, vars]),]
    if (nrow(dd2) < 100 || sum(dd2$y) < 10) {
      rows[[length(rows)+1]] <<- data.frame(model_id=paste0("HARM_",dbname,"_M",mi), database=dbname,
        note=paste0("not estimable (N=", nrow(dd2), ", events=", sum(dd2$y), ")"), stringsAsFactors=FALSE)
      next
    }
    obj <- fit_mp(dd2, fml, cluster)
    g <- grab(obj, "gv10")
    rd <- risk_diff(dd2, fml, cluster)
    sd_db <- sd(dd2$gv10*10, na.rm=TRUE)
    rows[[length(rows)+1]] <<- data.frame(model_id=paste0("HARM_",dbname,"_M",mi), database=dbname,
      outcome="post-landmark hospital mortality", model=paste0("modified Poisson M",mi),
      N=nrow(dd2), events=sum(dd2$y),
      RR_per10=g$RR, lo=g$lo, hi=g$hi, P=g$P,
      RR_perSD=exp(log(g$RR)*sd_db/10), lo_perSD=exp(log(g$lo)*sd_db/10), hi_perSD=exp(log(g$hi)*sd_db/10),
      sd_within_db=sd_db, rd_q75_q25=rd$rd, gv_q25=rd$q25, gv_q75=rd$q75,
      variance=if(is.null(cluster)) "HC1 robust" else paste0("cluster-robust (", cluster, ")"),
      stringsAsFactors=FALSE)
  }
  invisible(NULL)
}
run_db(mm, "MIMIC-IV", "charlson_without_diabetes + sofa", NULL, gvterm="gv10_common")
run_db(ee, "eICU-CRD", "apachescore", "hospitalid")
tab <- do.call(rbind, rows)
write.csv(tab, file.path(ROOT,"results","10_cross_database_results.csv"), row.names=FALSE)
print(tab[,c("model_id","database","model","N","events","RR_per10","lo","hi","P","rd_q75_q25","variance")])

# eICU hospital random-intercept logistic(敏感性)
ri <- tryCatch({
  suppressMessages(library(lme4))
  f <- glmer(y ~ gv10 + age + sex + procedure_category + diabetes + creatinine + apachescore + (1|hospitalid),
             data=ee, family=binomial, control=glmerControl(optimizer="bobyqa"))
  s <- summary(f)$coefficients["gv10",]
  data.frame(model_id="HARM_eICU_randinteract", database="eICU-CRD",
    outcome="post-landmark hospital mortality", model="hospital random-intercept logistic (sensitivity)",
    N=nrow(ee), events=sum(ee$y), OR_per10=exp(s["Estimate"]), lo=exp(s["Estimate"]-1.96*s["Std. Error"]),
    hi=exp(s["Estimate"]+1.96*s["Std. Error"]), P=s["Pr(>|z|)"], stringsAsFactors=FALSE)
}, error=function(e) data.frame(model_id="HARM_eICU_randinteract", note=conditionMessage(e)))
write.csv(ri, file.path(ROOT,"results","eicu_random_intercept_sensitivity.csv"), row.names=FALSE)
print(ri)
cat("PHASE8_DONE\n")
sink()
