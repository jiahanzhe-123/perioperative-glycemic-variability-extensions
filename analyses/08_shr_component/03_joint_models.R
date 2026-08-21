# [pgv] 路径经配置驱动(src/common/paths.R);不得在此写入本机绝对路径。
if (!exists("PGV")) {
  .candidates <- c("src/common/paths.R", "../src/common/paths.R", "../../src/common/paths.R", "../../../src/common/paths.R")
  .hit <- .candidates[file.exists(.candidates)]
  if (length(.hit) == 0) stop("src/common/paths.R not found; run from repository root")
  source(.hit[1])
}
#!/usr/bin/env Rscript
# 21_joint_models.R — SHR–GV 联合模块:J0/J1/N1/N2/N3/D0–D3/I1 + D1 池化联合检验 + PH + timing + IPW
# Lead secondary test = 30d Model J1 联合整体检验(D1 pooled Wald);
# key secondary test = 365d 同一检验。单系数无论 P 完整报告。
rm(list=ls()); options(stringsAsFactors=FALSE, scipen=999)
SEED <- 20260726L; set.seed(SEED)
.pgv_src <- c("src/common/paths.R","../src/common/paths.R","../../src/common/paths.R","../../../src/common/paths.R")
.pgv_hit <- .pgv_src[file.exists(.pgv_src)]
if (length(.pgv_hit) == 0) stop("src/common/paths.R not found; run from repository root")
source(.pgv_hit[1])
.rlib <- tryCatch(PGV("rlib_extra"), error=function(e) NULL)
if (!is.null(.rlib) && dir.exists(.rlib)) .libPaths(c(.rlib, .libPaths()))
suppressMessages({library(survival); library(jsonlite); library(rms); library(mice)})
ROOT <- PGV("mimic_record_work")
sink(file.path(ROOT,"logs/21_joint_models.log"), split=TRUE)

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
K  <- fromJSON(file.path(ROOT,"results","standardization_constants.json"))
JK <- fromJSON(file.path(ROOT,"results","shr_gv_joint_constants.json"))
COVS <- c("age_at_admission","gender","bmi","diabetes","charlson_without_diabetes",
          "procedure_cat6","lactate_postop_first","creat_postop_first","sofa_24h")
covs_fml <- paste(COVS, collapse=" + ")
tgt <- d[d$landmark_eligible==TRUE & !is.na(d$gv) & d$glucose_count>=2,]
j  <- tgt[tgt$hba1c_tier=="1-90 days pre-index" & tgt$eAG>0,]
stopifnot(nrow(j)==JK$N_joint)

prep <- function(dd){
  dd$shr_zf <- (dd$shr - JK$shr_mean)/JK$shr_sd
  dd$gv_zf  <- (dd$gv  - JK$gv_mean )/JK$gv_sd
  dd$shr_c  <- dd$shr - JK$shr_center
  dd$gv_c   <- dd$gv  - JK$gv_center
  dd$mean_c <- dd$mean_glu - median(j$mean_glu)
  dd$a1c_c  <- dd$hba1c_pct - median(j$hba1c_pct)
  dd
}
j <- prep(j)
kt_shr  <- paste(format(JK$shr_knots3,  digits=10), collapse=",")
kt_gv   <- paste(format(JK$gv_knots3,   digits=10), collapse=",")
kt_mean3<- paste(format(JK$mean_knots3, digits=10), collapse=",")
kt_mean4<- paste(format(JK$mean_knots4, digits=10), collapse=",")
kt_a1c  <- paste(format(JK$a1c_knots3,  digits=10), collapse=",")

FML <- list(
  J0 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_fml)),
  J1 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml)),
  N1 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(shr, c(", kt_shr, ")) + gv_zf + ", covs_fml)),
  N2 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ shr_zf + rcs(gv, c(", kt_gv, ")) + ", covs_fml)),
  N3 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(shr, c(", kt_shr, ")) + rcs(gv, c(", kt_gv, ")) + ", covs_fml)),
  D0 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ ", covs_fml)),
  D1_30 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + ", covs_fml)),
  D1_365 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean4, ")) + rcs(hba1c_pct, c(", kt_a1c, ")) + ", covs_fml)),
  D1_365_lin = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean4, ")) + hba1c_pct + ", covs_fml)),
  D2_30 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + gv_zf + ", covs_fml)),
  D2_365 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean4, ")) + rcs(hba1c_pct, c(", kt_a1c, ")) + gv_zf + ", covs_fml)),
  D2_365_lin = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean4, ")) + hba1c_pct + gv_zf + ", covs_fml)),
  D3_30 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean3, ")) + hba1c_pct + gv_zf + mean_c:a1c_c + ", covs_fml)),
  D3_365 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ rcs(mean_glu, c(", kt_mean4, ")) + hba1c_pct + gv_zf + mean_c:a1c_c + ", covs_fml)),
  I1 = function(tm,ev) as.formula(paste0("Surv(",tm,",",ev,") ~ shr_c*gv_c + ", covs_fml))
)

# ---- MICE 子集 ----
imp_full <- readRDS(file.path(ROOT,"results","mice_m50_object.rds"))
jidx <- which(imp_full$data$stay_id %in% j$stay_id)
cat("MICE 数据行:", nrow(imp_full$data), "; 联合队列行:", length(jidx), "\n")
# MICE 数据不含 shr/HbA1c(不插补变量),各插补集内回联原始值
jkeys <- j[, c("stay_id","shr","hba1c_pct","eAG","hba1c_tier","hba1c_days_from_surgery")]

# ---- D1 pooled Wald(Li–Raghunathan–Rubin) ----
d1_wald <- function(ests, terms, label=""){
  ests <- ests[!sapply(ests, is.null)]
  ok <- sapply(ests, function(f) !is.null(coef(f)) && all(terms %in% names(coef(f))))
  if (!all(ok)) {
    badf <- ests[[which(!ok)[1]]]
    miss <- setdiff(terms, names(coef(badf)))
    return(list(stat=NA, df1=length(terms), df2=NA, p=NA, method="D1 pooled Wald (Li-Raghunathan-Rubin)",
                m=length(ests), note=paste0(label, ": terms [", paste(miss, collapse=","),
                "] missing in ", sum(!ok), "/", length(ests), " fits")))
  }
  Q <- t(sapply(ests, function(f) coef(f)[terms]))
  Q <- matrix(as.numeric(Q), ncol=length(terms), dimnames=list(NULL, terms))
  if (any(!is.finite(Q))) return(list(stat=NA, df1=length(terms), df2=NA, p=NA, note=paste0(label, ": non-finite coef")))
  U <- lapply(ests, function(f) vcov(f)[terms, terms, drop=FALSE])
  m <- nrow(Q); qbar <- colMeans(Q); ubar <- Reduce("+", U)/m
  B <- if (m>1) cov(Q) else matrix(0, length(terms), length(terms))
  Tm <- ubar + (1+1/m)*B
  if (qr(Tm)$rank < length(terms))
    return(list(stat=NA, df1=length(terms), df2=NA, p=NA, note=paste0(label, ": singular pooled vcov")))
  W <- as.numeric(t(qbar) %*% solve(Tm) %*% qbar)
  k <- length(terms); stat <- W/k
  r <- (1+1/m) * sum(diag(B %*% solve(Tm))) / k
  nu <- if (r>0) (m-1)*(1 + 1/r^2) else 1e9
  list(stat=stat, df1=k, df2=nu, p=pf(stat, k, nu, lower.tail=FALSE),
       method="D1 pooled Wald (Li-Raghunathan-Rubin)", m=m)
}
rubin1 <- function(ests, term){
  b <- sapply(ests, function(f) coef(f)[term])
  v <- sapply(ests, function(f) vcov(f)[term, term])
  m <- length(b); qb <- mean(b); U <- mean(v); B <- var(b); Tm <- U + (1+1/m)*B
  df <- (m-1)*(1+U/((1+1/m)*B))^2
  list(beta=qb, se=sqrt(Tm), df=df, p=2*pt(abs(qb/sqrt(Tm)), df=df, lower.tail=FALSE),
       fmi=(1+1/m)*B/Tm)
}

rows_lin <- list(); rows_omni <- list(); rows_nl <- list(); rows_dec <- list(); rows_int <- list()
add_lin  <- function(...) rows_lin [[length(rows_lin )+1]]  <<- data.frame(..., stringsAsFactors=FALSE)
add_omni <- function(...) rows_omni[[length(rows_omni)+1]]  <<- data.frame(..., stringsAsFactors=FALSE)
add_nl   <- function(...) rows_nl  [[length(rows_nl  )+1]]  <<- data.frame(..., stringsAsFactors=FALSE)
add_dec  <- function(...) rows_dec [[length(rows_dec )+1]]  <<- data.frame(..., stringsAsFactors=FALSE)
add_int  <- function(...) rows_int [[length(rows_int )+1]]  <<- data.frame(..., stringsAsFactors=FALSE)

per_scale <- function(rb, mult){
  c(HR=exp(rb$beta*mult), lo=exp((rb$beta - 1.96*rb$se)*mult), hi=exp((rb$beta + 1.96*rb$se)*mult))
}

run_horizon <- function(hz){
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  oc <- paste0("mortality by index day ",hz," among day-1 landmark survivors")
  ests <- list(J0=list(), J1=list(), N1=list(), N2=list(), N3=list(),
               D0=list(), D1=list(), D2=list(), D3=list(), I1=list())
  n_fail <- 0L
  for (i in 1:50) {
    dd <- merge(complete(imp_full, i)[jidx,], jkeys, by="stay_id", sort=FALSE)
    dd <- prep(dd)
    get <- function(key){
      kk <- key
      if (key=="D1") kk <- if(hz=="30") "D1_30" else "D1_365"
      if (key=="D2") kk <- if(hz=="30") "D2_30" else "D2_365"
      if (key=="D3") kk <- if(hz=="30") "D3_30" else "D3_365"
      f <- tryCatch(coxph(FML[[kk]](tm,ev), data=dd), error=function(e) NULL)
      if (is.null(f) && kk=="D1_365") f <- tryCatch(coxph(FML$D1_365_lin(tm,ev), data=dd), error=function(e) NULL)
      if (is.null(f) && kk=="D2_365") f <- tryCatch(coxph(FML$D2_365_lin(tm,ev), data=dd), error=function(e) NULL)
      if (is.null(f)) n_fail <<- n_fail + 1L
      f
    }
    for (key in names(ests)) ests[[key]][[i]] <- get(key)
  }
  cat(sprintf("[%s] 模型拟合失败(总计): %d\n", hz, n_fail))
  N <- length(jidx); events <- sum(j[[ev]])
  # ---- J1 单系数 + 联合检验 ----
  for (trm in c("shr_zf","gv_zf")) {
    rb <- rubin1(ests$J1, trm)
    if (trm=="shr_zf") {
      p01 <- list(beta=rb$beta/JK$shr_sd*0.1, se=rb$se/JK$shr_sd*0.1)
      add_lin( model_id=paste0("J1_",hz,"d_shr_per01"), analysis="lead secondary module (J1)",
        outcome=oc, term="SHR per 0.1", N=N, events=events,
        HR=exp(p01$beta), lo=exp(p01$beta-1.96*p01$se), hi=exp(p01$beta+1.96*p01$se),
        P=2*pnorm(abs(p01$beta/p01$se), lower.tail=FALSE))
    }
    mult <- if(trm=="gv_zf") 10/JK$gv_sd else 1
    sc <- per_scale(rb, mult)
    add_lin( model_id=paste0("J1_",hz,"d_",trm), analysis="lead secondary module (J1)",
      outcome=oc, term=if(trm=="shr_zf") "SHR per fixed-cohort SD" else "GV per 10 mg/dL", N=N, events=events,
      HR=sc["HR"], lo=sc["lo"], hi=sc["hi"], P=rb$p, fmi=rb$fmi)
    if (trm=="gv_zf")
      add_lin( model_id=paste0("J1_",hz,"d_gv_perSD"), analysis="lead secondary module (J1)",
        outcome=oc, term="GV per fixed-cohort SD", N=N, events=events,
        HR=exp(rb$beta), lo=exp(rb$beta-1.96*rb$se), hi=exp(rb$beta+1.96*rb$se), P=rb$p, fmi=rb$fmi)
  }
  w <- d1_wald(ests$J1, c("shr_zf","gv_zf"))
  add_omni( model_id=paste0("OMNI_J1_vs_J0_",hz,"d"), analysis=if(hz=="30") "LEAD SECONDARY TEST" else "KEY SECONDARY TEST",
    outcome=oc, comparison="J1 (SHR+GV+clinical) vs J0 (clinical)", N=N, events=events,
    stat=w$stat, df1=w$df1, df2=w$df2, P=w$p, method=w$method, n_imputations=w$m)
  # ---- N1/N2/N3 整体与非线性 ----
  shr_spline_terms <- grep("shr", names(coef(ests$N1[[1]])), value=TRUE, fixed=TRUE)
  shr_spline_terms <- setdiff(shr_spline_terms, "shr_zf")
  gv_spline_terms  <- grep("gv",  names(coef(ests$N2[[1]])), value=TRUE, fixed=TRUE)
  gv_spline_terms  <- setdiff(gv_spline_terms, "gv_zf")
  w1 <- d1_wald(ests$N1, c(shr_spline_terms, "gv_zf"))
  w1n <- d1_wald(ests$N1, setdiff(shr_spline_terms, "shr"))
  w2 <- d1_wald(ests$N2, c("shr_zf", gv_spline_terms))
  w2n <- d1_wald(ests$N2, setdiff(gv_spline_terms, "gv"))
  add_nl( model_id=paste0("N1_",hz,"d"), analysis="secondary (nonlinear)", outcome=oc,
    model="RCS(SHR,3k) + linear GV", N=N, events=events,
    overall_stat=w1$stat, overall_df1=w1$df1, overall_p=w1$p,
    nonlinear_stat=w1n$stat, nonlinear_df1=w1n$df1, nonlinear_p=w1n$p, method=w1$method, knots=kt_shr)
  add_nl( model_id=paste0("N2_",hz,"d"), analysis="secondary (nonlinear)", outcome=oc,
    model="linear SHR + RCS(GV,3k)", N=N, events=events,
    overall_stat=w2$stat, overall_df1=w2$df1, overall_p=w2$p,
    nonlinear_stat=w2n$stat, nonlinear_df1=w2n$df1, nonlinear_p=w2n$p, method=w2$method, knots=kt_gv)
  # N3 稳定性门槛
  n3_ok <- all(sapply(ests$N3, function(f) !is.null(f)))
  if (hz=="365" || hz=="30") {
    if (n3_ok) {
      n3terms <- setdiff(grep("shr|gv", names(coef(ests$N3[[1]])), value=TRUE), c("shr_zf","gv_zf"))
      w3 <- d1_wald(ests$N3, n3terms)
      add_nl( model_id=paste0("N3_",hz,"d"), analysis=if(hz=="30") "exploratory (technical appendix; 30d EPV 受限)" else "secondary (nonlinear, stability checked)",
        outcome=oc, model="RCS(SHR,3k) + RCS(GV,3k) additive", N=N, events=events,
        overall_stat=w3$stat, overall_df1=w3$df1, overall_p=w3$p, method=w3$method,
        note=paste0("all 50 imputations converged; EPV=", round(events/(length(coef(ests$N3[[1]]))+1),1)))
    } else add_nl( model_id=paste0("N3_",hz,"d"), analysis="not run", outcome=oc,
      model="RCS(SHR,3k) + RCS(GV,3k)", N=N, events=events,
      note="not all imputations converged; event count did not support stable double-spline model")
  }
  # ---- D 分解 ----
  wD <- d1_wald(ests$D2, "gv_zf")
  rbD <- rubin1(ests$D2, "gv_zf")
  scD <- per_scale(rbD, 10/JK$gv_sd)
  add_dec( model_id=paste0("D2_vs_D1_",hz,"d"), analysis="secondary (component decomposition)",
    outcome=oc, comparison="D2 (flexible mean + HbA1c + GV) vs D1 (flexible mean + HbA1c)",
    N=N, events=events, stat=wD$stat, df1=wD$df1, df2=wD$df2, P=wD$p, method=wD$method,
    gv_HR_per10=scD["HR"], gv_lo=scD["lo"], gv_hi=scD["hi"],
    note="GV added-information test beyond acute mean glucose and HbA1c")
  wD3 <- d1_wald(ests$D3, "mean_c:a1c_c")
  rbD3 <- rubin1(ests$D3, "mean_c:a1c_c")
  add_int( model_id=paste0("D3_interact_",hz,"d"), analysis="exploratory",
    outcome=oc, comparison="D3 = D2 + centered mean×HbA1c",
    N=N, events=events, stat=wD3$stat, df1=wD3$df1, df2=wD3$df2, P=wD3$p, method=wD3$method,
    interaction_HR=exp(rbD3$beta), int_lo=exp(rbD3$beta-1.96*rbD3$se), int_hi=exp(rbD3$beta+1.96*rbD3$se))
  # ---- I1 SHR×GV ----
  wI <- d1_wald(ests$I1, "shr_c:gv_c")
  rbI <- rubin1(ests$I1, "shr_c:gv_c")
  add_int( model_id=paste0("I1_SHRxGV_",hz,"d"), analysis="exploratory (interaction)",
    outcome=oc, comparison="I1 = centered SHR + centered GV + SHR×GV + clinical",
    N=N, events=events, stat=wI$stat, df1=wI$df1, df2=wI$df2, P=wI$p, method=wI$method,
    interaction_HR=exp(rbI$beta), int_lo=exp(rbI$beta-1.96*rbI$se), int_hi=exp(rbI$beta+1.96*rbI$se),
    note="per 1 SHR-unit × 1 mg/dL GV (centered)")
  invisible(NULL)
}
run_horizon("30"); run_horizon("365")

rbind_fill <- function(lst){
  allc <- unique(unlist(lapply(lst, names)))
  lst2 <- lapply(lst, function(x){ for (cn in setdiff(allc, names(x))) x[[cn]] <- NA; x[, allc, drop=FALSE] })
  do.call(rbind, lst2)
}
lin <- rbind_fill(rows_lin); omni <- rbind_fill(rows_omni)
nl <- rbind_fill(rows_nl); dec <- rbind_fill(rows_dec); int <- rbind_fill(rows_int)
write.csv(lin, file.path(ROOT,"results","shr_gv_joint_linear_results.csv"), row.names=FALSE)
write.csv(omni, file.path(ROOT,"results","shr_gv_joint_omnibus_tests.csv"), row.names=FALSE)
write.csv(nl, file.path(ROOT,"results","shr_gv_joint_nonlinearity.csv"), row.names=FALSE)
write.csv(dec, file.path(ROOT,"results","shr_gv_component_decomposition.csv"), row.names=FALSE)
write.csv(int, file.path(ROOT,"results","shr_gv_interaction_results.csv"), row.names=FALSE)

# ---- PH(CC 上 J1/N1/N2/D2/I1) ----
cc <- j[complete.cases(j[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","shr","gv","mean_glu",COVS)]) & j$t_lm_30>0,]
cat("joint CC N =", nrow(cc), "\n")
ph_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  for (mk in c("J1","N1","N2","D2_30","I1")) {
    kk <- mk
    if (mk=="D2_30" && hz=="365") kk <- "D2_365"
    f <- tryCatch(coxph(FML[[kk]](tm,ev), data=cc, x=TRUE, y=TRUE), error=function(e) NULL)
    if (is.null(f)) next
    z <- tryCatch(cox.zph(f, transform="km"), error=function(e) NULL)
    if (is.null(z)) {
      ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("PH_",kk,"_",hz,"d"),
        term="(all)", chisq=NA, df=NA, p=NA,
        note="cox.zph failed (singular system; model may be degenerate)", stringsAsFactors=FALSE)
      next
    }; tab <- as.data.frame(z$table); tab$term <- rownames(tab)
    key_terms <- grep("shr|gv|GLOBAL", tab$term, value=TRUE)
    for (i2 in which(tab$term %in% key_terms | tab$term=="GLOBAL"))
      ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("PH_",kk,"_",hz,"d"),
        term=tab$term[i2], chisq=tab$chisq[i2], df=tab$df[i2], p=tab$p[i2], stringsAsFactors=FALSE)
    # 违反处理:连续暴露 log(time) 交互
    pv <- tab$p; names(pv) <- tab$term
    viol <- intersect(names(pv)[pv<0.05], c("shr_zf","gv_zf","shr_c","gv_c","shr","gv","hba1c_pct"))
    if (length(viol)) {
      tt_fml <- paste0("Surv(",tm,",",ev,") ~ shr_zf + gv_zf + ", covs_fml, " + ",
                       paste(paste0("tt(", viol, ")"), collapse=" + "))
      f_tt <- tryCatch(coxph(as.formula(tt_fml), data=cc, tt=function(x,t,...) x*log(t)), error=function(e) NULL)
      ph_rows[[length(ph_rows)+1]] <- data.frame(model_id=paste0("PH_",kk,"_",hz,"d_tt"),
        term=paste(viol, collapse="+"), chisq=NA, df=NA, p=NA,
        note=paste0("log(time) interaction added for: ", paste(viol, collapse=","),
                    if(is.null(f_tt)) " (fit failed)" else ""), stringsAsFactors=FALSE)
    }
  }
}
allc <- unique(unlist(lapply(ph_rows, names)))
ph_rows2 <- lapply(ph_rows, function(x){ for (cn in setdiff(allc, names(x))) x[[cn]] <- NA; x[, allc, drop=FALSE] })
ph <- do.call(rbind, ph_rows2)
write.csv(ph, file.path(ROOT,"results","shr_gv_ph_diagnostics.csv"), row.names=FALSE)

# ---- HbA1c timing 敏感性(CC,J1 同型) ----
timing_rows <- list()
for (tier in list(c("1-90 days pre-index","91-365 days pre-index"),
                  c("1-90 days pre-index","91-365 days pre-index","index calendar day"),
                  c("91-365 days pre-index"))) {
  dd <- tgt[tgt$hba1c_tier %in% tier & tgt$eAG>0,]
  dd <- prep(dd)
  dd <- dd[complete.cases(dd[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","shr","gv","mean_glu",COVS)]) & dd$t_lm_30>0,]
  label <- paste(tier, collapse=" + ")
  if (nrow(dd) < 200) {
    timing_rows[[length(timing_rows)+1]] <- data.frame(tiers=label, note=paste0("not estimable (N=", nrow(dd), ")"))
    next
  }
  for (hz in c("30","365")) {
    tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
    f <- coxph(FML$J1(tm,ev), data=dd)
    s <- summary(f)
    for (trm in c("shr_zf","gv_zf"))
      timing_rows[[length(timing_rows)+1]] <- data.frame(tiers=label, model_id=paste0("TIMING_",hz,"d_",trm),
        N=nrow(dd), events=sum(dd[[ev]]), term=trm,
        HR=s$coefficients[trm,"exp(coef)"], lo=s$conf.int[trm,"lower .95"],
        hi=s$conf.int[trm,"upper .95"], P=s$coefficients[trm,"Pr(>|z|)"],
        note="complete-case; joint-cohort fixed constants retained", stringsAsFactors=FALSE)
  }
}
allc <- unique(unlist(lapply(timing_rows, names)))
timing_rows2 <- lapply(timing_rows, function(x){ for (cn in setdiff(allc, names(x))) x[[cn]] <- NA; x[, allc, drop=FALSE] })
timing <- do.call(rbind, timing_rows2)
write.csv(timing, file.path(ROOT,"results","shr_gv_hba1c_timing_sensitivity.csv"), row.names=FALSE)

# ---- IPW 敏感性(复用主项目 propensity + 99% 截断) ----
ipw_d <- tgt[!is.na(tgt$hba1c_tier),]
ipw_d$A <- as.integer(!is.na(ipw_d$hba1c_pct))
ipw_d <- ipw_d[complete.cases(ipw_d[,c("t_lm_30","event_lm_30","gv","mean_glu",COVS)]) & ipw_d$t_lm_30>0,]
ps_fit <- glm(A ~ age_at_admission + gender + bmi + diabetes + charlson_without_diabetes +
              procedure_cat6 + lactate_postop_first + creat_postop_first + sofa_24h, data=ipw_d, family=binomial)
ipw_d$p_hat <- predict(ps_fit, type="response")
ipw_d$w_stab <- mean(ipw_d$A==1)/ipw_d$p_hat
j_ipw <- prep(j[complete.cases(j[,c("t_lm_30","t_lm_365","event_lm_30","event_lm_365","shr","gv","mean_glu",COVS)]) & j$t_lm_30>0,])
j_ipw <- merge(j_ipw, ipw_d[,c("stay_id","w_stab","p_hat")], by="stay_id")
q99 <- unname(quantile(j_ipw$w_stab, .99))
j_ipw$w99 <- pmin(j_ipw$w_stab, q99)
ess <- sum(j_ipw$w99)^2/sum(j_ipw$w99^2)
ipw_rows <- list()
for (hz in c("30","365")) {
  tm <- paste0("t_lm_",hz); ev <- paste0("event_lm_",hz)
  f <- coxph(FML$J1(tm,ev), data=j_ipw, weights=w99, robust=TRUE)
  s <- summary(f)
  for (trm in c("shr_zf","gv_zf"))
    ipw_rows[[length(ipw_rows)+1]] <- data.frame(model_id=paste0("IPW_J1_",hz,"d_",trm),
      N=nrow(j_ipw), events=sum(j_ipw[[ev]]), weighted_events=sum(j_ipw$w99*j_ipw[[ev]]),
      ess=ess, ess_ratio=ess/nrow(j_ipw), max_weight=max(j_ipw$w99), truncation="99%",
      term=trm, HR=s$coefficients[trm,"exp(coef)"], lo=s$conf.int[trm,"lower .95"],
      hi=s$conf.int[trm,"upper .95"], P=s$coefficients[trm,"Pr(>|z|)"], stringsAsFactors=FALSE)
}
ipw <- do.call(rbind, ipw_rows)
prev <- read.csv(file.path(ROOT,"results","ipw_balance_smd.csv"), stringsAsFactors=FALSE)
max_smd <- max(abs(prev$smd_weighted_99), na.rm=TRUE)
flag <- max_smd > 0.10 || ess/nrow(j_ipw) < 0.5
ipw$selection_correction <- if(flag) "selection correction inadequate" else "acceptable balance (sensitivity only)"
write.csv(ipw, file.path(ROOT,"results","shr_gv_ipw_sensitivity.csv"), row.names=FALSE)
cat("IPW: max|SMD|(全队列 99%)=", round(max_smd,3), "; ESS/N(joint)=", round(ess/nrow(j_ipw),3),
    "; ", if(flag) "selection correction inadequate" else "ok", "\n")

cat("PHASE21_DONE\n")
sink()
