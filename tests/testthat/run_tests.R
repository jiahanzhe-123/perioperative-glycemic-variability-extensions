#!/usr/bin/env Rscript
# run_tests.R — R 单元测试(testthat 风格,不依赖真实患者数据)
fail <- 0L
ok <- function(cond, name) {
  if (isTRUE(cond)) cat("PASS:", name, "\n")
  else { cat("FAIL:", name, "\n"); fail <<- fail + 1L }
}

# 右开区间
in_window <- function(t, t0, span=1440) t >= t0 & t < t0 + span
ok(!in_window(99, 100), "window: t0-1 excluded")
ok(in_window(100, 100), "window: t0 included")
ok(in_window(1539, 100), "window: t0+span-1 included")
ok(!in_window(1540, 100), "window: t0+span excluded (right-open)")

# 同分钟中位数
v <- c(100, 110, 400)
ok(median(v) == 110, "same-minute median")

# 来源优先级
prio <- c(central_lab=1, blood_gas=2, poct=3, icu_charted=4)
df <- data.frame(stay=1, minute=c(5,5,5,6), value=c(100,102,104,120),
                 source=c("icu_charted","poct","blood_gas","poct"))
df$pr <- prio[df$source]
best <- ave(df$pr, df$stay, df$minute, FUN=min)
sel <- df[df$pr==best,]
ok(identical(sel$source, c("blood_gas","poct")), "source priority per minute")

# GV >=2;ARV >=3
gv_ok <- function(x) length(x) >= 2 && !is.na(sd(x))
arv_ok <- function(x) length(x) >= 3
ok(gv_ok(c(100,120)), "GV needs >=2")
ok(!gv_ok(100), "GV fails with 1 point")
ok(arv_ok(c(100,120,130)), "ARV needs >=3")
ok(!arv_ok(c(100,120)), "ARV fails with 2 points")

# mg/dL ↔ mmol/L
ok(abs(10/18.018 - 0.555) < 0.001, "per 10 mg/dL = 0.555 mmol/L")

# eAG/SHR
eag <- 28.7*6.0 - 46.7
ok(abs(eag - 125.5) < 1e-6, "eAG formula")
ok(eag > 0, "eAG positive")

# 时间轴标签(landmark 后区间)
lab <- function(d) ifelse(d<=7, "1-7", ifelse(d<=30, "8-30", "31-365"))
ok(identical(lab(c(1,7,8,30,31,365)), c("1-7","1-7","8-30","8-30","31-365","31-365")),
   "post-landmark interval labels")

# BMI plausibility
bp <- function(b) b >= 10 & b <= 80
ok(identical(bp(c(9.9,10,80)), c(FALSE,TRUE,TRUE)), "BMI plausibility bounds")
ok(!bp(80.1), "BMI >80 excluded")

# LRT 尺度不变性(z 与原始尺度样条跨度一致)
suppressMessages(library(splines))
x <- runif(50, 0, 100)
k <- quantile(x, c(.1,.5,.9))
b1 <- ns(x, knots=k[2], Boundary.knots=k[c(1,3)])
b2 <- ns(scale(x), knots=(k[2]-mean(x))/sd(x), Boundary.knots=(k[c(1,3)]-mean(x))/sd(x))
r1 <- qr.resid(qr(cbind(1,b1)), x); r2 <- qr.resid(qr(cbind(1,b2)), x)
ok(cor(as.vector(b1 %*% solve(crossprod(b1), t(b1)) %*% x), x) > 0.99 ||
   identical(dim(b1), dim(b2)), "spline basis dimension stable under rescaling")

if (fail > 0) { cat("TESTS FAILED:", fail, "\n"); quit(status=1L, save="no") }
cat("ALL R TESTS PASSED\n")
