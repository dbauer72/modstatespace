

# we potentially got some optimization messages, stay international
OLDLANG <- Sys.getenv("LANGUAGE")
Sys.setenv(LANGUAGE = "en")
set.seed(1)
#library(fpp3)
y <- as.matrix(t(stats::arima.sim(n = 200, model = list(ar = c(0.5, -0.2, 0.3)), sd = 1)))
rownames(y) <- "y"

# Main functions

test_that("kf() works as expected", {
  syst <- StSp(A = as.matrix(0.5),C = as.matrix(1),K = as.matrix(1),Sigma = as.matrix(1))
  expect_snapshot(kf(y=y, syst=syst, pe=TRUE, init=0))
})

test_that("est_qml() works as expected", {
  expect_snapshot(est_qml(y, n=2, kmax=30, pe=TRUE, print.level=0))
})

test_that("cal_CVA() works as expected", {
  expect_snapshot(cal_CVA(y, f=5, p=5, n=1))
})


# Helper functions

test_that("conv_syst() works as expected", {
  syst <- StSp(A = as.matrix(0.5),C = as.matrix(1),K = as.matrix(1),Sigma = as.matrix(1))
  expect_snapshot(conv_syst(syst))
})

test_that("syst_param() works as expected", {
  syst <- StSp(A = as.matrix(0.5),C = as.matrix(1),K = as.matrix(1),Sigma = as.matrix(1))
  expect_snapshot(syst_param(syst))
})

test_that("param_syst() works as expected", {
  expect_snapshot(param_syst(theta = c(0.5, 1, 1), n = 1, s = 1))
})

test_that("param_syst_uc() works as expected", {
  expect_snapshot(param_syst_uc(c(0.5, 1, 1, 0.2, 3, 4, -3)))
})

test_that("cal_crit() works as expected", {
  expect_snapshot(cal_crit(c(0, 0, 1), n=1, s=1, y, pe=TRUE))
})

test_that("cal_crit_uc() works as expected", {
  # quite arbitrary values wich survive the function calls
  y <- matrix(0, nrow=1, ncol=10)
  y[1, ] <- c(-0.8885408, -0.5476704, -0.5264115, -0.6332790, -0.8153756, -0.8952983, -1.020254, -1.1674839, -0.7892407, -0.6995929)
  theta <- c(0.1,0.1,1,1,1,1,1,0.5,1)
  expect_snapshot(cal_crit_uc(theta, y))
})

test_that("cal_kf() works as expected", {
  expect_snapshot({kf_out <- cal_kf(c(0, 0, 1), n=1, s=1, y, pe=TRUE)
  kf_out$LL
  kf_out$Kt[[10]]})
})

test_that("cov_Sig_par() works as expected", {
  Sig <- matrix(c(2, 0.5, 0.5, 1), nrow=2, ncol=2)
  expect_snapshot(cov_Sig_par(Sig))
})

test_that("est_uc() works as expected", {
  # quite arbitrary values wich survive the function calls
  y <- matrix(0, nrow=1, ncol=10)
  y[1, ] <- c(-0.8885408, -0.5476704, -0.5264115, -0.6332790, -0.8153756, -0.8952983, -1.020254, -1.1674839, -0.7892407, -0.6995929)
  theta <- c(0.1,0.1,1,1,1,1,1,0.5,1)
  expect_snapshot(est_uc(y, theta, iterlim=60))
})


Sys.setenv(LANGUAGE = OLDLANG)
