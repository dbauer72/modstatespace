
test_that("modstatespace::plot_IR works", {
  #library(tidyverse)
  syst <- StSp(A=as.matrix(0.5), C=as.matrix(1), K=as.matrix(1), Sigma = as.matrix(1), type = "innov")
  # get some covariances ####
  expect_no_error(C1 <- plot_IR(syst))
  syst$A[1,1] <- -0.9
  expect_no_error(C2 <- plot_IR(syst))
  syst$A[1,1] <- 1
  expect_no_error(plot_IR(syst))
  syst$A[1,1] <- 0.999
  expect_no_error(C3 <- plot_IR(syst))


})
