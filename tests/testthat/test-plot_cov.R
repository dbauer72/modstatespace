
test_that("modstatespace::plot_cov works", {

  #library(tidyverse)
  syst_vec <- list()
  # get some covariances ####
  syst <- StSp(A=as.matrix(0.5), C=as.matrix(1),  K=as.matrix(1), Sigma = as.matrix(1), type= "innov")
  syst_vec[[1]] <- syst
  syst$A[1,1] <- 0.999
  syst_vec[[2]] <- syst
  syst$A[1,1] <- -0.6
  syst_vec[[3]] <- syst
  
  expect_no_error(C1 <- plot_cov(syst_vec))
})
