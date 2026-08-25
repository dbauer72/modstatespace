
test_that("modstatespace::covariance_function works", {

  syst <- StSp(A=as.matrix(0.5), C=as.matrix(1), K=as.matrix(1), Sigma = as.matrix(1), type = "innov")
  # get some covariances ####
  expect_no_error(C1 <- covariance_function.StSp(syst, M=20))
  syst$A[1,1] <- -0.9
  expect_no_error(C2 <- covariance_function.StSp(syst, M=20))
  syst$A[1,1] <- 1
  expect_warning(covariance_function.StSp(syst, M=20), 
               "Covariance_function: A has eigenvalues at z=1! Adjusting!",
               fixed = TRUE)
  syst$A[1,1] <- 0.999
  expect_no_error(C3 <- covariance_function.StSp(syst, M=20))

  # check covariances ####
  expect_snapshot(C1)
  expect_snapshot(C2)
  expect_snapshot(C3)

})
