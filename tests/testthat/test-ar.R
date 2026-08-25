
# to test simulate a simple system
set.seed(1)

A <- matrix(0.9,1,1)
y <- sim_ar(A,Sigma=1,T=100,stationary =TRUE)
y_mat <- tsibble_to_mat(y)

test_that("mlag() function works", {


  # test that mlag works
  expect_snapshot(mlag(t(y_mat), 3))
  # test errors on wrong lag size
  expect_error(mlag(t(y_mat), dim(y_mat)[1]),"Maximum lag length must be smaller than number of data points")
  expect_error(mlag(t(y_mat), 0),"Minimum lag length must be 1")


})

test_that("fit_ar() function works", {


  # test that fit_ar works in auto mode
  expect_snapshot(fit_ar(y, k=-1, kmax=12))
  # explicit k
  expect_snapshot(fit_ar(y, k=3))

  # test errors on wrong max lag size
  expect_error(fit_ar(y, k=-1, kmax=dim(y_mat)[1]),"Maximum lag length must be smaller than number of data points")

  expect_error(fit_ar(y, k=-1, kmax=0),"Minimum lag length must be 1")

  # test warning on zero lag length
  expect_snapshot(fit_ar(y, k=0))


})
