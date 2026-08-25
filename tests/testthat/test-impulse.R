test_that("impulse function works", {

  A=as.matrix(0.5) 
  C=as.matrix(1) 
  K=as.matrix(1)
  syst <- StSp(A=A,C=C,K=K,Sigma=K,type = 'innov')  
  
#  impulse response for (0.5, 1, 1)
im1 <- modstatespace::impulse(syst, M=20)
#  impulse response for (1, 1, 1)
syst$A[1,1] = 1
im2 <- modstatespace::impulse(syst, M=20)
#  impulse response for (-0.6, 1, 1)
syst$A[1,1] = -0.6

im3 <- modstatespace::impulse(syst, M=20)

expect_snapshot(
  im1
)

expect_snapshot(
  im2
)

expect_snapshot(
  im3
)

})
