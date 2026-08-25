#'
#' unobserved components model: Convert parameters to system matrices
#'
#' This function converts a parameter vector `theta` into system matrices for an unobserved components model
#' for a univariate quarterly time series with a doubly integrated trend, seasonal common cycles and a stationary
#' added univariate signal. 
#'
#' @param theta vector of parameters: C.. 7 params, A: 1 param, Sigma: 1 param. Rest fixed. 
#'
#' @return list of system matrices A, K, C, Sigma, Q
#'
#' @export
#'
param_syst_uc <- function(theta){
  
  # check parameters
  if(!is.vector(theta)){
    stop("theta must be a vector")
  }

  syst <- StSp()
  syst$type <- "UC"
  n <- 7
  syst$A <- matrix(0, nrow=n, ncol=n) # initialize A as zero matrix of size n x n
  
  # trend component
  syst$A[1,1] <- 1
  syst$A[1,2] <- 1
  syst$A[2,2] <- 1
  
  # common cycle 
  for (j in 1:3){
    syst$A[2+j, j+3] <- 1 # set A_{i,i+1} = 1
  }
  syst$A[6, 3] <- 1 # complete cycle 
  syst$A[7,7] <- theta[8]
  syst$C <- matrix(theta[1:7], nrow=1, ncol=7) # set C using parameters theta
  syst$Sigma <- matrix(theta[9]^2, nrow=1, ncol = 1)  # calculate Sigma as square of last parameter
  syst$Q <- diag(7) # set Q as diagonal matrix of size 7 x 7 (no state noise)
  syst$R <- syst$Sigma
  syst$S <- matrix(0,7,1)
  syst$type = "UC"
  
  return(syst)
}




#'
#' unobserved components model: Calculate criterion function / negative log likelihood
#'
#' This function calculates the criterion function, i.e. the negative log likelihood, for an unobserved components model given the parameter vector `theta` and observed data `y`.
#'
#' @param theta vector of parameters
#' @param y matrix of observations
#'
#' @return negative log likelihood
#'
#' @export
#'
cal_crit_uc <- function(theta, y){
  
  # check parameters
  if(!is.vector(theta)){
    stop("theta must be a vector")
  }
  if(!is.matrix(y)){
    stop("y must be a matrix")
  }
  
  s <- dim(y)[1] # get number of outputs
  T <- dim(y)[2] # number of observations 
  n <- 7
  syst <- param_syst_uc(theta) # get system matrices
  kf_out <- FKF::fkf(a0 = rep(0,n), P0 = diag(n), dt = matrix(0,n,T), ct =  matrix(0,s,T), Tt = syst$A, Zt = syst$C, HHt= syst$Q, GGt = syst$R, y= y)
  kf_out$LL <-  kf_out$logLik*(-2)
  return (kf_out$LL) # return negative log likelihood
}




#'
#' Estimate Unobserved Components Model
#'
#' This function estimates the parameters of an unobserved components model using a nonlinear least squares approach.
#'
#' @param y matrix of observations
#' @param theta_init vector of initial parameters
#' @param iterlim integer representing the maximum number of iterations for the nonlinear least squares estimation. Defaults to `100`.
#' @param initialize logical value indicating whether to use an initial estimate for the parameters. Defaults to `FALSE`.
#' @param print.level integer representing the level of `nlm` output to print. Defaults to `0`.
#'
#' @return list containing the estimated parameters.
#'
#'
#' @export
#'
est_uc <- function(y, theta_init, iterlim=100, initialize=FALSE, print.level=0){
  
  if (initialize == TRUE){
    syst <- param_syst_uc(theta_init) # get system matrices
    Sigma <- syst$Sigma # get innovation variance
    Q <- syst$Q # get state noise variance
    S <- syst$S
    kf_out <- kf(y=y, syst=syst, pe = FALSE, init=1) # calculate Kalman filter output
    # estimat H based on xh.
    xh <- kf_out$xhat # get estimated states
    Ch <- y %*% t(xh) %*% solve(xh %*% t(xh)) # calculate C using xh
    theta_init[1:7] <- matrix(Ch, ncol=1) # set H in theta_init
    # estimate Sig based on hh.
    Sig <- (y - Ch %*% xh) %*% t((y - Ch %*% xh))/dim(y)[2] # calculate innovation variance
    theta_init[9] <- sqrt(Sig) # extract parameters from innovation variance
  }
  theta_hat <- stats::nlm(cal_crit_uc, p=theta_init, y=y, print.level=print.level, iterlim=iterlim) # run nonlinear least squares estimation
  return(theta_hat)
}
