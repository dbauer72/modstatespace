#'
#' Extract parameters from system matrices
#'
#' This function extracts the parameters from the system matrices (F, H, K) and returns them in vector form.
#'
#' @param syst a \code{StSp} object
#'
#' @return vector containing the extracted F, H, K parameters from the system matrices.
#'
#'
#' @export
#'
syst_param <- function(syst){

  if (!inherits(syst,"StSp")){
    stop("syst must be an StSp object!")
  }

  # check parameters
  if(!syst$type == "innov"){
    stop("syst must be in innovation representation.")
  }

  n <- dim(syst$A)[1] # Get the number of rows in matrix A
  s <- dim(syst$C)[1] # Get the number of rows in matrix C
  theta <- matrix(0, nrow=2*n*s+s*(s+1)/2, ncol=1) # Initialize theta as a zero matrix with dimensions 2*n*s by 1
  syst <- conv_syst(syst) # Convert the system into a canonical form using 'conv_syst'

  if (n >= s){ # Check if the state dimension is greater than or equal to the output dimension
    theta[1:(n*s)] <- matrix(syst$A[(n-s+1):n, ], ncol=1) # Populate theta with elements from syst$F
    theta[(n*s+1):(2*s*n)] <- matrix(syst$K, ncol=1) # Populate theta with elements from syst$K
  } else { # If state dimension is less than the output dimension
    # parameters also in C
    theta[1:((s-n)*n)] <- matrix(syst$C[(n+1):s, ], ncol=1) # Populate theta with elements from syst$H
    theta[((s-n)*n)+c(1:n^2)] <- matrix(syst$A, ncol=1) # Populate theta with elements from syst$F
    theta[(n*s+1):(2*s*n)] <- matrix(syst$K, ncol=1) # Populate theta with elements from syst$K
  }

  # parameters for Sigma 
  theta[2*n*s+(1:(s*(s+1)/2))] <- cov_Sig_par(syst$Sigma)
  
  # Return the matrix theta
  return(theta)
}


### other way round: For given parameter theta,
# obtain system (based on parameter vector (F, K, H)
# in canonical, innovation representation form
# For this we need to know the state dimension n and the
# output dimension s.
#
# Adds also the lower triangular part of the Cholesky factor
# for the innovation variance sigma to the parameter vector.

#'
#' Obtain system from parameter vector
#'
#' For given parameter theta, obtain system, based on parameter vector (F, K, H), in canonical form
#' For this we need to know the state dimension n and the output dimension s.
#'
#' Adds also the lower triangular part of the Cholesky factor for the innovation variance sigma to
#' the parameter vector.
#'
#' @param theta list. Parameter vector.
#' @param n integer. State dimension.
#' @param s integer. Output dimension.
#' @param gradient logical; indicates whether the supplied vector are the parameters or the standard deviations. 
#' @return syst as a \code{StSp} object
#'
#'
#' @export
#'
param_syst <- function(theta, n, s, gradient =FALSE){ # function definition with default parameters

  # check parameters
  if(!is.vector(theta)){
    stop("theta must be a vector")
  }
  if(n <= 0 || s <= 0){
    stop("n, s must be positive integers")
  }

  if ( is.logical(gradient) != TRUE ){
    stop("gradient must be logical")
  }
  
  fact <- 1
  if (gradient){
    fact <- 0 # for a gradient, fill everything with zeros. 
  }
  syst <- StSp() # initialize empty list to store system matrices
  syst$A <- matrix(0, nrow=n, ncol=n) # initialize A matrix as zero matrix of size n x n
  syst$C <- matrix(0, nrow=s, ncol=n) # initialize C matrix as zero matrix of size s x n
  syst$K <- matrix(0, nrow=n, ncol=s) # initialize K matrix as zero matrix of size n x s

  if (n >= s){ # check if state dimension is greater than or equal to output dimension
    if (n-s > 0){ # check if there is a difference between state and output dimensions
      syst$A[1:(n-s), (s+1):n] <- diag(n-s)*fact # fill part of F with identity matrix
    }
    syst$A[(n-s+1):n, ] <- matrix(theta[1:(n*s)], ncol=n) # fill the remaining part of C with elements from theta
    syst$C[, 1:s] <- diag(s)*fact # fill C with identity matrix
    syst$K <- matrix(theta[(n*s+1):(2*s*n)], ncol=s) # fill K with elements from theta
  } else { # if state dimension is less than output dimension
    syst$C[1:n, 1:n] <- diag(n)*fact # fill part of C with identity matrix
    syst$C[(n+1):s, ] <- matrix(theta[1:((s-n)*n)], ncol=n) # fill the remaining part of C with elements from theta
    syst$A <- matrix(theta[((s-n)*n)+c(1:(n^2))], nrow=n) # fill A with elements from theta
    syst$K <- matrix(theta[(n*s+1):(2*s*n)], nrow=n) # fill K with elements from theta
  }

  syst$Sigma <- matrix(0, s, s) # initialize Sigma as zero matrix of size s x s
  L <- matrix(0, s, s) # initialize L as zero matrix of size s x s
  nums <- s*(s+1)/2 # calculate the number of elements in the lower triangular part
  theta_L <- theta[2*s*n+c(1:nums)] # extract corresponding elements from theta
  for (j in 1:s){ # loop over the number of rows
    L[j, 1:j] <- theta_L[1:j] # fill L with elements from theta_L
    theta_L <- theta_L[-c(1:j)] # remove used elements from theta_L
  }
  syst$Sigma <- L %*% t(L) # calculate Sigma as L multiplied by its transpose
  syst$Sigma_chol <- L
  
  return(syst) # return the list containing system matrices
}


#'
#' Obtain derived system from parameter vector
#'
#' For given parameter theta, obtain the derived system, based on parameter vector (F, K, H), in canonical form
#' For this we need to know the state dimension n and the output dimension s. 
#'
#' Adds also the lower triangular part of the Cholesky factor for the innovation variance sigma to
#' the parameter vector.
#'
#' @param theta list. Parameter vector.
#' @param n integer. State dimension.
#' @param s integer. Output dimension.
#' 
#' @return systs. A list of \code{StSp} objects. First element is the system, the remaining the various derivatives.  
#'
#'
#' @export
#'
param_syst_deriv <- function(theta, n, s){ # function definition with default parameters
  
  # check parameters
  if(!is.vector(theta)){
    stop("theta must be a vector")
  }
  if(n <= 0 || s <= 0){
    stop("n, s must be positive integers")
  }
  
  syst <- param_syst(theta,n,s) # calculate the system corresponding to the parameter vector
  
  # set up calculations for Sigma derivatives
  L <- matrix(0, s, s) # initialize L as zero matrix of size s x s
  nums <- s*(s+1)/2 # calculate the number of elements in the lower triangular part
  theta_L <- theta[2*s*n+c(1:nums)] # extract corresponding elements from theta
  for (j in 1:s){ # loop over the number of rows
    L[j, 1:j] <- theta_L[1:j] # fill L with elements from theta_L
    theta_L <- theta_L[-c(1:j)] # remove used elements from theta_L
  }
  
  dsyst <- StSp() # initialize the derivative
  dsyst$A <- matrix(0, nrow=n, ncol=n) # initialize A matrix as zero matrix of size n x n
  dsyst$C <- matrix(0, nrow=s, ncol=n) # initialize C matrix as zero matrix of size s x n
  dsyst$K <- matrix(0, nrow=n, ncol=s) # initialize K matrix as zero matrix of size n x s
  dsyst$Sigma <- matrix(0,nrow=s,ncol=s) # initiale Sigma matrix as zero matrix of size n x s.
  
  # prepare output 
  systs <- list()
  systs[[1]] <- syst
  
  for (coord in 1:(2*n*s)){ # iterate over coordinates. 
    dthn <- theta*0 # dthn is the derivative of the parameter vector with respect to the coord coordinate
    dthn[coord] <- 1 
    
    dsystn <- dsyst # initialize empty
    if (n >= s){ # check if state dimension is greater than or equal to output dimension
      dsystn$A[(n-s+1):n, ] <- matrix(dthn[1:(n*s)], ncol=n) # fill the remaining part of A with elements from theta
      dsystn$K <- matrix(dthn[(n*s+1):(2*s*n)], ncol=s) # fill K with elements from theta
    } else { # if state dimension is less than output dimension
      dsystn$C[(n+1):s, ] <- matrix(dthn[1:(s-n)*n], ncol=n) # fill the remaining part of C with elements from theta
      dsystn$A <- matrix(dthn[((s-n)*n)+c(1:n^2)], nrow=n) # fill A with elements from theta
      dsystn$K <- matrix(dthn[(n*s+1):(2*s*n)], nrow=n) # fill K with elements from theta
    }
    
    systs <- append(systs,list(dsystn))
  } 
  
  # now derivative with respect to parameters for Sigma. 
  for (j in 1:nums){
    dsystn <- dsyst # initialize empty
    dL <- matrix(0, s, s) # initialize L as zero matrix of size s x s
    dthn_L <- theta[2*s*n+c(1:nums)]*0 # extract corresponding elements from theta
    dthn_L[j] <- 1
    for (j in 1:s){ # loop over the number of rows
      dL[j, 1:j] <- dthn_L[1:j] # fill L with elements from theta_L
      dthn_L <- dthn_L[-c(1:j)] # remove used elements from theta_L
    }
    dsystn$Sigma <- L %*% t(dL) + dL %*% t(L)
    systs <- append(systs,list(dsystn))
  }
  
  return(systs) # return the list containing system matrices
}


#'
#' Kalman filter
#'
#' This function calculates the Kalman filter for system (A, K, C) given as matrices.
#' General form using Sigma (innovation variance), Q (variance of state noise), S (cross covariance).
#'
#' @param y matrix of observed data where rows represent variables and columns represent time points
#' @param syst a \code{StSp} object.
#' @param pe logical. If TRUE the prediction error formulation is used. Default is `TRUE`.
#' @param init numeric. If =0 stationary distribution is assumed, if >0 initial state error variance
#'   equals init*I. Default is `0`.
#'
#' @return list xhat, res, sigmat: list of matrices Sigmat, Kt: list of Kalman gains, ll: -2 times
#'   log likelihood.
#'
#' @export
#'
#'
#' @export
#'
kf <- function(y, syst, pe=TRUE, init=0){

  # check parameters
  if (!(is.matrix(y) || !inherits(syst,"StSp"))) {
    stop("y must be a matrix, syst a StSp object")
  }
  if (!is.logical(pe)) {
    stop("pe must be a logical value")
  }
  if (!is.numeric(init)) {
    stop("init must be a numeric value")
  }

  n <- dim(syst$A)[1]
  s <- dim(syst$C)[1]
  T <- dim(y)[2] # assumes y has observations as columns
  ll <- 0
  
  A <- syst$A
  C <- syst$C
  
  if (syst$type == "innov"){
    Q <- syst$K %*% syst$Sigma %*% t(syst$K)
    S <- syst$K %*% syst$Sigma
    R <- syst$Sigma
    K <- syst$K
  } else {
    Q <- syst$Q
    S <- syst$S
    R <- syst$R
    K <- diag(n)
  }

  if (pe == TRUE){
    if (syst$type != "innov"){
      syst <- conv_to_innov(syst)
      K <- syst$K
    }
  }
  
  # check for non-invertibility, if prediction error for innovation representation is used.
  if (pe==TRUE){
    Abar <- A - K %*% C
    if (dim(Abar)[1]>1){
      ev <- max(abs(eigen(Abar)$values))
    } else {
      ev <- abs(Abar)
    }

    if (ev>0.99){
      # correct non-invertibility
      Abar <- Abar /ev*0.99
      ll <- ll+10^6*(exp(ev-0.99)-1)
    }
  }

  # initialize
  if (pe==TRUE){
    # initialize the filter
    eps <- y*0
    x <- matrix(0, n, T)
    eps[, 1] <- y[, 1]- C %*% x[, 1]
    Sigmat <- syst$Sigma
    iSigmat <- solve(Sigmat)
    Sigt <- list()
    Sigt[[1]] <- Sigmat
    Kt <- K
    Kts <- list()
    Kts[[1]] <- Kt

    Pt <- K %*% Sigmat %*% t(Kt)
    ll <- ll + log(det(Sigmat)) + t(eps[, 1]) %*% solve(Sigmat) %*% eps[, 1]

    # enter the recursion
    for (t in 2:T){
      x[, t] <- Abar %*% x[, t-1] + Kt %*% y[, t-1]
      eps[, t] <- y[, t] - C %*% x[, t]
      Pt <- Abar %*% Pt %*% t(Abar)

      ll <- ll + log(det(Sigmat)) + t(eps[, t]) %*% iSigmat %*% eps[, t]
    }
  } else {
    # calculate stationary state variance
    if (init == 0){
      #A <- diag(n^2) - kronecker(syst$F, syst$F)
      Om <- K %*% syst$Sigma %*% t(K)
      #vP <- solve(A) %*% matrix(Om, ncol=1)
      ev <- max(abs(eigen(A)$values))
      if (ev>0.99){
        # correct non-stability
        A <- A /ev*0.99
        ll <- ll+10^6*(exp(ev-0.99)-1)
      }
      
      Pt <- solve_lyap(A=A, Q= Om)
    } else {
      Pt <- diag(n)*init # init with multiples of identity 
    }
    #S <- K %*% syst$Sigma
    #Q <- Om
    # initialize the filter
    eps <- y*0
    x <- matrix(0, n, T)
    eps[, 1] <- y[, 1]
    Sigmat <- C %*% Pt %*% t(C) + R
    Sigt <- list()
    Sigt[[1]] <- Sigmat
    Kt <- (A %*% Pt %*% t(C) + S) %*% solve(Sigmat)
    Kts <- list()
    Kts[[1]] <- Kt

    ll <- ll + log(det(Sigmat)) + t(eps[, 1]) %*% solve(Sigmat) %*% eps[, 1]

    # enter the recursion
    for (t in 2:T){
      x[, t] <- A %*% x[, t-1] + Kt %*% eps[, t-1]
      Pt <- A %*% Pt %*% t(A) + Q - Kt %*% Sigmat %*% t(Kt)
      eps[, t] <- y[, t] - C %*% x[, t]
      Sigmat <- C %*% Pt %*% t(C) + R
      iSigmat <- solve(Sigmat)
      Kt <- (A %*% Pt %*% t(C) + S) %*% iSigmat

      ll <- ll + log(det(Sigmat)) + t(eps[, t]) %*% iSigmat %*% eps[, t]
      Sigt[[t]] <- Sigmat
      Kts[[t]] <- Kt
    }
  }

  # add the constant 
  ll = ll + T*s*log(2*pi)
  return (list(xhat=x, res=eps, sigmat=Sigt, Kt=Kts, Pt =  Pt, LL=ll))
}


#'
#' Calculate the criterion function as a function of the parameter vector
#'
#' This function calculates the criterion function with help of the Kalman filter as a function of the parameter vector.
#' It takes the parameter vector theta, the number of states n, the number of outputs s, the data y and a Boolean pe as input.
#' It returns the negative log likelihood of the model parameters.
#'
#' @param theta vector of parameters
#' @param n number of states
#' @param s number of outputs
#' @param y matrix of observations
#' @param pe Boolean whether to use prediction error formulation. Default is `TRUE`.
#'
#' @return negative log likelihood of the model
#'
#'
#' @export
#'
cal_crit <- function(theta, n, s, y, pe=TRUE){

  # check parameters
  if(!is.vector(theta)){
    stop("theta must be a vector")
  }
  if(n <= 0){
    stop(sprintf("n must be a positive integer %d", n))
  }
  if(s <= 0){
    stop("s must be a positive integer")
  }
  if(!is.matrix(y)){
    stop("y must be a matrix")
  }
  if(!is.logical(pe)){
    stop("pe must be a logical value")
  }

  # calculate the system matrices based on the parameter vector
  syst <- param_syst(theta, n, s)

  # calculate Sigma representing the variance of the measurement noise based on the parameter vector
  Sigma <- syst$Sigma

  # calculate Q representing the variance of the state noise based on the parameter vector
  Q <- syst$K %*% Sigma %*% t(syst$K)

  # calculate the Kalman filter output
  # Compute Kalman filter output
  # kf_out <- kf(y=y, syst=syst, pe=pe)
  # alternative using package FKF 
  Ctilde = cbind(syst$C,diag(s))
  Atilde = rbind(cbind(syst$A,syst$K),matrix(0,s,n+s))
  Rtilde = rbind(matrix(0,n,s),diag(s))
  Qtilde = syst$Sigma
  
  Sig <- Rtilde %*% Qtilde %*% t(Rtilde)
  ll = 0 
  ev <- max(abs(eigen(Atilde)$values))
  if (ev>0.99){
    # correct non-stability
    Atilde <- Atilde /ev*0.99
    ll <- ll+10^6*(exp(ev-0.99)-1)
  }
  
  P0 <- solve_lyap(A=Atilde, Q=Sig)
  T <- dim(y)[2]
  kf_out <- FKF::fkf(a0 = rep(0,n+s), P0 = P0, dt = matrix(0,n+s,T), ct =  matrix(0,s,T), Tt = Atilde, Zt = Ctilde, HHt= Sig, GGt = diag(s)*0, y= y)
  
  # return the negative log likelihood
  ll <- ll+  (-2)*kf_out$logLik
  return (ll)
}





#'
#' Calculate the Kalman filter output
#'
#' This function calculates the Kalman filter output using the parameters theta, the number of states n, the number of outputs s, the data y and a Boolean pe as input.
#' It returns the Kalman filter output as a list containing the estimated states, the residuals, the innovation variance, the Kalman gain and the negative log likelihood.
#'
#' @param theta list of parameters
#' @param n number of states
#' @param s number of outputs
#' @param y matrix of observations
#' @param pe Boolean whether to use prediction error formulation. Default is `TRUE`.
#'
#' @return list of Kalman filter results
#'
#'
#' @export
#'
cal_kf <- function(theta, n, s, y, pe=TRUE){
  # Get system matrices based on parameters
  syst <- param_syst(theta, n, s)
  T <- dim(y)[2]
  # Extract measurement noise variance
  Sigma <- syst$Sigma

  # Calculate state noise variance
  Q <- syst$K %*% Sigma %*% t(syst$K)

  # Compute Kalman filter output
  # kf_out <- kf(y=y, syst=syst, pe=pe)
  # alternative using package FKF 
  Ctilde = cbind(syst$C,diag(s))
  Atilde = rbind(cbind(syst$A,syst$K),matrix(0,s,n+s))
  Rtilde = rbind(matrix(0,n,s),diag(s))
  Qtilde = syst$Sigma
  Sig <- Rtilde %*% Qtilde %*% t(Rtilde)
  ll = 0 
  ev <- max(abs(eigen(Atilde)$values))
  if (ev>0.99){
    # correct non-stability
    Atilde <- Atilde /ev*0.99
    ll <- ll+10^6*(exp(ev-0.99)-1)
  }
  
  P0 <- solve_lyap(A=Atilde, Q=Sig)
  T <- dim(y)[2]
  kf_out <- FKF::fkf(a0 = rep(0,n+s), P0 = P0, dt = matrix(0,n+s,T), ct =  matrix(0,s,T), Tt = Atilde, Zt = Ctilde, HHt= Sig, GGt = diag(s)*0, y= y)
  kf_out$LL <- kf_out$logLik*(-2)
  
  # Return Kalman filter results
  return (kf_out)
}




#'
#' Quasi Maximum Likelihood Estimation (qMLE)
#'
#' qML estimation for system of order n using CVA initial estimate
#'
#' @param y matrix of observed data where rows represent variables and columns represent time points
#' @param n integer representing the order of the system to estimate
#' @param kmax integer representing the maximum lag to use for the initial autoregressive model fit
#' @param pe logical. If TRUE, the prediction error likelihood is used, otherwise the full likelihood is used. Defaults to `TRUE`.
#' @param print.level integer representing the level of `nlm` output to print. Defaults to `0` (print nothing).
#' @param iterlim integer; maximal number of iterations.   
#'
#' @return list containing the estimated parameters; output of nlm..
#'
#'
#' @export
#'
est_qml <- function(y, n, kmax, pe=TRUE, print.level=0, iterlim = 100){

  y_mat <- y
  if (inherits(y, "tbl_ts")){
    y_mat <- tsibble_to_mat(y)
  } 
  if (inherits(y_mat,"matrix") == FALSE){
    stop('Input to est_qml needs to be either a matrix or an appropriate data structure.')
  }  
  
  if (n <= 0 || kmax <= 0) {
    stop("n, kmax must be positive integers")
  }
  if (!is.logical(pe)) {
    stop("pe must be a logical value")
  }
  if (print.level < 0 || print.level > 2) {
    stop("print.level must be an integer between zero and two")
  }

  ### get AR order using AIC.
  ar <- fit_ar(y_mat, kmax=kmax)  # Fit autoregressive model to determine optimal lag
  k <- ar$k  # Retrieve the optimal lag order
  if(k==0){ # avoid picking AR order 0.
    k <- 1  # Set minimum lag order to 1 if 0 is chosen
  }

  # get CVA estimate
  f <- 2*k  # Set future horizon based on lag order
  p <- 2*k  # Set past horizon based on lag order

  # check, if conflicting with n.
  s <- dim(y_mat)[1]  # Get the number of variables
  if (min(f*s, p*s<n)){  # Ensure horizons are large enough
    f <- n  # Adjust future horizon to state order
    p <- n  # Adjust past horizon to state order
  }
  syst <- cal_CVA(y_mat, f, p, n)  # Perform CVA to estimate system matrices

  # convert into canonical form
  theta_init <- syst_param(syst=syst)  # Convert system matrices to parameter vector

  # add parameters for Sigma.
  theta_hat <- stats::nlm(cal_crit, p=theta_init, n=n, s=dim(y_mat)[1], y=y_mat, pe=pe , print.level=print.level, iterlim = iterlim, hessian =TRUE)  # Optimize parameters using non-linear minimization

  # Return the optimized parameter estimates
  return(theta_hat)
}


#'
#' Present Quasi Maximum Likelihood Estimator (qMLE)
#'
#' summary of qMLE output  
#'
#' @param theta_hat list containing the estimated parameters; output of nlm. 
#' @param n integer, system dimension 
#' @param s integer; output dimension 
#'
#' @return output on screen.
#'
#' @export
#'
summary_qmle <- function ( theta_hat,n=1, s=1){
  
  syste <- param_syst(theta_hat$estimate, n=n, s=s)
  theta_sd <- sqrt(diag(solve(theta_hat$hessian)))
  
  syst_sd <- param_syst(theta_sd, n=n, s=s,gradient=TRUE)
  
  writeLines("State Space Model")
  writeLines("A =")
  print_est_sd(
    est_sd = cbind(syste$A, syst_sd$A)  
  )
  writeLines("K =")
  print_est_sd(
    est_sd = cbind(syste$K, syst_sd$K)
  )
  writeLines("C =")
  print_est_sd(
    est_sd = cbind(syste$C, syst_sd$C)  
  )
}
