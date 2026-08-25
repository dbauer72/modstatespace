

#'
#' Canonical Variate Analysis (CVA) subspace estimation
#'
#' see Bauer, D.; Buschmeier, R. Asymptotic Properties of Estimators for Seasonally Cointegrated
#'  State Space Models Obtained Using the CVA Subspace Method. Entropy 2021, 23, 436
#'
#' This function performs Canonical Variate Analysis (CVA) to estimate a state space model.
#' It calculates the system matrices in canonical form based on the input data and specified parameters.
#'
#' @param y matrix of observed data where rows represent variables and columns represent time points.
#' @param f integer representing the future time horizon.
#' @param p integer representing the past time horizon.
#' @param n integer representing the number of states to estimate (state order). Must be less than or equal to `f` and `p`.
#'
#' @return syst in \code{StSp} object form. 
#'
#' @examples
#' set.seed(123)
#' y <- as.matrix(t(stats::arima.sim(n = 100, model = list(ar = c(0.5, -0.2, 0.3)), sd = 1)))
#' cal_CVA(y, f=5, p=5, n=1)
#'
#' @export
#'
cal_CVA <- function(y, f, p, n){
  
  # check parameters
  if (!is.matrix(y)) {
    stop("y must be a matrix")
  }
  if (f <= 0 || p <= 0 || n <= 0) {
    stop("f, p, and n must be positive integers")
  }
  kmax <- f+p  # Calculate maximum lags
  s <- dim(y)[1]  # Number of variables in the system
  Teff <- dim(y)[2]-f-p  # Effective sample size after accounting for lags
  
  # set up matrices
  Yk <- mlag(y, kmax)  # Set up matrices for lagged observations
  Yf <- matrix(Yk[1:(s*f), ], ncol=Teff)  # Extract future observations matrix
  Yp <- matrix(Yk[s*f + c(1:(s*p)), ], ncol=Teff)  # Extract past observations matrix
  Hfp <- Yf %*% t(Yp)  # Cross-covariance matrix between future and past observations
  Wf <- chol(Yf %*% t(Yf))  # Cholesky decomposition of future covariance matrix
  Wp <- chol(Yp %*% t(Yp))  # Cholesky decomposition of past covariance matrix
  
  beta <- solve(t(Wf)) %*% Hfp %*% solve(Wp)  # Calculate beta using the inverse of the Cholesky factors
  
  # SVD, singual value decomposition
  usv <- svd(beta)  # Singular value decomposition of beta
  Khat <- t(usv$v[, 1:n]) %*% solve(t(Wp))  # Estimate state transition matrix Khat
  xhat <- Khat %*% Yp  # Estimate state sequence
  
  # estimate system matrices
  yeff <- matrix(Yk[(f-1)*s + (1:s), ], ncol=Teff)  # Effective observations
  Chat <- yeff %*% t(xhat) %*% solve(xhat %*% t(xhat))  # Estimate observation matrix Chat
  res <- yeff[, 1:(Teff-1)] - Chat %*% xhat[, 1:(Teff-1)]  # Residuals of the state estimation
  xh <- matrix(xhat[, 1:(Teff-1)], ncol=Teff-1)  # Lagged state matrices
  xh1 <- matrix(xhat[, 2:(Teff)], ncol=Teff-1) # Transition the state using the lagged state matrix xh
  Ahat <- xh1 %*% t(xh) %*% solve(xh %*% t(xh))  # Estimate state transition matrix Ahat
  Khat <- xh1 %*% t(res) %*% solve(res %*% t(res))  # Estimate disturbance matrix Khat
  Sigma <- res %*% t(res)/Teff  # Covariance matrix of residuals
  
  # write results into a StSp system 
  syst <- StSp(A=Ahat, C= Chat, K=Khat, Sigma = Sigma, type = "innov")
  
  # Convert to canonical form
  syst <- conv_syst(syst)
  
  # Return system matrices and covariance
  return(syst)
}


#'
#' Extract Parameters for Estimated Covariance Matrix
#'
#' This function extracts the parameters from the estimated covariance matrix `Sig` by performing a Cholesky decomposition
#' and returns them in vector form. The parameters are extracted from the lower triangular part of the Cholesky factor.
#'
#' @param Sig matrix representing the estimated covariance matrix of the innovations.
#'
#' @return vector containing the extracted parameters from the Cholesky factor of the covariance matrix.
#'
#' @examples
#' Sig <- matrix(c(2, 0.5, 0.5, 1), nrow=2, ncol=2)
#' cov_Sig_par(Sig)
#'
#' @export
#'
cov_Sig_par <- function(Sig){
  
  # check parameters
  if(!is.matrix(Sig)){
    stop("Sig must be a matrix")
  }
  
  s <- dim(Sig)[1] # get dimension of innovation variance
  L <- t(chol(Sig)) # get Cholesky factor of innovation variance
  theta_l <- c()
  for (j in 1:s){
    theta_l <- c(theta_l, L[j, 1:j]) # extract elements from Cholesky factor
  }
  return(theta_l)
}

