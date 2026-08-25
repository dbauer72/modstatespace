#'
#' Calculate the criterion function for minimum distance estimation as a function of the parameter vector
#'
#' This function calculates the criterion function as the squared distance between the directly estimated covariances 
#' and the corresponding covariances as a function of the parameter vector.
#' It takes the parameter vector theta, the number of states n, the number of outputs s, the data C, the lags contained in C and
#' a weighting matrix W as input.
#' It returns the corresponding squared difference. 
#'
#' @param theta vector of parameters
#' @param n number of states
#' @param s number of outputs
#' @param C list of estimated covariances. 
#' @param lags vector of integers containing the lags contained in C. 
#' @param W matrix for weighing the measurement errors in C. 
#'
#' @return criterion value of squared differences. 
#'
#'
#' @export
#'
cal_crit_min_dist <- function(theta, n, s, C, lags = NULL, W = NULL){
  
  if ((length(C) != length(lags))|(is.null(lags))){
    lags <- c(1:length(C))
  }
  
  # direct estimate
  pid <- matrix(0,length(lags)*s^2,1)
  for (j in 1:length(lags)){
    Cj <- C[[j]]
    pid[(j-1)*s^2+c(1:s^2)] <- matrix(Cj,ncol=1)
  }
  
  # covariance corresponding to parameter vector  
  pid_theta <-cov_seq_grad(theta,n=n,s=s,lags = lags,grad_bool = FALSE)
  
  # weights provided? 
  if (is.null(W)){
    W = diag(length(pid))
  }
  
  if (dim(W)[1] != length(pid_theta)){
    warning('cal_crit_min_dist: Weight matrix W does not have the right dimensions, adjusting.')
    W = diag(length(pid))
  }
  # calculate the criterion 
  dpi <- pid_theta - pid
  crit <- t(dpi) %*% W %*% dpi
  
  return (crit)
} 

#'
#' Calculate the solution to Ljapunov equations (and their derivatives) for a state space system. 
#'
#' This function calculates the solution to the Ljapunov equation P = F P F' + K Om K' as well as 
#' its derivative with respect to the parameter vector theta (for a system with output dimension s and state dimension n).
#' The output is given as a list of matrices. 
#'
#' @param theta vector of parameters
#' @param n number of states
#' @param s number of outputs
#' @param grad_bool if TRUE, also derivative is calculated
#'
#' @return list of matrices. First is solution, the rest the derivatives. 
#'
#' @export
#'
dLjapunov <- function(theta, n, s,grad_bool){
  
  systs <- param_syst_deriv(theta,n,s)
  syst <- systs[[1]]
  
  # calculate solution 
  n <- dim(syst$A)[1]
  Akron <- diag(n ^ 2) - kronecker(syst$A, syst$A)
  Om <- syst$K %*% syst$Sigma %*% t(syst$K)
  vP <- solve(Akron) %*% matrix(Om, ncol = 1)
  P0 <- matrix(vP, ncol = n)
  P0 <- (P0 + t(P0)) / 2 # ensure that P is symmetric
  
  dPs <- list()
  dPs[[1]] <- P0
  
  if (grad_bool){
    
    # now go for derivatives. 
    for (j in 1:length(theta)){
      dsyst <- systs[[j+1]]
      
      # collect matrices 
      A <- syst$A
      dA <- dsyst$A 
      K <- syst$K
      dK <- dsyst$K 
      Sig <- syst$Sigma
      dSig <- dsyst$Sigma
      
      # calculate the derivative of the right hand side 
      dQ <- dA %*% P0 %*% t(A) + A %*% P0 %*% t(dA) + dK %*% Sig %*% t(K) + K %*% dSig %*% t(K) + K %*% Sig %*% t(dK) 
      
      # solve the equation 
      dvP <- solve(Akron) %*% matrix(dQ, ncol = 1)
      dP <- matrix(dvP, ncol = n)
      # symmetrize 
      dP <- (dP + t(dP)) / 2 # ensure that P is symmetric
      
      # write result in list. 
      dPs <- append(dPs,list(dP))
    }
  }
  
  return( dPs ) 
} 


#'
#' Calculate covariances for lags given in 'lags'. 
#'
#' @param theta parameter vector for system. 
#' @param n integer state dimension
#' @param s integer output dimension 
#' @param lags vector of integers 
#' @param grad_bool Boolean indicating whether derivatives should be calculated as well. 
#'
#' @return matrix of vectorized covariances. First column actual values, remaining columns (if any) derivatives. 
#'
#' @export 
#'
cov_seq_grad <- function(theta,n,s,lags = 0,grad_bool = TRUE){

  lth <- 0 
  if (grad_bool){
    systs <- param_syst_deriv(theta,n,s)
    syst <- systs[[1]]
    dPs <- dLjapunov(theta,n,s,grad_bool=TRUE)
    P0 <- dPs[[1]]
    lth <- length(theta)
  } else {
    syst <- param_syst(theta,n,s)
    dPs <- dLjapunov(theta,n,s,grad_bool=FALSE)
    P0 <- dPs[[1]]
  }
    
  J <- length(lags)
  CV <- matrix(0,nrow = s*s*J,ncol = 1+lth)
  
  Mm <- syst$K %*%  syst$Sigma + syst$A %*% P0 %*% t(syst$C)
  
  ev <- eigen(syst$A)$values
  if (min(abs(ev - 1)) < 0.0001) {
    # warn for eigenvalues
    warning("cov_seq_grad: F has eigenvalues at z=1!")
  } 
    
  lags <- sort(lags) # to be sure, it is sorted.   
  cur <- 1 
  Fj <- diag(n)
  
  # first calculate the value of Ck 
  for (js in 1:J) {
    j <- lags[js]
    if (cur<j){
      for (jc in ((cur+1):j)){
        Fj <- Fj %*% syst$A
      }
      cur <- j
    }
    Ck <- matrix(0,nrow = s,ncol = s) 

    if (j == 0){ # lag 0 is different!
      Ck <- syst$Sigma + syst$C %*% P0 %*% t(syst$C)
    } else {
      Ck <- syst$C %*% Fj %*% Mm 
    }
    CV[(js-1)*s*s+c(1:(s*s)),1] <- matrix(Ck,ncol=1) 
  }
  
  
  # next go for the derivative
  if (grad_bool){
    
    for (jj in 1:lth){
      dP <- dPs[[jj+1]]
      dsyst <- systs[[jj+1]]
      dMm <- dsyst$K %*%  syst$Sigma + syst$K %*% dsyst$Sigma 
      dMm <- dMm + dsyst$A %*% P0 %*% t(syst$C) + syst$A %*% dP %*% t(syst$C) + syst$A %*% P0 %*% t(dsyst$C)
      
      cur <- 1 
      Fj <- diag(n)
      dFj <- Fj*0
      
      # first calculate the value of Ck 
      for (js in 1:J) {
        j <- lags[js]
        if (cur<j){
          for (jc in ((cur+1):j)){
            dFj <- dFj %*% syst$A + Fj %*% dsyst$A
            Fj <- Fj %*% syst$A
          }
          cur <- j
        }
        dCk <- matrix(0,nrow = s,ncol = s) 
        
        if (j == 0){ # lag 0 is different!
          dCk <- dsyst$Sigma + dsyst$C %*% P0 %*% t(syst$C) + syst$C %*% dP %*% t(syst$C) + syst$C %*% P0 %*% t(dsyst$C)
        } else {
          dCk <- dsyst$C %*% Fj %*% Mm + syst$C %*% dFj %*% Mm + syst$C %*% Fj %*% dMm 
        }
        CV[(js-1)*s*s+c(1:(s*s)),jj+1] <- matrix(dCk,ncol=1) 
      }
      
    }
  }
  
  return(CV)
}

#'
#' Find the largest state dimension that leads to an identifiable system
#' given the lags contained in the sequence.  
#'
#' @param s integer output dimension 
#' @param lags vector of integers 
#'
#' @return min_evs vector of minimal eigenvalues for different state dimensions.  
#'
#' @export 
#'
find_dim_ident <- function(s,lags = 0){
  
  J <- length(lags)
  
  # find the maximum possible dimension. 
  n_max <- floor( (J-1)*s/2)
    
  if (n_max<1){
    stop('find_dim_ident: Not even n=1 is identified!')
  }
  # iterate up to max dimension.  
  min_evs <- numeric(n_max)
  for (n in 1:n_max){
    npar <- 2*s*n + s*(s+1)/2
    theta <- stats::runif(npar)

    CV <-  cov_seq_grad(theta,n=n,s=s,lags = lags,grad_bool = TRUE) 
    G <- CV[,-1]
    GM <- t(G) %*% G
    evs <- eigen(GM)$values
    min_evs[n] <- sqrt(min(abs(evs)))
  }
  
  # return results 
  return(min_evs)
} 



#'
#' Minimum distance estimation (MD)
#'
#' MD estimation for system of order n matches estimated covariances. 
#'
#' @param n integer representing the order of the system to estimate
#' @param C list of estimated covariances 
#' @param lags integer vector containing the lgas contained in the list C.
#' @param W covariance matrix of estimated C. 
#' @param theta_init initial estimate. If not present, random initialisation is used. 
#' @param print.level integer representing the level of `nlm` output to print. Defaults to `0` (print nothing).
#'
#' @return list containing the estimated parameters and the criterion function. 
#'
#'
#' @export
#'
est_MD <- function(n, C, lags = 0, W=NULL , theta_init = NULL, print.level=0){
  
  if (print.level < 0 || print.level > 2) {
    stop("print.level must be an integer between zero and two")
  }
  
  C0 <- C[[1]]
  if (inherits(C0,'matrix') == FALSE){
    stop('C must contain matrices. ')
  }
  s <- dim(C0)[1]

  # number of parameters 
  npar <- 2*s*n + s*(s+1)/2; 
  
  # set W
  dim_pi <- s*s*length(C)
  if (is.null(W)){
    W <- diag(dim_pi)
  }
  # convert into canonical form
  
  if (is.null(theta_init)){
    theta_init <- stats::runif(npar) 
  }
  
  # estimate the system minimising minimum distance 
  theta_hat <- stats::nlm(cal_crit_min_dist, p=theta_init, n=n, s=s, C=C, lags = lags, W = W, print.level=print.level)  # Optimize parameters using non-linear minimization
  
  # Return the optimized parameter estimates
  return(theta_hat)
}
