#' sim_syst simulates time series data from a state space system in innovation form. 
#' 
#' Time series data of length T is simulated using Gaussian innovations from the state
#' space system supplied. 
#' The state space system must be in innovation form. 
#' 
#' @param syst \code{StSp} object, which represents a state space model, either in innovation representation or in general form.
#' @param T  integer; sample length. 
#' @param stationary  Boolean indicating whether the simulated process should be stationary.           
#'             
#' @return y real sxT matrix of observations. 
#' 
#' @export
sim_syst <- function(syst,T= 100,stationary = FALSE){

  if (!inherits(syst,"StSp")) stop("syst must be a StSp object.")
  s = dim(syst$C)[1]
  n = dim(syst$C)[2]
  
  # initial state: zero initialisation
  x <- matrix(0,n,1)
  
  
  # generate residuals. 
  if (syst$type == "innov"){
    Sigma = syst$Sigma
    L = chol(Sigma)
    v = L %*% matrix(stats::rnorm(T*s),s,T)
    w = v
    K = syst$K
  } else {
    Sig = matrix(0,n+s,n+s)
    Sig[1:s,1:s] <- syst$R
    Sig[s+(1:n),s+(1:n)] <- syst$Q
    Sig[s+(1:n),1:s] <- syst$S
    Sig[1:s,s+(1:n)] <- t(syst$S)
    L = chol(Sig)
    vw = L %*% matrix(stats::rnorm(T*(s+n)),s+n,T)
    v = vw[1:s,,drop=FALSE]
    w = vw[s+(1:n),,drop=FALSE]
    K = diag(n)
  }
  
  # or stationary initialisation?
  if (stationary){
    stab_cond <- check_stability(syst) 
    if (stab_cond$check){
    if (syst$type == "innov"){
      Q = syst$K %*% Sigma %*% t(syst$K)
      F = syst$A
      P = solve_lyap(syst$A,Q)
    } else {
      P = solve_lyap(syst$A,syst$Q)
    }
    L_p = t(chol(P))
    x = L_p %*% stats::rnorm(n)
    } else {
      warning('System is unstable. No stationary solution exists. Starting with zero initial state instead!')
    }
    
  }
  
  
  # start the recursion. 
  y = v*0
  
  # iterate over sample size
  for (t in 1:T){
    y[,t] = syst$C %*% x + v[,t]
    x = syst$A %*% x + K %*% w[,t]
  }
  
  labs <- list(s)
  labs_rep <- vector()
  for (j in 1:s){
    labs[[j]] <- paste('y',j,sep='')
    labs_rep <- append(labs_rep,rep(labs[[j]],T))
  }
  time <- seq(from =1, by = 1, to = T)
  time_rep <- rep(time,s)
  y_rep <- matrix(t(y),ncol=1)
  
  yts <- tsibble::tsibble(y = y_rep, time = time_rep, vars =labs_rep, index = time, key = 'vars')

  return(yts)
}
