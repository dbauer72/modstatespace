#'
#' Convert an ARMA system to an equivalent state space system.
#'
#' This function converts an ARMA system specified in a list with matrices 'a' and 'b' and innovation variance 'Sigma' to an innovation 
#' state space system (F,H,K) and Sigma. 
#'
#' @param arma  a list of matrices a and b and a variance matrix Sigma. 
#'
#' @return syst a \code{StSp} object.
#'
#' @export
#'
arma2StSp <- function(arma){
  
  a <- arma$a
  b <- arma$b
  
  s <- dim(a)[1]
  
  if (s == 1){ # univariate system 
    p <- dim(a)[2]-1
    q <- dim(b)[2]-1
    n <- p+q
    
    # C matrix
    C <- matrix(0,1,n)
    if (p>0){
      C[1,1:p] <- -a[,-c(1)]
    }
    
    if (q>0){
      C[1,(p+1):(n)] <- b[,-1]
    }
    # A matrix
    A <- matrix(0,n,n)
    
    if (p>0){
      A[1,] <- C[1,]
      if (p>1){
        A[2:(p),1:(p-1)] <- diag(p-1)
      }
    }
    if (q>1){
      A[(p+2):n,(p+1):(n-1)] <- diag(q-1)
    }
    
    # K matrix
    K <- matrix(0,n,1)
    K[1,1]<-1
    if (q>0){
      K[p+1,1]<-1
    }
    
  } else { # multivariate system
    
    ps <- dim(a)[2]-s
    qs <- dim(b)[2]-s
    n <- ps+qs
    
    # C matrix
    C <- matrix(0,s,n)
    if (ps>0){
      C[1:s,1:ps] <- -a[,-c(1:s)]
    }
    
    if (qs>0){
      C[1:s,(ps+1):(n)] <- b[,-c(1:s)]
    }
    # A matrix
    A <- matrix(0,n,n)
    
    if (ps>0){
      A[1:s,] <- C[1:s,]
      A[(s+1):(ps),1:(ps-s)] <- diag(ps-s)
    }
    if (qs>0){
      A[(ps+s+1):n,(ps+1):(n-s)] <- diag(qs-s)
    }
    
    # K matrix
    K <- matrix(0,n,s)
    K[1:s,]<-diag(s)
    K[ps+c(1:s),]<-diag(s)
  }
  
  Sigma <- arma$Sigma 
  if (is.null(Sigma)){
    Sigma <- diag(s)
  }
  
  syst <- StSp(A, C, K,Sigma = Sigma,type ="innov")
  
  # model reduction 
  syst <- model_reduce(syst)
 

  return (syst)
}
