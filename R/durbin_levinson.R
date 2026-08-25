#'
#' Durbin Levinson algorithm for the estimation of AR models recursive in the lag length. 
#'
#' `durbin_levinson` takes covariance sequence (list of covariance matrices) as input as well 
#' as a maximal lag length and calculates sequentially in lag a series of AR models fitted to the 
#' covariances using the Durbin Levinson recursion.
#' 
#' The multivariate version uses two intertwined recursions according to Hannan and Deistler (1988), section 5.7., p. 121. (5.6.13) 
#' (attributed to Whittle). 
#' 
#' Each AR models is characterised using (I) the estimated coefficients, (II) the innovation variance Sigma and (III) the PACC. 
#'
#' @param C list of covariances of size sxs. 
#' @param L integer. maximal lag length
#'
#' @return list of information for estimated AR systems. (I) (a_j(p)), (II) Sigma_p (III) rho_{p,p}.  
#'
#'
#' @export
#'
durbin_levinson <- function(C, L) {

  if (is.list(C)==FALSE){
    stop('durbin_levinson requires the covariance sequence as a list as input.')
  }
  
  p <- length(C)-1
  if (L>p){
    message('durbin_levinson: Covariance sequence needs at least as many entries as L. Correcting L.')
    L <- p 
  }
  
  s <- dim(C[[1]])[1]
  ar_systs <- list()
  ar_systs <- append(ar_systs,list(list(Sigmah = C[[1]], Phihj = diag(s), Phihh = diag(s),tPhihj = diag(s),tPhihh = diag(s),tSigmah = C[[1]])))
  
  tSigh <- C[[1]]
  Sigh <- C[[1]]
  Phihh <- diag(s)
  tPhihh <- diag(s)
  Phihjm1 <- diag(s)
  tPhihjm1 <- diag(s)
  
  for (h in 1:L){ # start recursion in k 
    # Deltah 
    Deltah <- C[[h+1]]
    if (h>1){
      for (j in 1:(h-1)){
        Deltah <- Deltah + Phihj[,j*s+c(1:s)] %*% C[[h+1-j]]
      }
    }
    
    # Phihh
    Phihh <- - Deltah %*% solve(tSigh)
    tPhihh <- - t(Deltah) %*% solve(Sigh)
    
    # Phihj
    Phihj <- matrix(0,s,(h+1)*s)
    Phihj[,1:s]<- diag(s)
    Phihj[,h*s+c(1:s)]<- Phihh

    tPhihj <- matrix(0, s,(h+1)*s)
    tPhihj[,1:s]<- diag(s)
    tPhihj[,h*s+c(1:s)]<- tPhihh
    
    if (h>1){
      for (j in 1:(h-1)){
        Phihj[,j*s+c(1:s)] <- Phihjm1[,j*s+c(1:s)] + Phihh %*% tPhihjm1[,(h-j)*s+c(1:s)]
        tPhihj[,j*s+c(1:s)] <- tPhihjm1[,j*s+c(1:s)] + tPhihh %*% Phihjm1[,(h-j)*s+c(1:s)]
      }
    }
    
    # Sigmah 
    Sigh <- (diag(s)-Phihh %*% tPhihh) %*% Sigh
    tSigh <- (diag(s)-tPhihh %*% Phihh) %*% tSigh
    
    # write results 
    ar_systs <- append(ar_systs,list(list(Sigmah = Sigh, Phihj = Phihj, Phihh = Phihh, tPhihj = tPhihj, tPhihh = tPhihh, tSigmah = tSigh)))
    
    # fill in last values 
    Phihjm1 <- Phihj
    tPhihjm1 <- tPhihj
  }
  return(ar_systs)
}

