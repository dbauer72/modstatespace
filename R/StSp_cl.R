
#' constructor method for StSp object
#'
#' @param syst StSp object
#' @param A   nxn matrix 
#' @param C   sxn matrix
#' @param K   nxs matrix
#' @param Sigma   sxs innovation variance matrix
#' @param Q   nxn matrix of noise variance in state equation
#' @param R   sxs matrix of noise covariance in observaiton equation
#' @param S   nxs matrix of cross covariance
#' @param type  string; default: "innov" for innovation representation. Alternative "UC§ for unobserved components. 
#'
#' @return syst StSp object
#'
#' @export
#'
new_StSp <- function(syst,A=as.matrix(0.5), C=as.matrix(1), K=as.matrix(1),Sigma = as.matrix(1),Q=as.matrix(1),R=as.matrix(1),S = as.matrix(1),type ="innov")
{
  stopifnot(is.list(syst))
  structure(syst, 
            class = "StSp",
            A=A, 
            C=C, 
            K=K,
            Sigma = Sigma,
            Sigma_chol = t(chol(Sigma)), 
            Q=Q,
            R=R,
            S=S
  )
}


#' validator method for verifying the created StSp object is valid.
#'
#' @param syst  StSp object
#'
#' @return syst StSp object
#' 
#' @export
#'
validate_StSp <- function(syst) {
  
  if (syst$type == 'innov') {
    if (is.matrix(syst$A)*is.matrix(syst$C)*is.matrix(syst$K) == FALSE){
      stop(
        "All entries of syst must be matrices!",
        call. = FALSE
      )
    }
    if (min(eigen(syst$Sigma)$values)<0){
      stop("Sigma must be positive semidefinite.")
      call. = FALSE
    }
    if (is.null(syst$Q)){
      syst$Q = syst$K %*% syst$Sigma %*% t(syst$K)
    }

    if (is.null(syst$R)){
      syst$R = syst$Sigma
    }
    if (is.null(syst$S)){
      syst$S = syst$K %*% syst$Sigma 
    }
  }
  
  return(syst)
}


#' user helper function: supports users to define the StSp object. 
#'
#' @param A   nxn matrix 
#' @param C   sxn matrix
#' @param K   nxs matrix
#' @param Sigma   sxs innovation variance matrix
#' @param Q   nxn matrix of noise variance in state equation
#' @param R   sxs matrix of noise covariance in observaiton equation
#' @param S   nxs matrix of cross covariance
#' @param type  string; default: "innov" for innovation representation. Alternative "UC§ for unobserved components. 
#'
#' @return syst StSp object
#' 
#' @export
#'
StSp <- function(A=as.matrix(0.5), C=as.matrix(1), K=as.matrix(1),Sigma = as.matrix(1),Q=as.matrix(1),R=as.matrix(1),S = as.matrix(1),type ="innov") {

  if (type == "innov"){
    syst <- list(A=A,C=C,K=K,Sigma=Sigma,Sigma_chol = t(chol(Sigma)),Q=NULL,R=NULL,S=NULL,type=type)
  }
  if (type == "UC"){
    syst <- list(A=A,C=C,K=NULL,Sigma=Sigma,Q=Q,R=R,S=S,type=type)
  }
  validate_StSp(new_StSp(syst,A,C,K,Sigma,Q,R,S))
}


#' print function: prints a StSp system object in a nice way. 
#'
#' @param x StSp object
#' @param ... additional arguments passed on. 
#'
#' @return on screen information
#'
#' @export
#' @method print StSp
#'
print.StSp <- function(x,...){
  
  if (is.null(names(x)) == FALSE){
    x_v <- list(1)
    x_v[[1]] <- x
  } else {
    x_v <- x
  }
  for (jx in 1:length(x_v)){
    syst <- x_v[[jx]]
    
    writeLines("State Space Model")
    cat('Type of state space system:',syst$type,'\n')
    dims <- dim(syst$C)
    cat(c('s:',dims[1],', n:',dims[2]), "\n")
    
    if (syst$type == "innov"){
      writeLines("A =")
      print_est_sd(
        est_sd = cbind(syst$A, syst$A*NA)  
      )
      writeLines("K =")
      print_est_sd(
        est_sd = cbind(syst$K, syst$K*NA)
      )
      writeLines("C =")
      print_est_sd(
        est_sd = cbind(syst$C, syst$C*NA)  
      )
      writeLines("Noise covariance: Sigma =")
      print_est_sd(
        est_sd = cbind(syst$Sigma, syst$Sigma*NA)  
      )
    } else {
      
      writeLines("A =")
      print_est_sd(
        est_sd = cbind(syst$A, syst$A*NA)  
      )
      writeLines("C =")
      print_est_sd(
        est_sd = cbind(syst$C, syst$C*NA)  
      )
      writeLines("State noise variance: \n Q =")
      print_est_sd(
        est_sd = cbind(syst$Q, syst$Q*NA)
      )
      
      writeLines("Covariance of noises: S =")
      print_est_sd(
        est_sd = cbind(syst$S, syst$S*NA)
      )    
      writeLines("Obervations noise variance: R =")
      print_est_sd(
        est_sd = cbind(syst$R, syst$R*NA)
      )
    }
  }
}

#' generic method for plotting the impulse response function 
#'
#' @param syst StSp object
#' @param M integer; number of lags to include
#' @param plots integer; 0... no plots, 1... static plots, 2... interactive plot. 
#'
#' @return on screen information
#'
#' @export
impulse <- function(syst, M, plots)
{
  #print("Calling the base impulse function")
  UseMethod("impulse")
  print(" this is not executed")
}


#'
#' Calculate impulse response for a system
#' System must be defined as a StSp object. 
#'
#' @param syst system matrices in StSp object form. 
#' @param M number of impulses to calculate
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return list of matrices with impulse response
#'
#' @export 
#' @method impulse StSp
#'
impulse.StSp <- function(syst, M=20,plots=0){
  
  IR <- list()
  
  if (is.null(names(syst)) == FALSE){
    syst_vec <- list(1)
    syst_vec[[1]] <- syst
  } else {
    syst_vec <- syst
  }
  J <- length(syst_vec)
  IR <- list()
  
  for (js in 1:J) {
    syst <- syst_vec[[js]]
    A=syst$A
    C=syst$C
    
    
    # dimensions
    s <- dim(C)[1]
    n <- dim(C)[2]
    
    # initialize IM
    IR[[js]] <- list()

    
    if (syst$type == "innov"){
      K=syst$K
      IR[[js]][[1]] <- diag(s)
    } else {
      K = t(chol(syst$Q))
      IR[[js]][[1]] <- matrix(0,s,s+dim(K)[2])
      IR[[js]][[1]][,1:s] <- diag(s)
    }
    
    # iterate to calculate impulse response
    cur <- diag(n)
    null <- matrix(0,s,s)
    if (syst$type == "innov"){   
      for (j in 1:M){
        IR[[js]][[j+1]] <- C %*% cur %*% K
        cur <- cur %*% A
      } 
    } else {
      for (j in 1:M){
        IR[[js]][[j+1]] <- cbind(null,C %*% cur %*% K)
        cur <- cur %*% A
      } 
    }
  }
  if (plots>0){
    #plot_inter(IR, main = 'Impulse response', ylab = 'IR[r,c,.]',interactive = (plots>1))
    p<- plot_inter(syst_vec, type= 'IR', interactive = (plots>1))
    p
  }
  return (IR)
}


#'generic method for covariance function 
#'
#' Calculate first M covariances for matrix system
#' Remark: F matrix needs to be stable. 
#'
#' @param syst system matrices in StSp object form. 
#' @param M number of impulses to calculate
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return list of matrices with covariance sequence
#'
#' @export
covariance_function <- function(syst, M,plots)
{
  #print("Calling the base covariance function")
  UseMethod("covariance_function")
  #print(" this is not executed")
}


#'
#' Calculate first M covariances for matrix system
#' F matrix needs to be stable. 
#'
#' @param syst system matrices in StSp object form. 
#' @param M number of impulses to calculate
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return list of matrices with covariance sequence
#'
#' @export 
#' @method covariance_function StSp
#'
covariance_function.StSp <- function(syst,M=20,plots=0){
  
  if (is.null(names(syst)) == FALSE){
    syst_vec <- list(1)
    syst_vec[[1]] <- syst
  } else {
    syst_vec <- syst
  }
  J <- length(syst_vec)
  CV <- vector("list",J)
  
  for (js in 1:J) {
    syst <- syst_vec[[js]]
    if (syst$type == "innov"){
      # state variance P from Ljapunov equation.
      ev <- eigen(syst$A)$values
      if (min(abs(ev - 1)) < 0.0001) {
        # warn for eigenvalues
        warning("Covariance_function: A has eigenvalues at z=1! Adjusting!")
        A <- syst$A/ev*0.99
      } else {
        A <- syst$A
      } 
      
      n <- dim(syst$A)[1]
      Akron <- diag(n ^ 2) - kronecker(A, A)
      Om <- syst$K %*% syst$Sigma %*% t(syst$K)
      vP <- solve(Akron) %*% matrix(Om, ncol = 1)
      P <- matrix(vP, ncol = n)
      P <- (P + t(P)) / 2 # ensure that P is symmetric
      # Variance
      CV[[js]] <- vector("list",M+1)
      CV[[js]][[1]] <- syst$Sigma + syst$C %*% P %*% t(syst$C)
        
      # M
      Mm <- syst$K %*%  syst$Sigma + syst$A %*% P %*% t(syst$C)
      cur <- diag(n)
      for (j in 1:M) {
        CV[[js]][[j + 1]] <- syst$C %*% cur %*% Mm
        cur <- cur %*% syst$A
      }
    }
    if (syst$type == "UC"){
      ev <- max(abs(eigen(syst$A)$values))
      if (ev > 0.9999) {
        # warn for eigenvalues
        warning("Covariance_function: A has eigenvalues on the unit circle! Adjusting.")
        A <- syst$A/ev*0.99
      } else {
        A <- syst$A
      } 

      n <- dim(syst$A)[1]
      Akron <- diag(n ^ 2) - kronecker(A, A)
      vP <- solve(Akron) %*% matrix(syst$Q, ncol = 1)
      P <- matrix(vP, ncol = n)
      P <- (P + t(P)) / 2 # ensure that P is symmetric
      # Variance
      CV[[js]] <- vector("list",M+1)
      CV[[js]][[1]] <- syst$R + syst$C %*% P %*% t(syst$C)
        
      # M
      Mm <- syst$S + syst$A %*% P %*% t(syst$C)
      cur <- diag(n)
      for (j in 1:M) {
        CV[[js]][[j + 1]] <- syst$C %*% cur %*% Mm
        cur <- cur %*% syst$A
      }
    }
  }
  
  attr(CV,"cM") <- Mm 
  if (plots>0){
    p<- plot_inter(syst_vec, type= 'CV', interactive = (plots>1))
    p
  }
  return(CV)
}


#'generic method for transfer function (including plotting) 
#'
#' @param syst system matrices in StSp object form. 
#' @param M integer; number of frequencies
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return list of matrices with transfer function at Fourier frequencies
#'
#' @export
transfer_function <- function(syst, M, plots)
{
  #print("Calling the base covariance function")
  UseMethod("transfer_function")
  #print(" this is not executed")
}

#'generic method for transfer function (including plotting) 
#'
#' @param syst system matrices in StSp object form. 
#' @param M integer; number of points where to evaluate transfer function (default: M=100)
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return list of matrices with transfer function at Fourier frequencies
#'
#' @export 
#' @method transfer_function StSp
#'
transfer_function.StSp <- function(syst,M=100,plots=0){
  
  if (is.null(names(syst)) == FALSE){
    syst_vec <- list(1)
    syst_vec[[1]] <- syst
  } else {
    syst_vec <- syst
  }
  J <- length(syst_vec)
  TF <- vector("list",J)
  freq <- seq(from=0,to = pi, by = pi/M)
  zj <- exp(1i*freq)
  
  for (js in 1:J) {
    syst <- syst_vec[[js]]
    TF[[js]] <- matrix(0,M+1,1) # norm of entry. 
    
    if (syst$type == "innov"){
      A <- syst$A
      C <- syst$C
      K <- syst$K 
      s <- dim(C)[1]
      Id <- diag(dim(A)[1])
      for (j in 1:length(freq)){
        zjj <- zj[j]
        TF[[js]][j] <- sum(abs((diag(s) + zjj* C %*% solve(Id- zjj*A) %*% K )))
      }
    } else {
      A <- syst$A
      evs <- max(abs(eigen(A)$values))
      if (evs>0.99999){
        message('TF: A is not stable, correcting!')
        A <- A/evs*0.99999
      }
      C <- syst$C
      Q <- syst$Q 
      K <- chol(syst$Q)
      s <- dim(C)[1]
      sin <- dim(Q)[2]+s
      Id <- diag(dim(A)[1])
      for (j in 1:length(freq)){
        zjj <- zj[j]
        TF[[js]][j] <- sum(abs( diag(s) ))
        for (jj in 1:sin){
          TF[[js]][j] <- TF[[js]][j] + sum(abs( zjj* C %*% solve(Id- zjj*A) %*% K ))
        }
      }
      
    }
  }
  
  if (plots>0){
    data <- data.frame(fr = freq,tf = TF[[1]], syst = 1)

    if (J>1){
      for (js in 1:J) {
        data <- rbind(data,data.frame(fr = freq,tf = TF[[js]], syst = js))
      } 
    }
    
    p <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data$fr, y = .data$tf, color = factor(syst))) +
      ggplot2::geom_line() + ggplot2::labs( x = "Frequency", y = "Log of modulus of transfer function", title = "TF plot") + ggplot2::scale_y_continuous(trans='log10')
    print(p) 
  }
  return(TF)
}


#'generic method for filtering  
#'
#' Filter using a state space system 
#'
#' @param syst system matrices in StSp object form. 
#' @param u time series to be filtered 
#' @param factor double; scaling the output 
#' @param xinit initial state (if supplied; default zero initialisation)
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return filtered series in tsibble form. 
#'
#' @export
lin_filter <- function(syst, u, xinit=NULL, factor = 1, plots)
{
  #print("Calling the base covariance function")
  UseMethod("lin_filter")
  #print(" this is not executed")
}


#'generic method for filtering  
#'
#' Filter using a state space system 
#'
#' @param syst system matrices in StSp object form. 
#' @param u time series to be filtered 
#' @param factor double; scaling the output 
#' @param xinit initial state (if supplied; default zero initialisation)
#' @param plots integer; 0.. no plots; 1.. static plots, 2.. interactive plots
#'
#' @return filtered series in tsibble form. 
#'
#' @export
#' @method lin_filter StSp
#'
lin_filter.StSp <- function(syst,u,xinit=NULL,factor = 1, plots=0){
  
  # extract time series 
  u_mat <- tsibble_to_mat(u)
  T <- dim(u_mat)[1]
  
  # initialize filter 
  A <- syst$A
  K <- syst$K
  C <- syst$C 
  n <- dim(A)[1]
  
  if (length(xinit) != n){
    xinit <- rep(0,n)
  }
  # run the recursion
  x_mat <- matrix(0,T,n)
  x_mat[1,] <- xinit
  
  y_hat <-u_mat*0 
  y_hat[1,] <- factor * (x_mat[1,] %*% t(C) + u_mat[1,])
  
  for (t in 2:T){
    x_mat[t,] <- x_mat[t-1,] %*% t(A) + u_mat[t-1,] %*% t(K)
    y_hat[t,] <- factor * (x_mat[t,] %*% t(C) + u_mat[t,])
  }
  
  # provide output 
  labs_rep <- vector()
  labs <- paste('y_hat')
  labs_rep <- rep(labs,T)

  y_hat_ts <- tsibble::tsibble(y = y_hat, time = u$time, vars =labs_rep, index = 'time', key = 'vars')
  
  # plot, if wanted 
  if (plots>0){
    comb <- dplyr::bind_rows(u,y_hat_ts)
    p <- ggplot2::autoplot(comb) +
      ggplot2::labs(
        title = "Time Series and filtered version",
        x = "Time",
        y = "Value",
        color = "System"
      ) +  ggplot2::theme_minimal()
    
    print(p) 
  }
  return( y_hat_ts )
    
}


#' solve_shift provides the estimate of F using the shift invariance structure of cO.
#' 
#' @param cO observability matrix
#' @param s integer, output dimension
#' 
#' @return A estimated matrix A
#' 
#' @export
#'
solve_shift <- function(cO,s=1){
  
  krow <- dim(cO)[1]
  OtO <- t(cO[1:(krow-s),]) %*% cO[1:(krow-s),]
  
  # regularize, if not of full rank 
  evs <- eigen(OtO)$values
  if (min(abs(evs)) < 0.00000001){
    OtO <- OtO + diag(dim(OtO)[1])
  }
  
  A <- solve( OtO) %*% (t(cO[1:(krow-s),]) %*% cO[(s+1):krow,])
  
  return(A)
}

#'
#' Reduces the model to a minimal system eliminating non-obervability and non-controllabilty. 
#' REMARK: Only works for systems in innovation form. 
#'
#' @param syst system matrices in StSp object form. 
#' @param tol  cut of threshhold: singular values below this threshhold are neglected. Default set to 0.0000001. 
#'
#' @return reduced order system
#'
#' @export
#'
model_reduce <- function(syst,tol=0.0000001){
  if (syst$type == "innov"){
    n <- dim(syst$A)[1]
    s <- dim(syst$C)[1]
    cH <- Hankel(syst,krow=4*n,kcol=4*n)
    qrH <- qr(cH, LAPACK = FALSE,tol=tol)
    cO <- qr.Q(qrH)
    cC <- qr.R(qrH)
    
    nred = qrH$rank
    syst$C = cO[1:s,1:nred, drop=FALSE]
    syst$K = cC[1:nred,1:s, drop=FALSE]
    syst$A = solve_shift(cO[,1:nred, drop=FALSE],s)
    
    # check the reduction error in Hankel matrix norm 
    cH_red <- Hankel(syst,krow=4*n,kcol=4*n)
    cat('Reduction error in Hankel norm: ',norm(cH-cH_red),'\n')
  } else {
    # the system is in general form. 
    A <- syst$A
    C <- syst$C
    Q <- syst$Q
    S <- syst$S
    
    # dimensions 
    n <- dim(A)[1]
    s <- dim(C)[1]
    
    # first reduction: eliminate non-observability 
    cO = matrix(0,4*n*s,n)
    cO[1:s,] <- C
    cur <- C %*% A
    for (j in 2:(4*n)){
      cO[(s*(j-1)+c(1:s)),] <- cur
      cur <- cur %*% A
    }
    
    qrcO <- qr(cO, LAPACK = FALSE, tol= tol)
    nred <- qrcO$rank
    cOr <- qr.Q(qrcO)
    cCr <- qr.R(qrcO)
    Tr <- cCr[1:nred,]

    Cr  <- cOr[1:s,1:nred]
    Qr  <- Tr %*% Q %*% t(Tr)
    Ar  <- solve_shift(cOr[,1:nred],s)
    Sr  <- Tr %*% S
    
    # second reduction: adjust for non-controllability 
    Pr <- solve_lyap(A= Ar, Q = Qr)
    sv <- svd(Pr)
    nred <- sum(sv$d>tol)
    U <- sv$u[,1:nred]
    
    # fill in system
    syst$Q <- t(U) %*% Qr %*% U
    syst$S <- t(U) %*% Sr
    syst$A <- t(U) %*% Ar %*% U
    syst$C <- Cr %*% U
    
  }
  
  return(syst)
} 



#' Calculates the Hankel matrix. 
#'
#' @param syst system matrices in StSp object form. 
#' @param krow integer; number of block rows
#' @param kcol integer; number of block columns. 
#'
#' @return reduced order system
#'
#' @export
#'
Hankel <- function(syst, krow = 10, kcol=10){
  
  s = dim(syst$C)[1]
  n = dim(syst$C)[2]
 
  # construct Hankel matrix of impulse response coefficients
  cH = matrix(0, nrow = krow*s, ncol = kcol*s)
  max.lag = krow+kcol-1
  A = syst$A
  C = syst$C
  K = syst$K 
  
  curCA = C
  for (lag in 1:max.lag){
    kj = curCA %*% K
    for (i in (1:krow)){ 
      j <- lag-i+1
      if ((j>0)&(j<=kcol)){
        cH[((i-1)*s+1):(i*s),((j-1)*s+1):((j)*s)] = kj
      }
    }
    curCA = curCA %*% A
  }
  
  return(cH)
} 

#' Calculates the Toeplitz matrix. 
#'
#' @param Cv list of covariance matrices  
#' @param krow integer; number of block rows
#'
#' @return matrix Gamma_p. 
#'
#' @export
#'
Toeplitz <- function(Cv,krow = 10){
  
  s = dim(Cv[[1]])[1]
  # construct Toeplitz matrix of covariances 
  Gammap = matrix(0, nrow = krow*s, ncol = krow*s)
  
  for (rj in 1:krow){
    for (cj in 1:krow){
      if (rj >= cj){
        Ck <- Cv[[rj-cj+1]]
      } else {
        Ck <- t(Cv[[cj-rj+1]])
      }
      Gammap[ ((rj-1)*s)+c(1:s),((cj-1)*s)+c(1:s)] <- Ck
    }
  }
  
  return(Gammap)
}


#'
#' check for stability of an AR system (via the eigenvalues of the companion form)
#'
#' @param object system, \code{StSp} object
#'
#' @return stab_cond list containing the entries 'check' Boolean (TRUE, if stable) and 'poles' (eigenvalues of companion matrix)
#'
#' @export
#'
check_stability.StSp <- function(object){  
  
  A <- object$A
  evs <- base::eigen(A)$values 
  poles <- 1/evs
  lambdamax <- max(base::abs(evs))
  check = TRUE 
  if (lambdamax > 1){ check = FALSE } 
  
  stab_cond <- list(check = check, poles = poles)
  return( stab_cond ) 
}



#'
#' Convert system to canonical form. Input and output are stored in a StSp object. 
#'
#' @param syst a \code{StSp} object
#'
#' @return syst a \code{StSp} object
#'
#' @examples
#' ones_mat <- as.matrix(1,1,1)
#' syst <- StSp(A = ones_mat*0.5,C=  ones_mat, K = ones_mat, Sigma = ones_mat, type = "innov")
#' conv_syst(syst)
#'
#' @export
#'
conv_syst <- function(syst){
  
  stopifnot(inherits(syst,"StSp"))
  stopifnot(syst$type == "innov")
  
  n <- dim(syst$A)[1]  # number of states
  s <- dim(syst$C)[1]  # number of observations
  
  # calculate observability.
  O <- matrix(0, s*n, n)  # initialize observability matrix
  cur <- diag(n)  # initialize identity matrix of size n
  for (j in 1:n){  # iterate over number of states
    O[(j-1)*s+c(1:s), ] <- syst$C %*% cur  # calculate current observability block
    cur <- cur %*% syst$A  # update current state transition
  }
  
  Trafo <- O[1:n,]  # extract transformation matrix
  iTrafo <- solve(Trafo)  # calculate inverse of transformation matrix
  
  At <- Trafo %*% syst$A %*% iTrafo  # transform F matrix
  Ct <- syst$C %*% iTrafo  # transform H matrix
  Kt <- Trafo %*% syst$K  # transform K matrix
  
  # write matrices into system
  syst$A <- At
  syst$C <- Ct
  syst$K <- Kt

  # return transformed system matrices
  return(syst)
}


#'
#' Convert system to innovation representation. 
#'
#' @param syst a \code{StSp} object
#'
#' @return syst a \code{StSp} object
#'
#'
#' @export
#'
conv_to_innov <- function(syst){
  
  stopifnot(inherits(syst,"StSp"))
  if (syst$type == "innov"){
    cat("System already in innovation representation. ")
    return(syst)
  }
  
  n <- dim(syst$A)[1]  # number of states
  s <- dim(syst$C)[1]  # number of observations
  
  # calculate covariances 
  Cv <- covariance_function(syst,M=4*n) # covariance function up to lag 4*n as default. 
  cM <- attr(Cv,"cM")
  
  syst_cov <- syst
  syst_cov$K <- cM 
  
  # calculate Hankel matrix and Gamma_p 
  twon = 2*n
  Hfp <- Hankel(syst_cov,twon, twon)
  Gammap <- Toeplitz(Cv[[1]], twon)
  
  # estimate beta 
  evs <- eigen(Gammap)$values
  if (min(evs)<0){
    beta <- Hfp %*% solve(Gammap+ (0.1-min(evs))*diag(dim(Gammap)[1]))
  } else {
    beta <- Hfp %*% solve(Gammap)
  }
  
  # factor beta 
  sv <- svd(beta)
  cO <- sv$u[,1:n] %*% diag(sv$d[1:n])
  cC <- t(sv$v[,1:n])
  
  # estimate matrices 
  K <- cC[,1:s,drop=FALSE]
  C <- cO[1:s,,drop=FALSE]
  A <- solve_shift(cO,s=s)
  
  if (norm(cC[,(2*n-s):(2*n)])>0.0001){
    cat(' krow not large enough -> double! ')
  }
  
  # calculate Sigma, Q, R, S. 
  
  Sigma <- Gammap[1:s,1:s] - beta[1:s,, drop=FALSE] %*% Gammap %*% t(beta[1:s,,drop=FALSE]) 
  
  Q <- K %*% Sigma %*% t(K)
  S <- K %*% Sigma 
  R <- Sigma
  
  # return system. 
  syst$A <- A
  syst$C <- C
  syst$K <- K
  syst$Q <- Q
  syst$R <- R
  syst$S <- S
  syst$Sigma <- Sigma 
  
  syst$type = "innov"
  
  # return transformed system matrices
  return(syst)
}


