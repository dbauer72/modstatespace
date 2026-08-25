globalVariables(c(".data"))


#'
#' Produce matrix containing `k` lags of matrix `y` stapled in block rows.
#'
#' `mlag` takes a matrix `y`, an integer `k` and returns a matrix,
#' where each column contains the `k` lags of the corresponding column in `y`.
#' First row of returned matrix start with column `k` of `y`, last row starts with column 1
#'
#' @param y matrix. Each column is one time period.
#' @param k integer. Lag length
#'
#' @return matrix. Matrix with `k` lags of `y`.
#'
#'
#' @export
#'
mlag <- function(y, k) {

  if(!is.matrix(y)) {
    stop("y must be a matrix")
  }

  s <- dim(y)[1] # number of rows (observations) in y
  T <- dim(y)[2] # number of columns (time points) in y

  if (k<=0) {
    stop(sprintf("Minimum lag length must be 1, but k=%s is smaller.", k))
  }
  if (k>=T) {
    stop(sprintf("Maximum lag length must be smaller than number of data points, but k=%s is not smaller than #time points: %s.", k, T))
  }

  Yk <- matrix(0, s * (k + 1), T - k) # initialize a matrix to store the lags
  Yk[1:s, ] <- y[, (k + 1):T] # copy the non-lagged values to the first s rows

  for (j in 1:k) { # loop over each lag
    Yk[(j * s) + 1:s, ] <- y[, (k - j) + 1:(T - k)] # fill in the lagged values
  }

  return(Yk) # return the matrix with lags
}

#'
#' Fit an AR(p) model
#'
#' Fit  model with optional AIC selection of optimal lag length
#' **REMARK**: Subsample is built using kmax lags.
#'
#' @param y time series data frame (data contained in columns)
#' @param k integer. Lag length. If negative, AIC is used to pick the lag length.
#' @param kmax integer. For AIC pick the largest possible lag.
#'
#' @return `ar_mod` object, with elements
#' \describe{
#'   \item{`A`}{matrix of AR coefficients.}
#'   \item{`AIC`}{AIC values for various lags.}
#'   \item{`k`}{optimal lag, if should be determined by AIC, else k like given}
#'   \item{`res`}{matrix of residuals}
#'   \item{`Sigma`}{estimated variance matrix of residuals}
#'   \item{`vcov`}{estimation variance matrix, corresponds to row-wise vectorization}
#' }
#'
#'
#' @export
#'
fit_ar <- function(y, k = -1, kmax = 10) {

  y_mat <- y
  if (inherits(y, "tbl_ts")){
    y_mat <- t(tsibble_to_mat(y))
    # get names of vars 
    key_data <- tsibble::key(y)
    if (length(key_data)==0){
      # univariate data without name 
      ser[[1]] <- "y"
    }
    if (length(key_data)==1){
      ser <- as.data.frame(unique(y[,key_data[[1]]]))[[1]]
    }
  } else {
    ser <- rownames(y_mat)
  }
  if (inherits(y_mat,"matrix") == FALSE){
    stop('Input to fit_ar needs to be either a matrix or an appropriate data structure.')
  }
  
  Teff <- dim(y_mat)[2] - kmax  # Effective number of time points after accounting for maximum lag
  s <- dim(y_mat)[1]  # Number of rows (observations) in y
  AIC <- rep(0, kmax + 1)  # Initialize AIC vector of length kmax+1 with zeros

  mod <- list(s)
  if (k==0) {
    warning("Lag length k should not be zero, else coefficients are zero.")
    A <- matrix(0, s, 0)
    return(list(A = A, AIC = AIC, k = k))  # Return the AR coefficients, AIC values, and given zero lag
  }

  # If k is negative, use AIC to determine the optimal lag length
  if (k < 0) {
    Yk <- mlag(y_mat, kmax)  # Generate matrix with kmax lags of y
    yeff <- matrix(Yk[1:s, ], ncol = Teff)  # Extract effective response matrix
    res <- matrix(Yk[1:s, ], ncol = Teff)  # Initialize residuals matrix
    Sigma <- res %*% t(res) / Teff  # Calculate initial covariance matrix of residuals
    AIC[1] <- log(det(Sigma))  # Calculate AIC for zero lag

    # Loop through all possible lags up to kmax
    for (k in 1:kmax) {
      X <- matrix(Yk[s + c(1:k * s), ], ncol = Teff)  # Generate lagged predictor matrix

      XX <- X %*% t(X)  # Calculate cross-product of predictor matrix
      A <- yeff %*% t(X) %*% solve(XX)  # Calculate AR coefficients using least squares

      res <- yeff - A %*% X  # Calculate residuals
      Sigma <- res %*% t(res) / Teff  # Update covariance matrix of residuals

      AIC[k + 1] <- log(det(Sigma)) + 2 * k * s ^ 2 / Teff  # Calculate AIC for current lag
    }

    k <- which(AIC == min(AIC)) - 1  # Determine the optimal lag length by finding minimum AIC

    # Refit model using optimal lag length if greater than zero
    if (k > 0) {
      Yk <- mlag(y_mat, k)  # Generate matrix with optimal number of lags
      Teff <- dim(y_mat)[2] - k  # Update effective number of time points
      X <- matrix(Yk[s + c(1:(k * s)), ], ncol = Teff)  # Generate lagged predictor matrix
      A <- Yk[1:s, ] %*% t(X) %*% solve(X %*% t(X))  # Calculate AR coefficients
      res <- Yk[1:s, ] - A %*% X  # Calculate residuals
      Sigma <- res %*% t(res) / Teff  # Update covariance matrix of residuals
    } else {
      A <- matrix(0, s, 0)  # If optimal lag is zero, set coefficients to zero
      res <- yeff  # Calculate residuals
      Sigma <- res %*% t(res) / Teff  # Update covariance matrix of residuals
    }
  } else {  # If a specific lag length is provided
    Yk <- mlag(y_mat, k)  # Generate matrix with specified number of lags
    Teff <- dim(y_mat)[2] - k  # Update effective number of time points
    X <- matrix(Yk[s + c(1:(k * s)), ], ncol = Teff)  # Generate lagged predictor matrix
    yeff <- matrix(Yk[1:s, ], ncol = Teff)  # Extract effective response matrix
    A <- matrix(Yk[1:s, ], ncol = Teff) %*% t(X) %*% solve(X %*% t(X))  # Calculate AR coefficients
    res <- Yk[1:s, ] - A %*% X  # Calculate residuals
    Sigma <- res %*% t(res) / Teff  # Calculate covariance matrix of residuals
    AIC[k + 1] <- Teff*(s*log(2*pi)+log(det(Sigma))+s) + 2 * k * s ^ 2    # Calculate AIC for provided lag length

  }

  # calculate variance 
  if (k>0){
    vcov <- kronecker(Sigma,solve(X %*% t(X))) # variance corresponds to row-wise vectorization!!
    rownames(A) <- ser
    labs <- list(k)
    for (j in 1:k){
      labs[[j]] <- paste('L',j,sep='')
    }
    col_names <- as.vector(outer(
      ser,
      labs,
      function(c, r) {paste(r, c, sep = ':')}
    ))
    
    colnames(A) <- col_names
    
    # generate output from standard function 'lm' 
    
    for ( j in 1:s){
      data <- data.frame( y = Yk[j,], X = t(X))
      colnames(data) <- c(as.character(ser[j]), colnames(A))
      mod[[j]] <- eval(parse(text= paste('lm(',as.character(ser[j]), '~ . - 1,data =data)')))
    }
    
    
  } else { 
    vv <- matrix(0,0,0)
  }
  rownames(res) <- ser

  ar_mod <- list(A = A, AIC = AIC, k = k,res= res, Sigma = Sigma, vcov = vcov, mods = mod)
  class(ar_mod) <- "ar_mod"
  return(ar_mod)  # Return the AR coefficients, AIC values, and optimal lag
}

#' Summary function provides a summary for an estimated AR(k) model.
#' 
#' @param object ... AR(k) model. 
#' @param ... additional arguments passed on. 
#' #' 
#' @return output on screen
#'
#' @export
#' @method summary ar_mod
#' 
summary.ar_mod <- function(object,...){
  
  Teff <- dim(object$res)[2]
  s <- dim(object$res)[1]
  ks <- dim(object$A)[2]
  k <- ks / s 
  
  v_beta <- diag(object$vcov)
  s_beta <- sqrt(v_beta)
  
  std_beta <- t(matrix(s_beta,k*s,s))
  
  
  ### build output
  structure(
    list(
      "T" = Teff,
      "k" = k,
      "vars" = rownames(object$A),
      "A" = object$A,
      "std" = std_beta, 
      "AIC" = object$AIC[k+1],
      "Sigma" = object$Sigma, 
      "LL" = -0.5*Teff*(s*log(2*pi)+log(det(object$Sigma))+s),
      "mods" =  object$mods
    ),
    class = c("summary.ar_mod", "list")
  )    
}


#' @noRd
#' @exportS3Method

print.summary.ar_mod <- function(x, ...) {
  
  writeLines("Summary of AR model:")
  
  ### Data information
  writeLines("")
  writeLines(paste("Sample size: T =", x$T))
  writeLines(paste("Lag length: k=", x$k))
  writeLines(paste('Log likelihood: ',x$LL,', AIC: ',x$AIC))
  writeLines("")

  ### compile matrix with estimates and stds for printing 
  s <- dim(x$A)[1]
  est_sd <- cbind(x$A,x$std)
  ### estimation info 
  print_est_sd(
    est_sd,
    row_names = rownames(x$A),
    col_names = colnames(x$A)
  )
  writeLines("")
  
  ### add summary from regressions 
  for (j in 1:s){
    print(summary(x$mods[[j]]))
    writeLines("")
  }
  
}


#' Print Estimated and Standard Deviations (taken from Rprobit)
#'
#' @description
#' This helper function prints a formatted version of a matrix with estimates
#' and standard deviations.
#'
#' @param est_sd
#' A \code{matrix} of dimension \code{m} times \code{2n}, where the first
#' \code{m} columns are estimates and the
#' last \code{m} columns are corresponding standard deviations.
#' @param row_names
#' A \code{character} vector of length \code{m} with row names.
#' Can be \code{NULL} (default) for no row names.
#' @param col_names
#' A \code{character} vector of length \code{n} with column names.
#' Can be \code{NULL} (default) for no column names.
#'
#' @return
#' No return value, prints a \code{matrix}.
#'
#' @examples
#' \dontrun{
#' est_sd <- matrix(1:12, nrow = 3, ncol = 4)
#' print_est_sd(est_sd, row_names = LETTERS[1:3], col_names = LETTERS[4:5])
#' }
#'
#' @keywords internal utils

print_est_sd <- function(est_sd, row_names = NULL, col_names = NULL) {
  
  ### compute dimensions
  nr <- dim(est_sd)[1]
  nc <- dim(est_sd)[2] / 2
  
  ### format the rownames
  if (is.null(row_names)) {
    row_names_print <- rep("", nr)
  } else {
    row_names_print <- row_names
    max_length <- max(nchar(row_names_print))
    for (i in 1:length(row_names_print)) {
      while (nchar(row_names_print[i]) < max_length) {
        row_names_print[i] <- paste(row_names_print[i], " ", sep = "")
      }
    }
  }
  
  ### format the entries to have equal length
  est_sd_print <- est_sd
  for (r in 1:nr) {
    for (c in 1:(2 * nc)) {
      est_sd_print[r, c] <- sprintf("%f", est_sd[r, c])
    }
  }
  for (c in 1:(2 * nc)) {
    max_length <- max(nchar(est_sd_print[, c]))
    for (r in 1:nr) {
      while (nchar(est_sd_print[r, c]) < max_length) {
        est_sd_print[r, c] <- paste(" ", est_sd_print[r, c], sep = "")
      }
    }
  }
  
  ### add column names
  if (!is.null(col_names)) {
    col_names_print <- as.character(col_names)
    str <- paste(rep(" ", nchar(row_names_print[1]) + 1), sep = "", collapse = "")
    for (c in 1:nc) {
      while (nchar(col_names_print[c]) < (nchar(est_sd_print[1, c]) + nchar(est_sd_print[1, c + nc]) + 4)) {
        col_names_print[c] <- paste(col_names_print[c], " ", sep = "")
      }
      if (nchar(col_names_print[c]) > (nchar(est_sd_print[1, c]) + nchar(est_sd_print[1, c + nc]) + 5)) {
        col_names_print[c] <- paste(substring(col_names_print[c], 1, (nchar(est_sd_print[1, c]) + nchar(est_sd_print[1, c + nc]) + 2)), ". ", sep = "")
      }
      
      str <- c(str, col_names_print[c])
    }
    cat(c(str, "\n"))
  }
  
  ### print the matrix row by row
  for (r in 1:nr) {
    str <- row_names_print[r]
    for (c in 1:nc) {
      str <- c(str, sprintf(" %s (%s)", est_sd_print[r, c], est_sd_print[r, c + nc]))
    }
    cat(c(str, "\n"))
  }
}

#'
#' Simulate autoregressive time series via companion form of AR parameters
#'
#' Receiving the parameters of an AR(p) process, this function creates the according companion matrix for it, and simulates the system.
#'
#' @param A s x sk matrix of coefficients as estimated in fit_ar. 
#' @param Sigma sxs covariance matrix of noise. Default: I(s). 
#' @param T integer; sample size 
#' @param stationary Boolean indicating whether initialisation should be done with a stationary distribution. 
#'
#' @return list.
#' \describe{
#'   \item{`ysim`}{simulated time series observations in df format}
#' }
#'
#' @export
#'
sim_ar <- function(A,Sigma,T=100,stationary = TRUE){
  
  # convert AR(k) to StSp system 
  syst <- ar_to_StSp(A,Sigma)

  # simulate state space system 
  y_sim <- sim_syst(syst,T= T,stationary=stationary)
  
  return(y_sim)
}


#'
#' Convert an AR(k) into a state space system using the companion form. 
#'
#' @param A s x sk matrix of coefficients as estimated in fit_ar. 
#' @param Sigma sxs covariance matrix of noise. Default: I(s). 
#'
#' @return syst StSp object
#'
#' @export
#'
ar_to_StSp <- function(A,Sigma){
  s = dim(A)[1]
  sk = dim(A)[2]
  k = sk/s # lag length equals k. 
  
  if (is.null(Sigma)){ Sigma = diag(s)}
  
  # transform AR to state space 
  # companion matrix A with coefficents in first row and diagonal matrix in following rows
  Pi <- matrix(0, nrow = sk, ncol = sk)
  if (k>1){
    Pi[1:s,] <- A
    Pi[(s+1):sk, 1:(sk-s)] <- diag(sk-s)
  } else {
    Pi <- A
  }
  
  # matrix B as row vector with one in first row and zeros in following rows
  K <- matrix(0, nrow = sk, ncol = s)
  K[1:s, 1:s] <- diag(s) 
  S <- K %*% Sigma 
  Q <- S %*% t(K)
  R <- Sigma
  syst <- StSp(A=Pi,C= A,K = K, Sigma = Sigma, Q=Q,R=R,S=S, type = "innov")
  
  return (syst)
} 

#'
#' check for stability of an AR system (via the eigenvalues of the companion form)
#'
#' @param object ar_model 
#'
#' @return stab_cond list containing the entries 'check' Boolean (TRUE, if stable) and 'poles' (eigenvalues of companion matrix)
#'
#' @export
check_stability <- function(object)
{
  #print("Calling the base function")
  UseMethod("check_stability")
  print(" this is not executed")
}


#'
#' check for stability of an AR system (via the eigenvalues of the companion form)
#'
#' @param A s x sk matrix of coefficients as estimated in fit_ar. 
#'
#' @return stab_cond list containing the entries 'check' Boolean (TRUE, if stable) and 'poles' (eigenvalues of companion matrix)
#'
#' @export
#'
check_stability_ar <- function(A){  
  
  s <- dim(A)[1]
  Sigma <- diag(s)
  syst <- ar_to_StSp(A,Sigma)
  
  stab_cond <- check_stability.StSp(syst)
  return( stab_cond ) 
}
  

