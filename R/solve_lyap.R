#'
#' solve the Lyapunov equation P = F P F' + Q. 
#'
#' @param A nxn transition matrix of the state space model.
#' @param Q nxn state noise variance
#'
#' @examples
#' # solution to Lyapunov equation 
#' P <- solve_lyap(A=as.matrix(0.5), Q=as.matrix(1))
#'
#' @return P nxn stationary state variance. 
#'
#' @export
#'
solve_lyap <- function(A=as.matrix(0.5), Q=as.matrix(1)){
  
  tol = 0.00001 # set the tolerance level
  
  # check stability of A
  ev <- eigen(A)$values
  if (max(abs(ev)) > 1-tol ) {
    stop("solve_lyap: F has unstable eigenvalues!")
  } 
  n <- dim(A)[1]
  Akron <- diag(n ^ 2) - kronecker(A, A)
  vP <- solve(Akron) %*% matrix(Q, ncol = 1)
  P <- matrix(vP, ncol = n)
  P <- (P + t(P)) / 2 # ensure that P is symmetric
    
  return(P)
}


