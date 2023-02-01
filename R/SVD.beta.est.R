#' OLS Estimates using SVD Decomposition
#'
#' @param X: A n by p design matrix (including intercept)
#' @param y: A numeric vector of length n
#'
#' @return A numeric vector
#' @export
#'
#' @examples
#' 
#' X <- matrix(c(16, 12, 28, 36, 12, 58, 77, 41, 28, 77, 138, 134, 36, 41, 134, 431), byrow =T, nrow=4)
#' y <- c(7, 18, 4, 11)
#' SVD.beta.est(X, y)
#' 
SVD.beta.est <- function(X, y){
  n <- length(y)
  SVD <- svd(X)
  sigma <- diag(SVD$d)
  U <- SVD$u
  V <- SVD$v
  inv_sigma <- solve(sigma)
  # right_hand = matvecprod(inv_sigma, t(U) , y, W=F)
  # beta_hat = matvecprod(A = as.matrix(diag(n)), B=as.matrix(V), v=as.vector(right_hand), W=F)
  beta_hat = V %*% (inv_sigma %*% (t(U) %*% y))
  return(beta_hat)
}
