#' OLS Beta Estimates using either Cholesky or QR Decomposition
#'
#' @param X: An n by p design matrix (including intercept)
#' @param y: A numeric vector of length n
#' @param W: A logical (boolean) object. If W = T, then OLS algorithm uses Cholesky, otherwise OLS algorithm uses QR 
#' @return A numeric vector
#' 
#' @export
#'
#' @examples
#' 
#' set.seed(1247)
#' X <- matrix(sample(1:20, 16, replace=T), ncol=4)
#' y <- sample(1:20, 4, replace=T)
#' W = T
#' beta.estimator(X, y, W)
#'
beta.estimator <- function(X, y, W){
  A = cov(X)
  if (W == T){
  if (any(eigen(A)$values <= 0)){
    warning("Error: Sigma must be positive definite")} else {
  L = chol(A)
  inv_A = solve(L %*% t(L))
  beta_hat = matvecprod(inv_A, t(X), y, W=F)
    }
  } else {
    QR <- qr(X)
    Q <- qr.Q(QR)
    R <- qr.R(QR)
    inv_R <- solve(R)
    beta_hat = matvecprod(inv_R, t(Q), y, W=F)
  }
  return(beta_hat)
}

