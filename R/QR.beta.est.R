#' OLS Estimates using QR Decomposition
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
#' QR.beta.estimator(X, y)
#' 
QR.beta.estimator <- function(X, y){
  QR <- qr(X)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  inv_R <- solve(R)
  beta_hat = matvecprod(inv_R, t(Q), y, W=F)
  return(beta_hat)
}

