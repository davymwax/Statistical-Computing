#' OLS Beta Estimates using Cholesky Decomposition
#'
#' @param X: An n by p design matrix (including intercept)
#' @param y: A numeric vector of length n
#' @return A numeric vector
#' 
#' @export
#'
#' @examples
#' 
#' X <- matrix(c(4, 0, 0, 0, 3, 7, 0, 0, 7, 8, 5, 0, 9, 2, 11, 15), byrow=T, nrow=4)
#' y <- c(7, 18, 4, 11)
#' chol.beta.estimator(X, y)
#' 
chol.beta.estimator <- function(X, y){
  A = t(X) %*% X
  if (any(eigen(A)$values <= 0)){
    warning("Error: Sigma must be positive definite")} else {
  L = chol(A)
  inv_A = solve(L %*% t(L))
  beta_hat = matvecprod(inv_A, t(X), y, W=F)
  return(beta_hat)
    }
}
