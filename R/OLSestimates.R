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
#' sim.data <- read.csv("https://github.com/miningroup/your-stat-230-r-package-davymwax/blob/main/homework2_regression.csv")
#' source("https://github.com/miningroup/your-stat-230-r-package-davymwax/blob/main/R/matvecprod.R")
#' dim(sim.data)
#' head(sim.data)
#' X <- as.matrix(sim.data[, c(3:6)])
#' y <- as.vector(sim.data[, 1])
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
  return(beta_hat)
    }
  } else {
    
  }
}

