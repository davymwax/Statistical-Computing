#' Simulating a Multivariate Normal
#' 
#' This package simulates N realizations multivariate normal given 
#'
#' @param mu : a numeric vector
#' @param sigma : an n by n positive definite matrix
#' @param N : numeric integer greater than or equal 1.
#'
#' @return a list
#' @export
#'
#' @examples
#' 
#' set.seed(1247)
#' N = 100
#' n = 4
#' mu = rep(0, n)
#' X = matrix(data=sample(1:100, replace = F, n*n), nrow=n)
#' sigma = cov(X)
#' multvnorm(mu, sigma, N)
#' 
#' @comparison
#' Comparing the results with the implementation in the MASS package
#' 
#' @examples
#' library(MASS)
#' mean(sapply(multvnorm(mu, sigma, N), mean))
#' mean(mvrnorm(n = N, mu, sigma))
#' 
#' sigma
#' var(mvrnorm(n = N, mu, sigma))
multvnorm <- function(mu, sigma, N){
  n = length(mu)
  z = vector(mode = "list", length = N)
  L = chol(sigma)
  x = vector(mode = "list", length = N)
  if (is.vector(mu) == F){
    warning("Error: The First Argument must be a vector")
  } else {
    if (is.matrix(sigma) == F){
      warning("Error: The Second Argument must be a matrix")
    } else {
      if (any(eigen(sigma)$values <= 0)){
        warning("Error: Sigma must be positive definite")
      } else {
        for (i in 1:N){
          z[[i]] <- rnorm(n, mean=0, sd=1)
        }
        for (j in 1:4){
          for (i in 1:N)
            x[[i]][j] <- sum(L[j, ] %*% z[[i]], mu[j])
        }
    }
    }
  }
  return(x)
}

