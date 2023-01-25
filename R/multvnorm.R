#' Simulating a Multivariate Normal
#' 
#' This package simulates N realizations multivariate normal given 
#'
#' @param mu : a numeric vector
#' @param sigma : an m by m positive (semi)definite matrix
#' @param N : numeric integer greater than or equal 1.
#'
#' @return a list
#' @export
#'
#' @examples
#' N = 100
#' m = 4
#' n = 4
#' mu = rep(0, n)
#' X = matrix(data=sample(1:20, m*n), nrow=m)
#' sigma = X %*% t(X)
#' multvnorm(mu, sigma, N)
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
      if (dim(sigma)[1] != dim(sigma)[2]){
        warning("Error: Sigma must be positive (semi)definite")
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

