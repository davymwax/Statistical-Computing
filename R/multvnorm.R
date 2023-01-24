#' Simulating a Multivariate Normal
#' 
#' This package simulates N realizations  multivariate normal given 
#'
#' @param mu 
#' @param sigma 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
#' N = 100
#' mu = rep(0, n)
#' X = matrix(data=sample(1:20, n*n), ncol=n)
#' sigma = X %*% t(X)
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

