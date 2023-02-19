#' Logistic MLE
#'
#' @param Y: n times 1 matrix of binary responses
#' @param X: n times p design matrix
#' @param init: vector p of initial values
#' @param method: 1 for Newton-Raphson (Fisher Scoring), 2 for Gradient Descent
#'
#' @return a vector of MLE estimates of a logistic regression model
#' @export
#'
#' @examples Y <- c(1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
#' X <- matrix(data = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3.89, 3.80, 3.54, 3.63, 3.16, 3.50, 3.34, 3.02, 2.87, 3.38, 3.56, 2.91), ncol=2, byrow=F)
#' init <- c(0, 0.95)
logistic_MLE <- function(Y, X, init, method)
{
  if (method == 1){
  n = length(Y)
  p = dim(X)[2]
  betas <- matrix(data=NA, nrow = p, ncol=10)
  g <- function(x)
  {
    p <- exp(x)/(1 + exp(x))
    return(p)
  }
  prob <- matrix(data=NA, nrow = n, ncol=10)
  b <- matrix(data=NA, nrow=n, ncol=10)
  I <- c(rep(1, n))
  log_lklhd <- matrix(data=NA, nrow=1, ncol=10)
  theta <- matrix(data = NA, nrow = 12, ncol=10)
  weights <- list()
  W <- list()
  info_number <- list()
  score_func <- list()
  output <- list()
  betas[, 1] <- init
  theta[,1] <- g(sum(betas[, 1]))
  prob[,1] <- theta[,1]
  b[, 1] <- -log(1 - prob[ ,1])
  weights[[1]] <- prob[, 1]*(1-prob[, 1])
  W[[1]] <- diag(weights[[1]])
  log_lklhd[, 1] <- t(Y) %*% (X %*% betas[, 1]) - (t(b[, 1])%*%I)
  score_func[[1]] <- t(X) %*% (Y - prob[, 1])
  info_number[[1]] <- (t(X) %*% W[[1]]) %*% X
  betas[, 2] <- betas[, 1] + (solve(info_number[[1]]) %*% score_func[[1]])
  for (i in 2:10){
  for (j in 1:n){
  theta[j, i] <- X[j, ] %*% betas[, i]
  prob[j, i] <- g(theta[j, i])
  b[j, i] <- -log(1 - prob[j , i])
  }
  weights[[i]] <- prob[, i]*(1-prob[, i])
  W[[i]] <- diag(weights[[i]])
  log_lklhd[, i] <- t(Y) %*% (X %*% betas[, i]) - (t(b[, i])%*%I)
  score_func[[i]] <- t(X) %*% (Y - prob[, i])
  info_number[[i]] <- (t(X) %*% W[[i]]) %*% X
  betas[, i+1] <- betas[, i] + (solve(info_number[[i]]) %*% score_func[[i]])
  }
  output[[1]] <- betas
  output[[2]] <- info_number
  output[[3]] <- log_lklhd
}
  return(output)
}




