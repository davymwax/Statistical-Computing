#' Logistic MLE
#'
#' @param Y: n times 1 matrix of binary responses
#' @param X: n times p-1 design matrix (these are your k = (p - 1) predictors; function will add intercept term)
#' @param init: vector p of initial values (include initial value for intercept term)
#' @param method: 1 for Newton-Raphson (Fisher Scoring), 2 for Gradient Descent
#'
#' @return a vector of MLE estimates of a logistic regression model
#' @export
#'
#' @examples 
#' rm(list = ls())
#' Y <- c(1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
#' X <- matrix(data = c(3.89, 3.80, 3.54, 3.63, 3.16, 3.50, 3.34, 3.02, 2.87, 3.38, 3.56, 2.91), ncol=1, byrow=F)
#' init <- c(0, 0.95)
#'
#'
#' rm(list = ls())
#' mydata <- read.csv("/Users/dmwakima/Documents/09_R/07_Winter_23/STATS230/firstpkg/data/SAheart.csv")
#' Y <- as.numeric(mydata[, c("chd")])
#' X <- as.matrix(mydata[, c("sbp", "ldl", "alcohol")])
#' init <- c(0, 0, 0, 0)
#' iter <- 100
#' model <- glm(chd~sbp+ldl+alcohol,family = binomial(link="logit"), data=mydata)
#' model$coefficients
#' 
logistic_MLE <- function(Y, X, init, method)
{
  if (method == 1){
  n = length(Y)
  intercept <- as.numeric(c(rep(1, n)))
  X <- cbind(intercept, X)
  p = dim(X)[2]
  g <- function(x)
    {
    p <- exp(x)/(1 + exp(x))
    return(p)
  }
  betas <- matrix(data=NA, nrow = p, ncol=2*iter)
  prob <- matrix(data=NA, nrow = n, ncol=iter)
  b <- matrix(data=NA, nrow=n, ncol=iter)
  I <- c(rep(1, n))
  log_lklhd <- matrix(data=NA, nrow=1, ncol=iter)
  weights <- list()
  W <- list()
  info_number <- list()
  as_variance <- list()
  score_func <- list()
  output <- list()
  betas[, 1] <- init
  prob[, 1] <- g(sum(betas[, 1]))
  b[, 1] <- -log(1 - prob[ ,1])
  weights[[1]] <- prob[, 1]*(1-prob[, 1])
  W[[1]] <- diag(weights[[1]])
  log_lklhd[, 1] <- t(Y) %*% (X %*% betas[, 1]) - (t(b[, 1])%*%I)
  score_func[[1]] <- t(X) %*% (Y - prob[, 1])
  info_number[[1]] <- (t(X) %*% W[[1]]) %*% X
  Sigma = diag(svd(info_number[[1]])$d)
  inv_Sigma = solve(Sigma)
  V = svd(info_number[[1]])$v
  U = svd(info_number[[1]])$u
  as_variance[[1]] <- V %*% (inv_Sigma %*% t(U))
  betas[, 2] <- betas[, 1] + (as_variance[[1]] %*% score_func[[1]])
  for (i in 2:iter){
      # no idea why log odds is not a function of X
      # for (j in 1:n){
      #   theta[j, i] <- X[j, ] %*% betas[, i]
      #   prob[j, i] <- g(theta[j, i])
      #   b[j, i] <- -log(1 - prob[j , i])
      #   }
    prob[, i] <- g(sum(betas[, i]))
    b[, i] <- -log(1 - prob[, i])
    weights[[i]] <- prob[, i]*(1 - prob[, i])
    W[[i]] <- diag(weights[[i]])
    log_lklhd[, i] <- t(Y) %*% (X %*% betas[, i]) - (t(b[, i])%*%I)
    score_func[[i]] <- t(X) %*% (Y - prob[, i])
    info_number[[i]] <- (t(X) %*% W[[i]]) %*% X
      # SVD calculation of the variance-covariance matrix failed.
      # Sigma = diag(svd(info_number[[i]])$d)
      # inv_Sigma = solve(Sigma)
      # V = svd(info_number[[i]])$v
      # U = svd(info_number[[i]])$u
      # use ginv() in MASS package
    as_variance[[i]] <- ginv(info_number[[i]])
    betas[, i + 1] <- betas[, i] + (as_variance[[i]] %*% score_func[[i]])
  }
  output[[1]] <- betas
  output[[2]] <- as_variance
  output[[3]] <- log_lklhd
  }
  if (method == 2){
    n = length(Y)
    intercept <- as.numeric(c(rep(1, n)))
    X <- cbind(intercept, X)
    p = dim(X)[2]
    g <- function(x)
    {
      p <- exp(x)/(1 + exp(x))
      return(p)
    }
    betas <- matrix(data=NA, nrow = p, ncol=2*iter)
    prob <- matrix(data=NA, nrow = n, ncol=iter)
    b <- matrix(data=NA, nrow=n, ncol=iter)
    I <- c(rep(1, n))
    log_lklhd <- matrix(data=NA, nrow=1, ncol=iter)
    score_func <- list()
    A <- list()
    output <- list()
    betas[, 1] <- init
    prob[, 1] <- g(sum(betas[, 1]))
    b[, 1] <- -log(1 - prob[ ,1])
    log_lklhd[, 1] <- t(Y) %*% (X %*% betas[, 1]) - (t(b[, 1])%*%I)
    score_func[[1]] <- t(X) %*% (Y - prob[, 1])
    A[[1]] <- 0.00005*diag(p)
    betas[, 2] <- betas[, 1] - (A[[1]] %*% score_func[[1]])
    for (i in 2:iter){
      # no idea why log odds is not a function of X
      # for (j in 1:n){
      #   theta[j, i] <- X[j, ] %*% betas[, i]
      #   prob[j, i] <- g(theta[j, i])
      #   b[j, i] <- -log(1 - prob[j , i])
      #   }
      prob[, i] <- g(sum(betas[, i]))
      b[, i] <- -log(1 - prob[, i])
      log_lklhd[, i] <- t(Y) %*% (X %*% betas[, i]) - (t(b[, i])%*%I)
      score_func[[i]] <- t(X) %*% (Y - prob[, i])
      A[[i]] <- 0.00005*diag(p)
      betas[, i + 1] <- betas[, i] - (A[[i]] %*% score_func[[i]])
    }
    output[[1]] <- betas
    output[[2]] <- as_variance
    output[[3]] <- log_lklhd
  }
  return(output)
}






