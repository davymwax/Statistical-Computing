#' Logistic MLE
#'
#' @param Y: n times 1 matrix of binary responses
#' @param X: n times p-1 design matrix (these are your k = (p - 1) predictors; function will add intercept term)
#' @param init: vector p of initial values (include initial value for intercept term)
#' @param method: An integer 1 for Newton-Raphson (Fisher Scoring), 2 for Gradient Descent
#' @param iter: how many iterations are you willing to do
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
#' iter <- 200
#' 
#' model <- glm(chd~sbp+ldl+alcohol,family = binomial(link="logit"), data=mydata)
#' model$coefficients
#'
logistic_MLE <- function(Y, X, init, iter, method){
  if (method == 1){
  n = length(Y)
  intercept <- as.numeric(c(rep(1, n)))
  X <- cbind(intercept, X)
  p = dim(X)[2]
  g <- function(x)
    {
    p_i <- exp(x)/(1 + exp(x))
    return(p_i)
  }
  epsilon <- function(a, b){
    sqrt(sum((a - b)^2))
  }
  betas <- matrix(data=NA, nrow = p, ncol=2*iter)
  thetas <- matrix(data=NA, nrow = n, ncol=iter)
  prob <- matrix(data=NA, nrow = n, ncol=iter)
  b <- matrix(data=NA, nrow=n, ncol=iter)
  I <- c(rep(1, n))
  log_lklhd <- matrix(data=NA, nrow=1, ncol=iter)
  weights <- list()
  W <- list()
  info_number <- list()
  as_variance <- list()
  score_func <- list()
  betas[, 1] <- init
  thetas[,1] <- X %*% betas[, 1]
  prob[, 1] <- g(thetas[,1])
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
  steps <- 0
  for (i in 2:iter){
    thetas[,i] <- X %*% betas[, i]
    prob[, i] <- g(thetas[,i])
    b[, i] <- -log(1 - prob[, i])
    weights[[i]] <- prob[, i]*(1 - prob[, i])
    W[[i]] <- diag(weights[[i]])
    log_lklhd[, i] <- t(Y) %*% (X %*% betas[, i]) - (t(b[, i])%*%I)
    score_func[[i]] <- t(X) %*% (Y - prob[, i])
    info_number[[i]] <- (t(X) %*% W[[i]]) %*% X
    as_variance[[i]] <- solve(info_number[[i]])
    betas[, i + 1] <- betas[, i] + (as_variance[[i]] %*% score_func[[i]])
    steps <- steps + 1
    if (epsilon(betas[, i + 1], betas[, i]) < 0.001)
    {break}
    if (steps > iter){
        print(paste("Failed to converge after", steps, "steps"))
        break
      }
    }
  coeffs <- betas[, steps+1]
  variance <- diag(as_variance[[steps+1]])
  se <- sqrt(variance/n)
  CIs <- cbind(LB=coeffs - 1.96*se, UB=coeffs + 1.96*se)
  log_lklhd <- log_lklhd[1:(steps+1)]
  }
  if (method == 2){
  n = length(Y)
  alpha = 0.05
  intercept <- as.numeric(c(rep(1, n)))
  X <- cbind(intercept, X)
  p = dim(X)[2]
  g <- function(x)
  {
    p_i <- exp(x)/(1 + exp(x))
    return(p_i)
  }
  epsilon <- function(a, b){
    sqrt(sum((a - b)^2))
  }
  betas <- matrix(data=NA, nrow = p, ncol=2*iter)
  thetas <- matrix(data=NA, nrow = n, ncol=iter)
  prob <- matrix(data=NA, nrow = n, ncol=iter)
  b <- matrix(data=NA, nrow=n, ncol=iter)
  I <- c(rep(1, n))
  log_lklhd <- matrix(data=NA, nrow=1, ncol=iter)
  weights <- list()
  W <- list()
  info_number <- list()
  as_variance <- list()
  score_func <- list()
  betas[, 1] <- init
  thetas[,1] <- X %*% betas[, 1]
  prob[, 1] <- g(thetas[,1])
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
  betas[, 2] <- betas[, 1] + ((alpha*as_variance[[1]]) %*% score_func[[1]])
  steps <- 0
  for (i in 2:iter){
    thetas[,i] <- X %*% betas[, i]
    prob[, i] <- g(thetas[,i])
    b[, i] <- -log(1 - prob[, i])
    weights[[i]] <- prob[, i]*(1 - prob[, i])
    W[[i]] <- diag(weights[[i]])
    log_lklhd[, i] <- t(Y) %*% (X %*% betas[, i]) - (t(b[, i])%*%I)
    score_func[[i]] <- t(X) %*% (Y - prob[, i])
    info_number[[i]] <- (t(X) %*% W[[i]]) %*% X
    as_variance[[i]] <- solve(info_number[[i]])
    betas[, i + 1] <- betas[, i] + (alpha*as_variance[[i]] %*% score_func[[i]])
    steps <- steps + 1
    if (epsilon(betas[, i + 1], betas[, i]) < 0.001)
    {break}
    if (steps >= iter){
      print(paste("Failed to converge after", steps, "steps"))
      break
    }
  }
  coeffs <- betas[, steps+1]
  variance <- diag(as_variance[[steps+1]])
  se <- sqrt(variance/n)
  CIs <- cbind(LB=coeffs - 1.96*se, UB=coeffs + 1.96*se)
  log_lklhd <- log_lklhd[1:(steps+1)]
}
  return(list(coeffs, CIs, log_lklhd))
}






