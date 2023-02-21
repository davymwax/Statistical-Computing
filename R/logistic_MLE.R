#' Logistic MLE
#'
#' @param Y: n times 1 vector of binary responses
#' @param X: n times p-1 design matrix (these are your k = (p - 1) predictors; function will add intercept term)
#' @param init: vector p of initial values (include initial value for intercept term)
#' @param method: An integer 1 for Newton-Raphson (Fisher Scoring), 2 for Gradient Descent
#' @param iter: An integer for how many iterations are you willing to do
#'
#' @return a list of MLE estimates, their corresponding asymptotic confidence intervals, and vector of log-likelihoods
#' @export
#'
#' @examples Y <- c(1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
#' X <- matrix(data = c(3.89, 3.80, 3.54, 3.63, 3.16, 
#' 3.50, 3.34, 3.02, 2.87, 3.38, 3.56, 2.91), 
#' ncol=1, byrow=F)
#' iter <- 200
#' init <- c(0, 0)
#' logistic_MLE(Y, X, init, iter, method = 1)
#' logistic_MLE(Y, X, init, iter, method = 2)

logistic_MLE <- function(Y, X, init, iter, method){
  if (is.vector(Y) == F)
  {warning("Wrong class of object for Y")} else {
    if (is.matrix(X) == F | (length(Y) != dim(X)[1]))
    {warning("Wrong class of object for X. Check the dimensions of your objects")}}
  if (missing(init) == T | length(init) != (dim(X)[2] + 1))
  {warning("Provide suitable values for the optimization. Check dimensions of your objects")} else {
    if (missing(method) == T | method !=1 & method != 2) 
    {warning("Choose a suitable method. 1 for Fisher Scoring or 2 for Gradient Descent")}}
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
  output <- list(coeffs, CIs, log_lklhd)
  names(output) <- c("MLE Coefficients", "Asymptotic 95% C.I.s", "Log-Likelihood Values")
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
  output<-list(coeffs, CIs, log_lklhd)
  names(output)<-c("MLE Coefficients", "Asymptotic 95% C.I.s", "Log-Likelihood Values")
}
  return(output)
}







