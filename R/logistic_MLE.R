#' Logistic MLE
#'
#' @param Y: n times 1 matrix of binary responses
#' @param X: n times p design matrix
#' @param method 
#'
#' @return a vector of MLE estimates of a logistic regression model
#' @export
#'
#' @examples
logistic_MLE <- function(Y, X, method)
{
  prob_Y <- mean(Y)
  log_odds <- function(prob_Y){
    prob_Y/(1 - prob_Y)
  }
  log_odds <- tr(X)*B
  while(distance.l2(B[i + 1] - B[i]<0.01)){
    
  }
}