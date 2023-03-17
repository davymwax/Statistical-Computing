#' Rosenbluth-Hastings Gamma MCMC
#'
#' @param curr_val This is the current state of the R-H Markov chain
#' @param tuning_param This is a tuning parameter of the Normal random walk
#'
#' @return Numeric vector of length 2. The first element of the vector contains the next state of the M-H Markov chain. The second element contains 0 if the proposed values was rejected and 1 otherwise. 
#' @export 
#'
#' @examples
#' lnorm_rw_next(1, 1.5)
lnorm_rw_next = function (curr_val, tuning_param){
  return_val = c(curr_val, 0)
  prop_val = curr_val*exp(rnorm(1, 0, tuning_param))
  q_prop_val = dlnorm(prop_val, exp(0), exp(tuning_param))
  q_curr_val = dlnorm(curr_val, exp(0), exp(tuning_param))
  pi_prop_val = dgamma(prop_val, shape =4, scale=2)
  pi_curr_val = dgamma(curr_val, shape =4, scale=2)
  if (runif(1) < (pi_prop_val*q_prop_val)/(pi_curr_val*q_curr_val)){
    return_val = c(prop_val, 1)
  }
  return(return_val)
}






