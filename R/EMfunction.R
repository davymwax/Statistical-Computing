rm(list=ls())
test_data <- read.table("~/Documents/09_R/07_Winter_23/STATS230/firstpkg/data/mixture_data.txt", 
                        header = F)
y <- as.matrix(test_data)
theta <- c(0.5, 4, 8)
iter = 200
#' MLE of Normal Mixture Models using the E-M Algorithm
#'
#' @param y: a vector of observed data
#' @param theta: a vector of initial values for all the parameters
#' @param iter: an integer for the number of iterations of the optimization task
#'
#' @return a vector of the Maximum Likelihood Estimates of your parameters
#' @export
#'
#' @examples rm(list=ls())
#' test_data <- read.table("mixture_data.txt", header = F)
#' y <- as.matrix(test_data)
#' theta = c(0.5, 4, 8)
#' iter = 200
#' EMfunction(y, theta, iter)
#' 
EMfunction <- function (y, theta, iter){
  #storing the output
  n=dim(y)[1]
  alpha <- matrix(data=NA, nrow=1, ncol=2*iter)
  mu1 <- matrix(data=NA, nrow=1, ncol=2*iter)
  mu2 <- matrix(data=NA, nrow=1, ncol=2*iter)
  N_1 <- list() 
  N_2 <- list()
  P_1 <- list()
  P_2 <- list()
  theta_MLE <- list()
  
  #checking convergence
  epsilon <-  epsilon <- function(a, b){
    sqrt(sum((a - b)^2))
  }
 
  #count number of steps
  steps <- 1
  #Initialization
  alpha[, 1] <- theta[1]
  mu1[, 1] <- theta[2]
  mu2[, 1] <- theta[3]
  theta_MLE[[1]] <- c(alpha[, 1], mu1[, 1], mu2[, 1])
  
  #E-step
  N_1[[1]] <- 1/sqrt(2*pi)*exp(-(y - mu1[, 1])^2/2)
  N_2[[1]] <- 1/sqrt(1.6*pi)*exp(-(y - mu2[, 1])^2/1.6)
  P_1[[1]] <- (alpha[, 1] * N_1[[1]])/(alpha[, 1]*N_1[[1]] +  (1 - alpha[, 1])*N_2[[1]])
  P_2[[1]] <- ((1 - alpha[, 1])*N_2[[1]])/(alpha[,1]*N_1[[1]] +  (1 - alpha[,1])*N_2[[1]])
  
  #M-step
  alpha[, 2] <- sum(P_1[[1]])/(sum(P_1[[1]]) + sum(P_2[[1]]))
  mu1[, 2] <- sum(P_1[[1]]*y)/sum(P_1[[1]])
  mu2[, 2] <- sum(P_2[[1]]*y)/sum(P_2[[1]])
  theta_MLE[[2]] <- c(alpha[, 2], mu1[, 2], mu2[,2])
 
 #iteration part
  for (i in 2:iter){
      #E-step
    N_1[[i]] <- 1/sqrt(2*pi)*exp(-(y - mu1[, i])^2/2)
    N_2[[i]] <- 1/sqrt(1.6*pi)*exp(-(y - mu2[, i])^2/1.6)
    P_1[[i]] <- (alpha[, i]*N_1[[i]])/(alpha[, i]*N_1[[i]]+(1 - alpha[, i])*N_2[[i]])
    P_2[[i]] <- ((1-alpha[, 1])*N_2[[i]])/(alpha[,i]*N_1[[i]]+(1 - alpha[,i])*N_2[[i]])
      
      #M-step
    alpha[, i + 1] <- sum(P_1[[i]])/(sum(P_1[[i]]) + sum(P_2[[i]]))
    mu1[, i + 1] <- sum(P_1[[i]]*y)/sum(P_1[[i]])
    mu2[, i + 1] <- sum(P_2[[i]]*y)/sum(P_2[[i]])
    theta_MLE[[i + 1]] <- c(alpha[, i+1], mu1[, i+1], mu2[,i+1])
    
     #iterate
    steps <- steps + 1
    
     #stopping rule
    if (epsilon(theta_MLE[[i + 1]], theta_MLE[[i]]) < 0.001)
        {break}
    if (steps >= iter){
      print(paste("Failed to converge after", steps, "steps"))
      break
      }
  }
return(theta_MLE[[steps]])
}


























