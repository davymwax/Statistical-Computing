rm(list=ls())
test_data <- read.table("~/Documents/09_R/07_Winter_23/STATS230/firstpkg/data/mixture_data.txt", 
                        header = F)
y <- as.matrix(test_data)
init.w <- 0.5
init.mu <- c(4, 8)
iter = 200
EMfunction <- function (y, init.w, init.mu, iter){
  #storing the output
  n=dim(y)[1]
  alpha <- matrix(data=NA, nrow=1, ncol=2*iter)
  mu1 <- matrix(data=NA, nrow=1, ncol=2*iter)
  mu2 <- matrix(data=NA, nrow=1, ncol=2*iter)
  N_1 <- list() 
  N_2 <- list()
  P_1 <- list()
  P_2 <- list()
  theta <- list()
  
  #checking convergence
  epsilon <-  epsilon <- function(a, b){
    sqrt(sum((a - b)^2))
  }
 
  #count number of steps
 steps <- 1
 #Initialization
 alpha[, 1] <- init.w
 mu1[, 1] <- init.mu[1]
 mu2[, 1] <- init.mu[2]
 theta[[1]] <- c(alpha[, 1], mu1[, 1], mu2[, 1])
 N_1[[1]] <- 1/sqrt(2*pi)*exp(-(y - mu1[, 1])^2/2)
 N_2[[1]] <- 1/sqrt(1.6*pi)*exp(-(y - mu2[, 1])^2/1.6)
 P_1[[1]] <- (alpha[, 1] * N_1[[1]])/(alpha[, 1]*N_1[[1]] +  (1 - alpha[, 1])*N_2[[1]])
 P_2[[1]] <- ((1 - alpha[, 1])*N_2[[1]])/(alpha[,1]*N_1[[1]] +  (1 - alpha[,1])*N_2[[1]])
 alpha[, 2] <- sum(P_1[[1]])/(sum(P_1[[1]]) + sum(P_2[[1]]))
 mu1[, 2] <- sum(P_1[[1]]*y)/sum(P_1[[1]])
 mu2[, 2] <- sum(P_2[[1]]*y)/sum(P_2[[1]])
 theta[[2]] <- c(alpha[, 2], mu1[, 2], mu2[,2])
 
 #iteration part
  for (i in 2:iter){
      #M-step
    N_1[[i]] <- 1/sqrt(2*pi)*exp(-(y - mu1[, i])^2/2)
    N_2[[i]] <- 1/sqrt(1.6*pi)*exp(-(y - mu2[, i])^2/1.6)
    P_1[[i]] <- (alpha[, i] * N_1[[i]])/(alpha[, i]*N_1[[i]] +  (1 - alpha[, i])*N_2[[i]])
    P_2[[i]] <- ((1 - alpha[, 1])*N_2[[i]])/(alpha[,i]*N_1[[i]] +  (1 - alpha[,i])*N_2[[i]])
      
      #E-step
    alpha[, i + 1] <- sum(P_1[[i]])/(sum(P_1[[i]]) + sum(P_2[[i]]))
    mu1[, i + 1] <- sum(P_1[[i]]*y)/sum(P_1[[i]])
    mu2[, i + 1] <- sum(P_2[[i]]*y)/sum(P_2[[i]])
    theta[[i + 1]] <- c(alpha[, i+1], mu1[, i+1], mu2[,i+1])
    
    steps <- steps + 1
    if (epsilon(theta[[i + 1]], theta[[i]]) < 0.001)
        {break}
    if (steps >= iter){
      print(paste("Failed to converge after", steps, "steps"))
      break
      }
  }
return(theta[[steps]])
}


























