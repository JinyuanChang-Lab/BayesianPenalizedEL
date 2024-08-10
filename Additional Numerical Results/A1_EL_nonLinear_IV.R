library(MASS)
library(gmm)

############# EE for IV model #############
IV.fun <- function(beta, X) {
  y <- X[,1]
  x <- X[,2:3]
  z <- X[,4:7]
  
  residual <- y - sin(x %*% beta)
  residual <- as.vector(residual)
  residual <- diag(residual)
  
  G.ee<- residual %*% z
  
  return(G.ee)
}


############# rep.fun #############
rep.fun <- function(ind)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  beta_el_all <- matrix(0, 2, init_num)
  
  for (j in c(1:init_num)) {
    
    init <- init.data[, j]
    init <- as.matrix(init)
    # print(init)
    
    cat('---- init_', j, '----')
    
    cat("...........begin.EL............\n")     ##### EL
    
    res <- gel(IV.fun, data, init, type="EL")
    
    beta_el <- res$coefficients
    beta_el_all[, j] <- beta_el
    
    cat("...........end.EL............\n")
  
  }
  
  return(list(beta_el_all = beta_el_all))
}

