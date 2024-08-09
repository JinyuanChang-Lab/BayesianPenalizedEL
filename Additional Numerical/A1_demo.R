library(MASS)

############# rep.tuning_1 #############
rep.tuning_1 <- function(nu, init, y, x, z) #BIC
{
  r <- ncol(z)
  n <- nrow(z)
  
  res_MH <- MH_1_arma(init, y, x, z, nu)
  beta_3_MH <- res_MH$beta_mean_3
  
  Gee_3_MH <- IVfun_arma(beta_3_MH, y, x, z)
  lambda_3_MH <- optimal_lam(Gee_3_MH, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MH <- sum(abs(lambda_3_MH) > 1e-6)
  BIC_3_MH <- log(sum((colMeans(Gee_3_MH))^2) / r) + Rn_3_MH * log(n) / n
  
  return(list(beta_3_MH = beta_3_MH, 
              BIC_3_MH = BIC_3_MH))
}


############# rep.fun_1 #############
rep.fun_1 <- function(ind, nu.seq)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  tun.num <-length(nu.seq)
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  BEL_res_3_MH <- matrix(0, init_num, tun.num*p)
  BIC_res_3_MH <- matrix(0, init_num, tun.num)
  MSE_res_3_MH <- matrix(0, init_num, tun.num)
  
  
  for (i in c(1:init_num)) {
    cat('---- init_', i, '----')
    
    for (j in c(1:tun.num)){
      
      init <- init.data[, j]
      init <- as.matrix(init)
      nu <- nu.seq[j]
      cat("...........begin.nu_",nu,"............\n")
      
      res <- rep.tuning_1(nu, init, y, x, z)
      
      BEL_res_3_MH[i, (j*p-1):(j*p)] <- res$beta_3_MH
      BIC_res_3_MH[i, j] <- res$BIC_3_MH
      MSE_res_3_MH[i,j] <- sum((res$beta_3_MH - c(0.5, 0.5))^2)
    }
  }
  
  optimal_num_3_MH <- apply(BIC_res_3_MH,1, which.min)
  ind_3_MH_1 <- cbind(1:init_num, optimal_num_3_MH*p-1)
  ind_3_MH_2 <- cbind(1:init_num, optimal_num_3_MH*p)
  ind_3_MH_3 <- cbind(1:init_num, optimal_num_3_MH)
  BEL_3_MH <- cbind(BEL_res_3_MH[ind_3_MH_1], BEL_res_3_MH[ind_3_MH_2])
  BIC_3_MH <- BIC_res_3_MH[ind_3_MH_3]
  MSE_3_MH <- MSE_res_3_MH[ind_3_MH_3]
  
  return(list(BEL_3_MH = BEL_3_MH, BIC_3_MH = BIC_3_MH, MSE_3_MH = MSE_3_MH))
}



############# rep.tuning_2 #############
rep.tuning_2 <- function(nu, init, y, x, z) #BIC
{
  r <- ncol(z)
  n <- nrow(z)
  
  res_MH <- MH_2_arma(init, y, x, z, nu)
  beta_3_MH <- res_MH$beta_mean_3
  
  Gee_3_MH <- IVfun_arma(beta_3_MH, y, x, z)
  lambda_3_MH <- optimal_lam(Gee_3_MH, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MH <- sum(abs(lambda_3_MH) > 1e-6)
  BIC_3_MH <- log(sum((colMeans(Gee_3_MH))^2) / r) + Rn_3_MH * log(n) / n
  
  return(list(beta_3_MH = beta_3_MH, 
              BIC_3_MH = BIC_3_MH))
}


############# rep.fun_2 #############
rep.fun_2 <- function(ind, nu.seq)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  tun.num <-length(nu.seq)
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  BEL_res_3_MH <- matrix(0, init_num, tun.num*p)
  BIC_res_3_MH <- matrix(0, init_num, tun.num)
  MSE_res_3_MH <- matrix(0, init_num, tun.num)
  
  
  for (i in c(1:init_num)) {
    cat('---- init_', i, '----')
    
    for (j in c(1:tun.num)){
      
      init <- init.data[, j]
      init <- as.matrix(init)
      nu <- nu.seq[j]
      cat("...........begin.nu_",nu,"............\n")
      
      res <- rep.tuning_2(nu, init, y, x, z)
      
      BEL_res_3_MH[i, (j*p-1):(j*p)] <- res$beta_3_MH
      BIC_res_3_MH[i, j] <- res$BIC_3_MH
      MSE_res_3_MH[i,j] <- sum((res$beta_3_MH - c(0.5, 0.5))^2)
    }
  }
  
  optimal_num_3_MH <- apply(BIC_res_3_MH,1, which.min)
  ind_3_MH_1 <- cbind(1:init_num, optimal_num_3_MH*p-1)
  ind_3_MH_2 <- cbind(1:init_num, optimal_num_3_MH*p)
  ind_3_MH_3 <- cbind(1:init_num, optimal_num_3_MH)
  BEL_3_MH <- cbind(BEL_res_3_MH[ind_3_MH_1], BEL_res_3_MH[ind_3_MH_2])
  BIC_3_MH <- BIC_res_3_MH[ind_3_MH_3]
  MSE_3_MH <- MSE_res_3_MH[ind_3_MH_3]
  
  return(list(BEL_3_MH = BEL_3_MH, BIC_3_MH = BIC_3_MH, MSE_3_MH = MSE_3_MH))
}


############# rep.tuning_3 #############
rep.tuning_3 <- function(nu, init, y, x, z) #BIC
{
  r <- ncol(z)
  n <- nrow(z)
  
  res_MH <- MH_3_arma(init, y, x, z, nu)
  beta_3_MH <- res_MH$beta_mean_3
  
  Gee_3_MH <- IVfun_arma(beta_3_MH, y, x, z)
  lambda_3_MH <- optimal_lam(Gee_3_MH, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MH <- sum(abs(lambda_3_MH) > 1e-6)
  BIC_3_MH <- log(sum((colMeans(Gee_3_MH))^2) / r) + Rn_3_MH * log(n) / n
  
  return(list(beta_3_MH = beta_3_MH, 
              BIC_3_MH = BIC_3_MH))
}


############# rep.fun_3 #############
rep.fun_3 <- function(ind, nu.seq)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  tun.num <-length(nu.seq)
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  BEL_res_3_MH <- matrix(0, init_num, tun.num*p)
  BIC_res_3_MH <- matrix(0, init_num, tun.num)
  MSE_res_3_MH <- matrix(0, init_num, tun.num)
  
  
  for (i in c(1:init_num)) {
    cat('---- init_', i, '----')
    
    for (j in c(1:tun.num)){
      
      init <- init.data[, j]
      init <- as.matrix(init)
      nu <- nu.seq[j]
      cat("...........begin.nu_",nu,"............\n")
      
      res <- rep.tuning_3(nu, init, y, x, z)
      
      BEL_res_3_MH[i, (j*p-1):(j*p)] <- res$beta_3_MH
      BIC_res_3_MH[i, j] <- res$BIC_3_MH
      MSE_res_3_MH[i,j] <- sum((res$beta_3_MH - c(0.5, 0.5))^2)
    }
  }
  
  optimal_num_3_MH <- apply(BIC_res_3_MH,1, which.min)
  ind_3_MH_1 <- cbind(1:init_num, optimal_num_3_MH*p-1)
  ind_3_MH_2 <- cbind(1:init_num, optimal_num_3_MH*p)
  ind_3_MH_3 <- cbind(1:init_num, optimal_num_3_MH)
  BEL_3_MH <- cbind(BEL_res_3_MH[ind_3_MH_1], BEL_res_3_MH[ind_3_MH_2])
  BIC_3_MH <- BIC_res_3_MH[ind_3_MH_3]
  MSE_3_MH <- MSE_res_3_MH[ind_3_MH_3]
  
  return(list(BEL_3_MH = BEL_3_MH, BIC_3_MH = BIC_3_MH, MSE_3_MH = MSE_3_MH))
}



############# rep.tuning_4 #############
rep.tuning_4 <- function(nu, init, y, x, z) #BIC
{
  r <- ncol(z)
  n <- nrow(z)
  
  res_MAMIS <- MAMIS_1_arma(init, y, x, z, nu)
  beta_3_MAMIS <- res_MAMIS$beta_3
  
  Gee_3_MAMIS <- IVfun_arma(beta_3_MAMIS, y, x, z)
  lambda_3_MAMIS <- optimal_lam(Gee_3_MAMIS, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MAMIS <- sum(abs(lambda_3_MAMIS) > 1e-6)
  BIC_3_MAMIS <- log(sum((colMeans(Gee_3_MAMIS))^2) / r) + Rn_3_MAMIS * log(n) / n
  
  return(list(beta_3_MAMIS = beta_3_MAMIS, 
              BIC_3_MAMIS = BIC_3_MAMIS))
}


############# rep.fun_4 #############
rep.fun_4 <- function(ind, nu.seq)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  tun.num <-length(nu.seq)
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  BEL_res_3_MAMIS <- matrix(0, init_num, tun.num*p)
  BIC_res_3_MAMIS <- matrix(0, init_num, tun.num)
  MSE_res_3_MAMIS <- matrix(0, init_num, tun.num)
  
  
  for (i in c(1:init_num)) {
    cat('---- init_', i, '----')
    
    for (j in c(1:tun.num)){
      
      init <- init.data[, j]
      init <- as.matrix(init)
      nu <- nu.seq[j]
      cat("...........begin.nu_",nu,"............\n")
      
      res <- rep.tuning_4(nu, init, y, x, z)
      
      BEL_res_3_MAMIS[i, (j*p-1):(j*p)] <- res$beta_3_MAMIS
      BIC_res_3_MAMIS[i, j] <- res$BIC_3_MAMIS
      MSE_res_3_MAMIS[i,j] <- sum((res$beta_3_MAMIS - c(0.5, 0.5))^2)
    }
  }
  
  optimal_num_3_MAMIS <- apply(BIC_res_3_MAMIS,1, which.min)
  ind_3_MAMIS_1 <- cbind(1:init_num, optimal_num_3_MAMIS*p-1)
  ind_3_MAMIS_2 <- cbind(1:init_num, optimal_num_3_MAMIS*p)
  ind_3_MAMIS_3 <- cbind(1:init_num, optimal_num_3_MAMIS)
  BEL_3_MAMIS <- cbind(BEL_res_3_MAMIS[ind_3_MAMIS_1], BEL_res_3_MAMIS[ind_3_MAMIS_2])
  BIC_3_MAMIS <- BIC_res_3_MAMIS[ind_3_MAMIS_3]
  MSE_3_MAMIS <- MSE_res_3_MAMIS[ind_3_MAMIS_3]
  
  return(list(BEL_3_MAMIS = BEL_3_MAMIS, BIC_3_MAMIS = BIC_3_MAMIS, MSE_3_MAMIS = MSE_3_MAMIS))
}



############# rep.tuning_5 #############
rep.tuning_5 <- function(nu, init, y, x, z) #BIC
{
  r <- ncol(z)
  n <- nrow(z)
  
  res_MAMIS <- MAMIS_2_arma(init, y, x, z, nu)
  beta_3_MAMIS <- res_MAMIS$beta_3
  
  Gee_3_MAMIS <- IVfun_arma(beta_3_MAMIS, y, x, z)
  lambda_3_MAMIS <- optimal_lam(Gee_3_MAMIS, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MAMIS <- sum(abs(lambda_3_MAMIS) > 1e-6)
  BIC_3_MAMIS <- log(sum((colMeans(Gee_3_MAMIS))^2) / r) + Rn_3_MAMIS * log(n) / n
  
  return(list(beta_3_MAMIS = beta_3_MAMIS, 
              BIC_3_MAMIS = BIC_3_MAMIS))
}


############# rep.fun_5 #############
rep.fun_5 <- function(ind, nu.seq)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  tun.num <-length(nu.seq)
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  BEL_res_3_MAMIS <- matrix(0, init_num, tun.num*p)
  BIC_res_3_MAMIS <- matrix(0, init_num, tun.num)
  MSE_res_3_MAMIS <- matrix(0, init_num, tun.num)
  
  
  for (i in c(1:init_num)) {
    cat('---- init_', i, '----')
    
    for (j in c(1:tun.num)){
      
      init <- init.data[, j]
      init <- as.matrix(init)
      nu <- nu.seq[j]
      cat("...........begin.nu_",nu,"............\n")
      
      res <- rep.tuning_5(nu, init, y, x, z)
      
      BEL_res_3_MAMIS[i, (j*p-1):(j*p)] <- res$beta_3_MAMIS
      BIC_res_3_MAMIS[i, j] <- res$BIC_3_MAMIS
      MSE_res_3_MAMIS[i,j] <- sum((res$beta_3_MAMIS - c(0.5, 0.5))^2)
    }
  }
  
  optimal_num_3_MAMIS <- apply(BIC_res_3_MAMIS,1, which.min)
  ind_3_MAMIS_1 <- cbind(1:init_num, optimal_num_3_MAMIS*p-1)
  ind_3_MAMIS_2 <- cbind(1:init_num, optimal_num_3_MAMIS*p)
  ind_3_MAMIS_3 <- cbind(1:init_num, optimal_num_3_MAMIS)
  BEL_3_MAMIS <- cbind(BEL_res_3_MAMIS[ind_3_MAMIS_1], BEL_res_3_MAMIS[ind_3_MAMIS_2])
  BIC_3_MAMIS <- BIC_res_3_MAMIS[ind_3_MAMIS_3]
  MSE_3_MAMIS <- MSE_res_3_MAMIS[ind_3_MAMIS_3]
  
  return(list(BEL_3_MAMIS = BEL_3_MAMIS, BIC_3_MAMIS = BIC_3_MAMIS, MSE_3_MAMIS = MSE_3_MAMIS))
}


############# rep.tuning_6 #############
rep.tuning_6 <- function(nu, init, y, x, z) #BIC
{
  r <- ncol(z)
  n <- nrow(z)
  
  res_MAMIS <- MAMIS_3_arma(init, y, x, z, nu)
  beta_3_MAMIS <- res_MAMIS$beta_3
  
  Gee_3_MAMIS <- IVfun_arma(beta_3_MAMIS, y, x, z)
  lambda_3_MAMIS <- optimal_lam(Gee_3_MAMIS, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MAMIS <- sum(abs(lambda_3_MAMIS) > 1e-6)
  BIC_3_MAMIS <- log(sum((colMeans(Gee_3_MAMIS))^2) / r) + Rn_3_MAMIS * log(n) / n
  
  return(list(beta_3_MAMIS = beta_3_MAMIS, 
              BIC_3_MAMIS = BIC_3_MAMIS))
}


############# rep.fun_6 #############
rep.fun_6 <- function(ind, nu.seq)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  tun.num <-length(nu.seq)
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  BEL_res_3_MAMIS <- matrix(0, init_num, tun.num*p)
  BIC_res_3_MAMIS <- matrix(0, init_num, tun.num)
  MSE_res_3_MAMIS <- matrix(0, init_num, tun.num)
  
  
  for (i in c(1:init_num)) {
    cat('---- init_', i, '----')
    
    for (j in c(1:tun.num)){
      
      init <- init.data[, j]
      init <- as.matrix(init)
      nu <- nu.seq[j]
      cat("...........begin.nu_",nu,"............\n")
      
      res <- rep.tuning_6(nu, init, y, x, z)
      
      BEL_res_3_MAMIS[i, (j*p-1):(j*p)] <- res$beta_3_MAMIS
      BIC_res_3_MAMIS[i, j] <- res$BIC_3_MAMIS
      MSE_res_3_MAMIS[i,j] <- sum((res$beta_3_MAMIS - c(0.5, 0.5))^2)
    }
  }
  
  optimal_num_3_MAMIS <- apply(BIC_res_3_MAMIS,1, which.min)
  ind_3_MAMIS_1 <- cbind(1:init_num, optimal_num_3_MAMIS*p-1)
  ind_3_MAMIS_2 <- cbind(1:init_num, optimal_num_3_MAMIS*p)
  ind_3_MAMIS_3 <- cbind(1:init_num, optimal_num_3_MAMIS)
  BEL_3_MAMIS <- cbind(BEL_res_3_MAMIS[ind_3_MAMIS_1], BEL_res_3_MAMIS[ind_3_MAMIS_2])
  BIC_3_MAMIS <- BIC_res_3_MAMIS[ind_3_MAMIS_3]
  MSE_3_MAMIS <- MSE_res_3_MAMIS[ind_3_MAMIS_3]
  
  return(list(BEL_3_MAMIS = BEL_3_MAMIS, BIC_3_MAMIS = BIC_3_MAMIS, MSE_3_MAMIS = MSE_3_MAMIS))
}