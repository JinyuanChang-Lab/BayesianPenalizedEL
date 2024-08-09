library(MASS)

############# hat_beta for given data #############
hat_beta <- function(ind, nu)
{
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  cat("...........begin.search............\n") 
  
  beta.1 <- seq(-0.5,1.5,0.02)
  beta.2 <- seq(-0.5,1.5,0.02)
  
  n_box <- length(seq(-0.5,1.5,0.02))
  
  ee_beta_value <- matrix(0, n_box, n_box)
  for (i in c(1:n_box)) {
    for (j in c(1:n_box)) {
      beta <- c(beta.1[i], beta.2[j])
      ee_beta_value[i,j] <- ee.beta(beta, y, x, z, nu)
    }
  }
  
  hat_beta.1 <- beta.1[which(ee_beta_value == min(ee_beta_value), arr.ind = TRUE)[1]]
  hat_beta.2 <- beta.2[which(ee_beta_value == min(ee_beta_value), arr.ind = TRUE)[2]]
  hat_beta <- c(hat_beta.1, hat_beta.2)
  
  cat("...........end.search............\n")
  
  return(hat_beta)
  
}



############# rep.fun #############
rep.fun <- function(ind, nu)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  n <- nrow(z)     #n should be a global variable 
  r <- ncol(z)     #r should be a global variable
  
  init_num <- 49
  init.data <- matrix(0, 2, init_num)
  
  for (i in c(1:7)) {
    init.data[1, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)[i]
    init.data[2, (7*i-6):(7*i)] = seq(-3,4, length.out = 7)
  }
  
  beta_MH_all_1 <- matrix(0, p, init_num)
  beta_MH_all_2 <- matrix(0, p, init_num)
  beta_MH_all_3 <- matrix(0, p, init_num)
  accept_rio_all <- rep(0,init_num)
  
  beta_MAMIS_all_1 <- matrix(0, p, init_num)
  beta_MAMIS_all_2 <- matrix(0, p, init_num)
  beta_MAMIS_all_3 <- matrix(0, p, init_num)
  
  beta_optim_all <- matrix(0, p, init_num)
  beta_nlm_all <- matrix(0, p, init_num)
  
  
  for (j in c(1:init_num)) {
    
    init <- init.data[, j]
    init <- as.matrix(init)
    # print(init)
    
    cat('---- init_', j, '----')
    
    cat("...........begin.MH............\n")    ##### MH
    
    res_MH <- MH_arma(init, y, x, z, nu)
    
    beta_MH_1 <- res_MH$beta_mean_1
    beta_MH_2 <- res_MH$beta_mean_2
    beta_MH_3 <- res_MH$beta_mean_3
    accept_rio <- res_MH$accept_rio
    
    beta_MH_all_1[, j] <- beta_MH_1
    beta_MH_all_2[, j] <- beta_MH_2
    beta_MH_all_3[, j] <- beta_MH_3
    accept_rio_all[j] <- accept_rio
    
    cat("...........end.MH...........\n")
    
    
    cat("...........begin.MAMIS............\n")    ##### MAMIS
    
    res_MAMIS <- MAMIS_arma(init, y, x, z, nu)
    
    beta_MAMIS_1 <- res_MAMIS$beta_1
    beta_MAMIS_2 <- res_MAMIS$beta_2
    beta_MAMIS_3 <- res_MAMIS$beta_3
    beta_MAMIS_all_1[, j] <- beta_MAMIS_1
    beta_MAMIS_all_2[, j] <- beta_MAMIS_2
    beta_MAMIS_all_3[, j] <- beta_MAMIS_3
    
    cat("...........end.MAMIS...........\n")
    
    
    cat("...........begin.optim............\n")    ##### optim
    
    res_optim <- optim(init, ee.beta, y = y, x = x, z = z, nu = nu)
    
    beta_optim <- res_optim$par
    beta_optim_all[, j] <- beta_optim
    
    cat("...........end.optim............\n")
    
    
    cat("...........begin.nlm............\n")     ##### nlm
    
    res_nlm <- nlm(ee.beta, init, y = y, x = x, z = z, nu =nu)
    
    beta_nlm <- res_nlm$estimate
    beta_nlm_all[, j] <- beta_nlm
    
    cat("...........end.nlm............\n")
  
  }
  
  return(list(beta_MH_all_1 = beta_MH_all_1, beta_MH_all_2 = beta_MH_all_2, beta_MH_all_3 = beta_MH_all_3, accept_rio_all = accept_rio_all,
              beta_MAMIS_all_1 = beta_MAMIS_all_1, beta_MAMIS_all_2 = beta_MAMIS_all_2, beta_MAMIS_all_3 = beta_MAMIS_all_3,
              beta_optim_all = beta_optim_all, beta_nlm_all = beta_nlm_all))
}

