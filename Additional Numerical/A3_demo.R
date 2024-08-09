library(MASS)

############# EE for IV model #############
IV.fun<-function(beta, y, x, z) {
  
  residual<- y - x %*% beta
  
  residual<- as.vector(residual)
  residual<- diag(residual)
  
  G.ee<- residual %*% z
  
  return(G.ee)
}


############ optimization of inner layer via arma  #############

optim.lambda <- function(G.ee, nu) {
  
  G.ee <- as.matrix(G.ee)
  r <- ncol(G.ee)
  n <- nrow(G.ee)
  
  lambda <- optimal_lam(G.ee, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  
  return(lambda)
}

############ objective fun of outer layer ############# 
ee.beta <- function(beta, y, x, z, nu) {
  
  G.ee <- IV.fun(beta, y, x, z)
  hat_lambda <- optim.lambda(G.ee, nu)
  
  n <- nrow(x)
  
  obj <- sum(log(1.0 + G.ee %*% hat_lambda)) - n * nu * sum(abs(hat_lambda))
  
  return(obj)
}

############# hat_beta for given data #############
hat_res <- function(y, x, z, nu)
{
  
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
  
  G.ee <- IV.fun(hat_beta, y, x, z)
  hat.lambda <- optim.lambda(G.ee, nu)
  supp.hat.lambda <- which(abs(hat.lambda) > 1e-6)
  
  p <- ncol(x)
  r <- ncol(z)
  n <- nrow(z)
  
  len.lambda <- length(supp.hat.lambda)
  
  grad.vect <- matrix(0, nrow = len.lambda, ncol = p)
  
  
  for(i in c(1:p))
  {
    grad.mat <- diag(as.vector(x[, i]*(-1))) %*% z[ , supp.hat.lambda]
    grad.vect[,i] <- colMeans(grad.mat)
  }
  
  
  temp.mat <- t(G.ee[, supp.hat.lambda]) %*% G.ee[, supp.hat.lambda]/n
  Omega <- ginv(temp.mat)
  
  hat_H <- t(grad.vect) %*% Omega %*% grad.vect
  
  return(list(hat_beta = hat_beta, hat_H = hat_H))
  
}


############# rep.tuning #############
rep.tuning <- function(nu, init, y, x, z) 
{
  r <- ncol(z)
  n <- nrow(z)
  
  hat_res <- hat_res(y, x, z, nu)
  hat_beta <- hat_res$hat_beta
  hat_Sig <- solve(n*hat_res$hat_H)
  
  res_MH <- MH_arma(init, y, x, z, nu)
  samples <- res_MH$beta_samples
  
  return(list( hat_beta = hat_beta, hat_Sig = hat_Sig, samples = samples))
}


############# rep.fun #############
rep.fun <- function(ind, nu, init)
{  
  data.IV <- Data[[ind]]
  data <- data.IV[[1]]
  y <- data[,1]
  x <- data[,2:3]
  z <- data[,-c(1:3)]
  
  p <- ncol(x)
  r <- ncol(z)    
  
  res_20 <- rep.tuning(nu, init, y[1:20], x[1:20,], z[1:20,]) 
  res_30 <- rep.tuning(nu, init, y[1:30], x[1:30,], z[1:30,]) 
  res_40 <- rep.tuning(nu, init, y[1:40], x[1:40,], z[1:40,]) 
  res_50 <- rep.tuning(nu, init, y[1:50], x[1:50,], z[1:50,]) 
  res_60 <- rep.tuning(nu, init, y[1:60], x[1:60,], z[1:60,]) 
  res_70 <- rep.tuning(nu, init, y[1:70], x[1:70,], z[1:70,]) 
  res_80 <- rep.tuning(nu, init, y[1:80], x[1:80,], z[1:80,]) 
  res_90 <- rep.tuning(nu, init, y[1:90], x[1:90,], z[1:90,]) 
  res_100 <- rep.tuning(nu, init, y[1:100], x[1:100,], z[1:100,]) 
  res_110 <- rep.tuning(nu, init, y[1:110], x[1:110,], z[1:110,]) 
  res_120 <- rep.tuning(nu, init, y[1:120], x[1:120,], z[1:120,]) 
  
  
  hat_beta_20 <- res_20$hat_beta
  hat_Sig_20 <- res_20$hat_Sig
  samples_20 <- res_20$samples
  
  hat_beta_30 <- res_30$hat_beta
  hat_Sig_30 <- res_30$hat_Sig
  samples_30 <- res_30$samples
  
  hat_beta_40 <- res_40$hat_beta
  hat_Sig_40 <- res_40$hat_Sig
  samples_40 <- res_40$samples
  
  hat_beta_50 <- res_50$hat_beta
  hat_Sig_50 <- res_50$hat_Sig
  samples_50 <- res_50$samples
  
  hat_beta_60 <- res_60$hat_beta
  hat_Sig_60 <- res_60$hat_Sig
  samples_60 <- res_60$samples
  
  hat_beta_70 <- res_70$hat_beta
  hat_Sig_70 <- res_70$hat_Sig
  samples_70 <- res_70$samples
  
  hat_beta_80 <- res_80$hat_beta
  hat_Sig_80 <- res_80$hat_Sig
  samples_80 <- res_80$samples
  
  hat_beta_90 <- res_90$hat_beta
  hat_Sig_90 <- res_90$hat_Sig
  samples_90 <- res_90$samples
  
  hat_beta_100 <- res_100$hat_beta
  hat_Sig_100 <- res_100$hat_Sig
  samples_100 <- res_100$samples
  
  hat_beta_110 <- res_110$hat_beta
  hat_Sig_110 <- res_110$hat_Sig
  samples_110 <- res_110$samples
  
  hat_beta_120 <- res_120$hat_beta
  hat_Sig_120 <- res_120$hat_Sig
  samples_120 <- res_120$samples
  
  res_all <- rbind(cbind(hat_beta_20, hat_Sig_20, samples_20), cbind(hat_beta_30, hat_Sig_30, samples_30), 
                   cbind(hat_beta_40, hat_Sig_40, samples_40), cbind(hat_beta_50, hat_Sig_50, samples_50),
                   cbind(hat_beta_60, hat_Sig_60, samples_60), cbind(hat_beta_70, hat_Sig_70, samples_70), 
                   cbind(hat_beta_80, hat_Sig_80, samples_80), cbind(hat_beta_90, hat_Sig_90, samples_90),
                   cbind(hat_beta_100, hat_Sig_100, samples_100), cbind(hat_beta_110, hat_Sig_110, samples_110), 
                   cbind(hat_beta_120, hat_Sig_120, samples_120))
  
  return(list(res_all = t(res_all)))
}
