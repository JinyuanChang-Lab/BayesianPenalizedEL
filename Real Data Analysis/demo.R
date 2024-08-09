library(MASS)

############# rep.tuning #############
rep.tuning <- function(init, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF, nu) #BIC
{
  r <- ncol(y)
  n <- nrow(y)
  
  res_MAMIS <- MAMIS_arma(init, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF, nu)
  res_MH <- MH_arma(init, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF, nu)
  
  beta_3_MAMIS <- res_MAMIS$beta_3
  beta_3_MH <- res_MH$beta_mean_3
  
  Gee_3_MAMIS <- GEE_arma(beta_3_MAMIS, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF)
  lambda_3_MAMIS <- optimal_lam(Gee_3_MAMIS, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MAMIS <- sum(abs(lambda_3_MAMIS) > 1e-6)
  BIC_3_MAMIS <- log(sum((colMeans(Gee_3_MAMIS))^2) / r) + Rn_3_MAMIS * log(n) / n
  
  
  Gee_3_MH <- GEE_arma(beta_3_MH, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF)
  lambda_3_MH <- optimal_lam(Gee_3_MH, nu, 1e-9, 1e-3, 5000, FALSE, "identity")
  Rn_3_MH <- sum(abs(lambda_3_MH) > 1e-6)
  BIC_3_MH <- log(sum((colMeans(Gee_3_MH))^2) / r) + Rn_3_MH * log(n) / n
  
  return(list(beta_3_MAMIS = beta_3_MAMIS, BIC_3_MAMIS = BIC_3_MAMIS, 
    beta_3_MH = beta_3_MH, BIC_3_MH = BIC_3_MH))
}


############# rep.fun #############
rep.fun <- function(ind, nu, init_all, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF)
{  
  
  r <- ncol(y)
  n <- nrow(y)
  p <- 5

  init <- init_all[, ind]
  
  res <- rep.tuning(init, y, Ran.v, Ran.a, Ran.h, X_bar, N_nF, nu)
  
  BEL_3_MAMIS <- res$beta_3_MAMIS
  BIC_3_MAMIS <- res$BIC_3_MAMIS
  
  BEL_3_MH <- res$beta_3_MH
  BIC_3_MH <- res$BIC_3_MH
  
  return(list(BEL_3_MAMIS = BEL_3_MAMIS, BIC_3_MAMIS = BIC_3_MAMIS, 
    BEL_3_MH = BEL_3_MH, BIC_3_MH = BIC_3_MH))
}