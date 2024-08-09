library(MASS)
library(CVXR)


############# EE for IV model #############
IV.fun<-function(beta, y, x, z) {
  
  residual<- y - sin(x %*% beta)
  
  residual<- as.vector(residual)
  residual<- diag(residual)
  
  G.ee<- residual %*% z
  
  return(G.ee)
}


############ optimization of inner layer via CVXR #############

optim.lambda <- function(G.ee, nu) {
  
  G.ee <- as.matrix(G.ee)
  r <- ncol(G.ee)
  n <- nrow(G.ee)
  
  lambda <- Variable(r)
  
  obj <- sum(log(1.0 + G.ee %*% lambda)) - n * nu * p_norm(lambda, 1)
  
  prob <- Problem(Maximize(obj))
  
  res1 <- solve(prob, solver = "ECOS_BB")
  
  if(res1$status=="solver_error")
  {
    res <- solve(prob, solver = "SCS")
  }else {
    res <- res1
  }
  
  lambda <- res$getValue(lambda)
  
  return(lambda)
}


############ objective fun of outer layer ############# 
ee.beta <- function(beta, y, x, z, nu) {
  
  if(max(abs(beta))<10)
  {
    G.ee <- IV.fun(beta, y, x, z)
    hat_lambda <- optim.lambda(G.ee, nu)
    
    n <- nrow(x)
    
    obj <- sum(log(1.0 + G.ee %*% hat_lambda)) - n * nu * sum(abs(hat_lambda))
  }else{
    obj <- 0
  }

  return(obj)
}


