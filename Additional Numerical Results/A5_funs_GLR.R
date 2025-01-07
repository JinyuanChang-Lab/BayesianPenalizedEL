library(abc)
library(BSL)


simulate_data <- function(theta, x, n) {
  x <- x[1:n]
  p <- 1/(1+exp(-theta * x))
  y <- rbinom(n, 1, p)
  return(y)
}

compute_summary <- function(data, x) {
  model <- glm(data ~ x - 1, family=binomial)
  coef <- coefficients(model)
  return(coef[1])
}


sim_bsl <- function(y, z, n) {
  M <- 7000  
  cov_rw <- matrix(0.1, 1, 1)
  y_obs <- y[1:n]
  
  model <- newModel(
    fnSim = function(theta) simulate_data(theta, x=z[1:n], n=n), 
    fnSum = function(data) compute_summary(data, x=z[1:n]),                    
    theta0 = c(0.1)                             
  )
  
  result_bsl <- bsl(
    y = y[1:n],
    n = 10,
    M = M,          
    model = model,                                
    covRandWalk = cov_rw,                         
    method = "BSL"                                
  )
  
  return(result_bsl@theta[2001:7000])
}


sim_abc <- function(y, z, n, M, tol) {
  
  simulated_thetas <- runif(M, -0.5, 0.9)
  param_names <- c("theta")
  
  summary_names <- c("slope")
  param_matrix <- matrix(simulated_thetas, ncol = 1, dimnames = list(NULL, param_names))
  
  simulated_summaries <- matrix(NA, nrow = length(simulated_thetas), ncol = 1)
  for (i in 1:length(simulated_thetas)) {
    simulated_data <- simulate_data(simulated_thetas[i], x=z[1:n], n=n)
    simulated_summaries[i, ] <- compute_summary(data=simulated_data, x=z[1:n])
  }
  
  sumstat_matrix <- simulated_summaries
  colnames(sumstat_matrix) <- summary_names
  
  abc_result <- abc(
    target = compute_summary(data = y[1:n], x=z[1:n]),       
    param = param_matrix,            
    sumstat = sumstat_matrix,        
    tol = tol,                 
    method = "rejection"             
  )
  
  return(abc_result$unadj.values)
}