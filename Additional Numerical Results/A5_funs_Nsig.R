library(abc)
library(BSL)


simulate_data <- function(theta, n) {
  y <- rnorm(n, mean = 0, sd = abs(theta))
  return(y)
}

compute_summary <- function(data) {
  sd_stat <- sd(data)
  return(sd_stat)
}


sim_bsl <- function(x, n) {
  M <- 7000  
  cov_rw <- matrix(0.1, 1, 1)
  y_obs <- x[1:n]
  
  model <- newModel(
    fnSim = function(theta) simulate_data(theta, n=n), 
    fnSum = compute_summary,                    
    theta0 = c(0.7)                             
  )
  
  result_bsl <- bsl(
    y = x[1:n],
    n = 10,
    M = M,          
    model = model,                                
    covRandWalk = cov_rw,                         
    method = "BSL"                                
  )
  
  return(result_bsl@theta[2001:7000])
}


sim_abc <- function(x, n, M, tol) {
  
  simulated_thetas <- runif(M, 0, 2)
  param_names <- c("theta")
  
  summary_names <- c("sd")
  param_matrix <- matrix(simulated_thetas, ncol = 1, dimnames = list(NULL, param_names))
  
  simulated_summaries <- matrix(NA, nrow = length(simulated_thetas), ncol = 1)
  for (i in 1:length(simulated_thetas)) {
    simulated_data <- simulate_data(simulated_thetas[i], n=n)
    simulated_summaries[i, ] <- compute_summary(data=simulated_data)
  }
  
  sumstat_matrix <- simulated_summaries
  colnames(sumstat_matrix) <- summary_names
  
  abc_result <- abc(
    target = compute_summary(data = x[1:n]),       
    param = param_matrix,            
    sumstat = sumstat_matrix,        
    tol = tol,                 
    method = "rejection"             
  )
  
  return(abc_result$unadj.values)
}