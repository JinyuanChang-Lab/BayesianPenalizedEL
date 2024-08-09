Rcpp::sourceCpp("MH.cpp")

library(parallel)
library(MASS)
library(R.matlab)
Data_TT <- readMat("data_SES.mat")
Data <- Data_TT$data.test.time


rep.num <- 500   # repeat times for simulation

start <- Sys.time()

res.mat <- mclapply(1:rep.num, function(rep.sim){
  
  cat('---- Simulation Time = ', rep.sim, '----- \rep.sim')
  
  iter_num_pen <- rep(0,20)
  accept_pen <- rep(0,20)
  iter_num_nopen <- rep(0,20)
  accept_nopen <- rep(0,20)
  
  for (i in c(1:20)) {
    cat('---- r dimension time = ', i, '-----')
    data <- Data[[i]]
    data <- data[[1]]
    y <- data[,1]
    x <- data[,2:3]
    z <- data[,-c(1:3)]
    
    res_pen <- MH_arma(c(0.3, 0.3), y, x, z, 0.03)
    res_nopen <- MH_arma(c(0.3, 0.3), y, x, z, 0)
    
    iter_num_pen[i] <- res_pen$iter_num
    accept_pen[i] <- res_pen$accept
    
    iter_num_nopen[i] <- res_nopen$iter_num
    accept_nopen[i] <- res_nopen$accept
  }
  
  res_time <- rbind(iter_num_pen, accept_pen, iter_num_nopen, accept_nopen)

  write.csv(res_time, paste0('result/res_time_', rep.sim, '.csv'))
  
  return(res_time)
}
, mc.cores = 100)

end <- Sys.time()
print(end-start)

write.csv(res.mat, paste0('result/rcpp_res_time_', rep.num, '.csv'))



################### data processing ######################

res <- read.csv(paste0('result/rcpp_res_time_500.csv')
                    , header = F)

res <- res[-1,]
res <- res[,-1]

iter_pen_res <- as.numeric(res[1,])
iter_pen_mat <- matrix(iter_pen_res, nrow=500, byrow = T)
iter_pen <- colMeans(iter_pen_mat)

iter_nopen_res <- as.numeric(res[3,])
iter_nopen_mat <- matrix(iter_nopen_res, nrow=500, byrow = T)
iter_nopen <- colMeans(iter_nopen_mat)


##################### plot ###############################

library(ggplot2)

dim <- c(seq(50,1000,by=50),seq(50,1000,by=50))
Iterations <- c(iter_nopen, iter_pen)
label <- c(rep("without penalty",20), rep("with penalty",20))
Res <- data.frame(dim, Iterations, label)
p <- ggplot(Res, aes(x = dim, y = Iterations, colour = label))
p + geom_line(aes(color = label))+geom_point(aes(shape = label, color = label), size=2)+ scale_shape_manual(values = c(15,16)) + theme_bw() + 
  theme(legend.position = c(0.2,0.85), panel.grid = element_blank(), legend.text = element_text(size = 15, family = "serif"), legend.title = element_blank()) + 
  scale_x_continuous(breaks=seq(0,1000, by=200)) + labs(x=expression(italic(r))) + 
  theme(axis.title.x=element_text(vjust=1, size=18, family = "serif"), axis.text.x=element_text(size=16, family = "serif"),
        axis.title.y=element_text(vjust=0, size=18, family = "serif"), axis.text.y=element_text(size=16, family = "serif"))

