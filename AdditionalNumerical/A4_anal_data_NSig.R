library(transport)
library(MASS)
library(ggplot2)
######################### r = 50 #########################
res <- read.csv(paste0('result/NSigma_all_nu_0.03_r_50.csv')
                , header = F)

data <- res
data <- data[-1,]
data <- data[, -1]
data <- as.matrix(data)


dis_1 <- rep(0,500)
dis_2 <- rep(0,500)
for (i in seq(31,16500,33)) {
  cat('---- Simulation Time = ', i, '----')
  Pos_samples <- as.numeric(data[, i])
  PEL_samples <- as.numeric(data[, (i+1)])
  EL_samples <- as.numeric(data[, (i+2)])
  
  dis_1[ceiling(i/33)] <- wasserstein1d(Pos_samples, PEL_samples)
  dis_2[ceiling(i/33)] <- wasserstein1d(Pos_samples, EL_samples)
}

dis <- cbind(dis_1, dis_2)

write.csv(dis, paste0('result/NSigma_dis_120_r_50.csv'))



######################### r = 70 #########################
res <- read.csv(paste0('result/NSigma_all_nu_0.03_r_70.csv')
                , header = F)

data <- res
data <- data[-1,]
data <- data[, -1]
data <- as.matrix(data)


dis_1 <- rep(0,500)
dis_2 <- rep(0,500)
for (i in seq(31,16500,33)) {
  cat('---- Simulation Time = ', i, '----')
  Pos_samples <- as.numeric(data[, i])
  PEL_samples <- as.numeric(data[, (i+1)])
  EL_samples <- as.numeric(data[, (i+2)])
  
  dis_1[ceiling(i/33)] <- wasserstein1d(Pos_samples, PEL_samples)
  dis_2[ceiling(i/33)] <- wasserstein1d(Pos_samples, EL_samples)
}

dis <- cbind(dis_1, dis_2)

write.csv(dis, paste0('result/NSigma_dis_120_r_70.csv'))


###################### plot r = 50 ############################
r = 50 

Data <- NULL

for (j in seq(20,120,10)){
  res <- read.csv(paste0('result/NSigma_dis_', j,'_r_', r, '.csv'), header = F)
  res <- res[-1,]
  res <- res[, -1]
  res <- as.matrix(res)
  
  Data <- cbind(Data, res)
}
data <- matrix(as.numeric(Data), nrow=500)


dim <- c(seq(20,120, 10), seq(20,120, 10))
dis <- c(colMeans(data)[seq(1,22,2)], colMeans(data)[seq(2,22,2)])
label <- c(rep("with penalty",11), rep("without penalty", 11))
Res <- data.frame(dim, dis, label)
p <- ggplot(Res, aes(x = dim, y = dis, colour = label))

p + geom_line(aes(color = label))+geom_point(aes(shape = label, color = label), size=2)+ ylim(0.030, 1.02) + scale_shape_manual(values = c(15,16)) + theme_bw() + 
  theme(legend.position = c(0.3,0.75), panel.grid = element_blank(), legend.text = element_text(size = 20, family = "serif"), legend.title = element_blank()) + 
  scale_x_continuous(breaks= seq(20, 120, 20)) + labs(x = 'sample size', y = 'Wasserstein distance') + 
  theme(axis.title.x=element_text(vjust=1, size=18, family = "serif"), axis.text.x=element_text(size=16, family = "serif"),
        axis.title.y=element_text(vjust=0, size=18, family = "serif"), axis.text.y=element_text(size=16, family = "serif"))



###################### plot r = 70 ############################
r = 70 

Data <- NULL

for (j in seq(20,120,10)){
  res <- read.csv(paste0('result/NSigma_dis_', j,'_r_', r, '.csv'), header = F)
  res <- res[-1,]
  res <- res[, -1]
  res <- as.matrix(res)
  
  Data <- cbind(Data, res)
}
data <- matrix(as.numeric(Data), nrow=500)


dim <- c(seq(20,120, 10), seq(20,120, 10))
dis <- c(colMeans(data)[seq(1,22,2)], colMeans(data)[seq(2,22,2)])
label <- c(rep("with penalty",11), rep("without penalty", 11))
Res <- data.frame(dim, dis, label)
p <- ggplot(Res, aes(x = dim, y = dis, colour = label))

p + geom_line(aes(color = label))+geom_point(aes(shape = label, color = label), size=2)+ ylim(0.030, 1.02) + scale_shape_manual(values = c(15,16)) + theme_bw() + 
  theme(legend.position = c(0.3,0.75), panel.grid = element_blank(), legend.text = element_text(size = 20, family = "serif"), legend.title = element_blank()) + 
  scale_x_continuous(breaks= seq(20, 120, 20)) + labs(x = 'sample size', y = 'Wasserstein distance') + 
  theme(axis.title.x=element_text(vjust=1, size=18, family = "serif"), axis.text.x=element_text(size=16, family = "serif"),
        axis.title.y=element_text(vjust=0, size=18, family = "serif"), axis.text.y=element_text(size=16, family = "serif"))





