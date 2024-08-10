library(transport)
library(MASS)

######################### r = 80 #########################
res <- read.csv(paste0('result/ks_test_all_nu_0.03_r_80.csv')
                , header = F)

data <- res
data <- data[-1,]
data <- data[, -1]
data <- as.matrix(data)

dis_20 <- rep(0,500)
for (i in seq(1,11000,22)) {
  hat_beta_20 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_20 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_20, hat_Sig_20)
  samples_20 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_20[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_20), p = 2)
}
write.csv(dis_20, paste0('result/dis_20.csv'))


dis_30 <- rep(0,500)
for (i in seq(3,11000,22)) {
  hat_beta_30 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_30 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_30, hat_Sig_30)
  samples_30 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_30[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_30), p = 2)
}
write.csv(dis_30, paste0('result/dis_30.csv'))


dis_40 <- rep(0,500)
for (i in seq(5,11000,22)) {
  hat_beta_40 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_40 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_40, hat_Sig_40)
  samples_40 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_40[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_40), p = 2)
}
write.csv(dis_40, paste0('result/dis_40.csv'))


dis_50 <- rep(0,500)
for (i in seq(7,11000,22)) {
  hat_beta_50 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_50 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_50, hat_Sig_50)
  samples_50 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_50[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_50), p = 2)
}
write.csv(dis_50, paste0('result/dis_50.csv'))


dis_60 <- rep(0,500)
for (i in seq(9,11000,22)) {
  hat_beta_60 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_60 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_60, hat_Sig_60)
  samples_60 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_60[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_60), p = 2)
}
write.csv(dis_60, paste0('result/dis_60.csv'))


dis_70 <- rep(0,500)
for (i in seq(11,11000,22)) {
  hat_beta_70 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_70 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_70, hat_Sig_70)
  samples_70 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_70[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_70), p = 2)
}
write.csv(dis_70, paste0('result/dis_70.csv'))


dis_80 <- rep(0,500)
for (i in seq(13,11000,22)) {
  hat_beta_80 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_80 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_80, hat_Sig_80)
  samples_80 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_80[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_80), p = 2)
}
write.csv(dis_80, paste0('result/dis_80.csv'))


dis_90 <- rep(0,500)
for (i in seq(15,11000,22)) {
  hat_beta_90 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_90 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_90, hat_Sig_90)
  samples_90 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_90[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_90), p = 2)
}
write.csv(dis_90, paste0('result/dis_90.csv'))


dis_100 <- rep(0,500)
for (i in seq(17,11000,22)) {
  hat_beta_100 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_100 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_100, hat_Sig_100)
  samples_100 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_100[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_100), p = 2)
}
write.csv(dis_100, paste0('result/dis_100.csv'))


dis_110 <- rep(0,500)
for (i in seq(19,11000,22)) {
  hat_beta_110 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_110 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_110, hat_Sig_110)
  samples_110 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_110[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_110), p = 2)
}
write.csv(dis_110, paste0('result/dis_110.csv'))


dis_120 <- rep(0,500)
for (i in seq(21,11000,22)) {
  hat_beta_120 <- as.numeric(data[1, i:(i+1)])
  hat_Sig_120 <- matrix(as.numeric(data[2:3, i:(i+1)]), 2 ,2)
  N_samples <- mvrnorm(n=5000, hat_beta_120, hat_Sig_120)
  samples_120 <- matrix(as.numeric(data[4:nrow(data), i:(i+1)]), (nrow(data)-3), 2)
  dis_120[ceiling(i/22)] <- wasserstein(pp(N_samples), pp(samples_120), p = 2)
}

write.csv(dis_120, paste0('result/dis_120.csv'))


##################################################
dis_20 <- read.csv(paste0('result/r=80/dis_20.csv')
                , header = F)
dis_20 <- dis_20 [-1,]
dis_20 <- dis_20 [, -1]
dis_20 <- as.numeric(dis_20)


dis_30 <- read.csv(paste0('result/r=80/dis_30.csv')
                   , header = F)
dis_30 <- dis_30 [-1,]
dis_30 <- dis_30 [, -1]
dis_30 <- as.numeric(dis_30)


dis_40 <- read.csv(paste0('result/r=80/dis_40.csv')
                   , header = F)
dis_40 <- dis_40 [-1,]
dis_40 <- dis_40 [, -1]
dis_40 <- as.numeric(dis_40)


dis_50 <- read.csv(paste0('result/r=80/dis_50.csv')
                   , header = F)
dis_50 <- dis_50 [-1,]
dis_50 <- dis_50 [, -1]
dis_50 <- as.numeric(dis_50)


dis_60 <- read.csv(paste0('result/r=80/dis_60.csv')
                   , header = F)
dis_60 <- dis_60 [-1,]
dis_60 <- dis_60 [, -1]
dis_60 <- as.numeric(dis_60)


dis_70 <- read.csv(paste0('result/r=80/dis_70.csv')
                   , header = F)
dis_70 <- dis_70 [-1,]
dis_70 <- dis_70 [, -1]
dis_70 <- as.numeric(dis_70)


dis_80 <- read.csv(paste0('result/r=80/dis_80.csv')
                   , header = F)
dis_80 <- dis_80 [-1,]
dis_80 <- dis_80 [, -1]
dis_80 <- as.numeric(dis_80)


dis_90 <- read.csv(paste0('result/r=80/dis_90.csv')
                   , header = F)
dis_90 <- dis_90 [-1,]
dis_90 <- dis_90 [, -1]
dis_90 <- as.numeric(dis_90)


dis_100 <- read.csv(paste0('result/r=80/dis_100.csv')
                   , header = F)
dis_100 <- dis_100 [-1,]
dis_100 <- dis_100 [, -1]
dis_100 <- as.numeric(dis_100)


dis_110 <- read.csv(paste0('result/r=80/dis_110.csv')
                   , header = F)
dis_110 <- dis_110 [-1,]
dis_110 <- dis_110 [, -1]
dis_110 <- as.numeric(dis_110)


dis_120 <- read.csv(paste0('result/r=80/dis_120.csv')
                   , header = F)
dis_120 <- dis_120 [-1,]
dis_120 <- dis_120 [, -1]
dis_120 <- as.numeric(dis_120)


x=seq(20,120,by=10)
y=c(mean(dis_20), mean(dis_30), mean(dis_40), mean(dis_50), mean(dis_60), mean(dis_70), 
    mean(dis_80), mean(dis_90), mean(dis_100), mean(dis_110), mean(dis_120))

plot(x, y, xlab='r=80', ylab='Wasserstein distance')
