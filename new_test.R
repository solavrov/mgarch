# experiment how good go-garch in
# finding params estimates of pure go-garch process

options('stringsAsFactors'=FALSE)
options("max.print"=50)

library(mgarch)

OMEGA <- rbind(c(-2, 1, 1),
               c(1, -1, 2),
               c(3, 0, 1))

# a <- c(0.1, 0.2, 0.05)
# b <- c(0.8, 0.6, 0.9)

a <- c(0.25, 0.49, 0.38)
b <- c(0.7, 0.5, 0.6)

n <- 1250

SIGMA <- OMEGA %*% t(OMEGA)

factor_process <-  gogarch.make_factor_process_sample(a, b, n)
R <- OMEGA %*% factor_process$X

SIGMA_est <- gogarch.find_SIGMA_estim(R)

OMEGA_est <- gogarch.find_OMEGA_estim(R, SIGMA_est, df=5)

OMEGA_est <- gogarch.permutate_colunms(OMEGA_est, OMEGA)

print(OMEGA)
print(OMEGA_est)

X_est <- gogarch.find_X_process_estim(R, OMEGA_est)
gp <- gogarch.find_garch_param_estim(X_est)

print(gp)

SD2_est <- gogarch.find_SD2_estim(X_est, gp)


stock1 <- 1
stock2 <- 2
stock3 <- 3

corr_12 <- c()
corr_12_est <- c()
corr_23 <- c()
corr_23_est <- c()
corr_13 <- c()
corr_13_est <- c()

var_1 <- c()
var_1_est <- c()
var_2 <- c()
var_2_est <- c()
var_3 <- c()
var_3_est <- c()

plot_range <- 500:1000

for (t in plot_range) {
  SIGMA_t <- gogarch.find_conditional_SIGMA(factor_process$SD2, OMEGA, t)
  SIGMA_t_est <- gogarch.find_conditional_SIGMA(SD2_est, OMEGA_est, t)

  corr_12[t - plot_range[1] + 1] <- cov2cor(SIGMA_t)[stock1,stock2]
  corr_12_est[t - plot_range[1] + 1] <- cov2cor(SIGMA_t_est)[stock1,stock2]

  corr_23[t - plot_range[1] + 1] <- cov2cor(SIGMA_t)[stock2,stock3]
  corr_23_est[t - plot_range[1] + 1] <- cov2cor(SIGMA_t_est)[stock2,stock3]

  corr_13[t - plot_range[1] + 1] <- cov2cor(SIGMA_t)[stock1,stock3]
  corr_13_est[t - plot_range[1] + 1] <- cov2cor(SIGMA_t_est)[stock1,stock3]

  var_1[t - plot_range[1] + 1] <- SIGMA_t[stock1,stock1]
  var_1_est[t - plot_range[1] + 1] <- SIGMA_t_est[stock1,stock1]

  var_2[t - plot_range[1] + 1] <- SIGMA_t[stock2,stock2]
  var_2_est[t - plot_range[1] + 1] <- SIGMA_t_est[stock2,stock2]

  var_3[t - plot_range[1] + 1] <- SIGMA_t[stock3,stock3]
  var_3_est[t - plot_range[1] + 1] <- SIGMA_t_est[stock3,stock3]
}

plot(corr_12, type='l', col='gray', ylim=c(min(corr_12), max(corr_12)), lwd=4)
lines(corr_12_est, type='l', col='red', lwd=1)

plot(corr_23, type='l', col='gray', ylim=c(min(corr_23), max(corr_23)), lwd=4)
lines(corr_23_est, type='l', col='red', lwd=1)

plot(corr_13, type='l', col='gray', ylim=c(min(corr_13), max(corr_13)), lwd=4)
lines(corr_13_est, type='l', col='red', lwd=1)

plot(var_1, type='l', col='gray', lwd=4)
lines(var_1_est, type='l', col='red', lwd=1)

plot(var_2, type='l', col='gray', lwd=4)
lines(var_2_est, type='l', col='red', lwd=1)

plot(var_3, type='l', col='gray', lwd=4)
lines(var_3_est, type='l', col='red', lwd=1)

