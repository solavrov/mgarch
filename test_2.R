library(mgarch)

OMEGA <- rbind(c(-2, 1, 1),
               c(1, -1, 2),
               c(3, 0, 1))

a <- c(0.1, 0.2, 0.05)
b <- c(0.8, 0.6, 0.9)

n <- 10^6

SIGMA <- OMEGA %*% t(OMEGA)

factor_process <-  gogarch.make_factor_process_sample(a, b, n)
R <- OMEGA %*% factor_process$X

SIGMA_est <- gogarch.find_SIGMA_estim(R)

OMEGA_est <- gogarch.find_OMEGA_estim(R, SIGMA_est, df=5)

print(OMEGA)
print(OMEGA_est)
print(gogarch.permutate_colunms(OMEGA_est, OMEGA))

X_est <- gogarch.find_X_process_estim(R, OMEGA_est)
gp <- gogarch.find_garch_param_estim(X_est)

print(gp)

SD2_est <- gogarch.find_SD2_estim(X_est, gp)


stock1 <- 2
stock2 <- 3
corr_12 <- c()
corr_12_est <- c()
var_1 <- c()
var_1_est <- c()
var_2 <- c()
var_2_est <- c()

plot_range <- 6000:6500

for (t in plot_range) {
  SIGMA_t <- gogarch.find_conditional_SIGMA(factor_process$SD2, OMEGA, t)
  SIGMA_t_est <- gogarch.find_conditional_SIGMA(SD2_est, OMEGA_est, t)
  corr_12[t - plot_range[1] + 1] <- cov2cor(SIGMA_t)[stock1,stock2]
  corr_12_est[t - plot_range[1] + 1] <- cov2cor(SIGMA_t_est)[stock1,stock2]
  var_1[t - plot_range[1] + 1] <- SIGMA_t[stock1,stock1]
  var_1_est[t - plot_range[1] + 1] <- SIGMA_t_est[stock1,stock1]
  var_2[t - plot_range[1] + 1] <- SIGMA_t[stock2,stock2]
  var_2_est[t - plot_range[1] + 1] <- SIGMA_t_est[stock2,stock2]
}

plot(corr_12, type='l', col='gray', ylim=c(0.4,0.8), lwd=4)
lines(corr_12_est, type='l', col='red', lwd=1)

plot(var_1, type='l', col='gray', lwd=4)
lines(var_1_est, type='l', col='red', lwd=1)

plot(var_2, type='l', col='gray', lwd=4)
lines(var_2_est, type='l', col='red', lwd=1)

