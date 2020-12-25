options('stringsAsFactors'=FALSE)
options("max.print"=50)

library(mgarch)

OMEGA <- rbind(c(-2, 1, 1),
               c(1, -1, 2),
               c(3, 0, 1))

a <- c(0.2, 0.3, 0.05)
b <- c(0.3, 0.6, 0.9)

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

t <- n/2

SIGMA_t <- gogarch.find_conditional_SIGMA(factor_process$SD2, OMEGA, t)
SIGMA_t_est <- gogarch.find_conditional_SIGMA(SD2_est, OMEGA_est, t)

print(SIGMA_t)
print(SIGMA_t_est)
print(cov2cor(SIGMA_t))
print(cov2cor(SIGMA_t_est))

