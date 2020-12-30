# experiment how go-garch is good for VaR for pure go-garch
# in comparision with true cov matrices and
# simple method for big m


options('stringsAsFactors'=FALSE)
options("max.print"=50)

library(mgarch)
library(beepr)

m <- 20
params <- gogarch.make_rand_params(m)
print(params)

OMEGA <- params$OMEGA
a <- params$a
b <- params$b

n <- 3250
calc_win <- 1250

fp <- gogarch.make_factor_process_sample(a, b, n)
R <- OMEGA %*% fp$X
los <- gogarch.make_los_sample(OMEGA %*% fp$X)

cov_data_true <- gogarch.make_sample_true_cov_data(fp,
                                                   OMEGA,
                                                   calc_win + 1)
cov_data <- gogarch.calc_cov_data(los, calc_win)

w <- rep(1, m) / m

sd_true <- gogarch.calc_sd(cov_data_true, w)
sd_garch <- gogarch.calc_sd(cov_data, w)
sd_simple_250 <- calc_sd_port_simple(los, w, 250)
sd_simple_1250 <- calc_sd_port_simple(los, w, 1250)

cls <- c(0.99, 0.95, 0.90)

for (p in cls) {
  print_var_analytics(sd_true, p, ticker = 'true')
  print_var_analytics(sd_garch, p, ticker = 'gogarch')
  print_var_analytics(sd_simple_250, p, ticker = 'simple_250')
  print_var_analytics(sd_simple_1250, p, ticker = 'simple_1250')
}

beep(5)

