options('stringsAsFactors'=FALSE)
options("max.print"=50)

library(mgarch)
library(beepr)

OMEGA <- rbind(c(-2, 1, 1),
               c(1, -1, 2),
               c(3, 0, 1))

a <- c(0.1, 0.2, 0.05)
b <- c(0.8, 0.6, 0.9)

n <- 4250
calc_win <- 1250

fp <- gogarch.make_factor_process_sample(a, b, n)
los <- gogarch.make_los_sample(OMEGA %*% fp$X)

cov_data_true <- gogarch.make_sample_true_cov_data(fp, OMEGA)
cov_data <- gogarch.calc_cov_data(los, calc_win)

w <- c(0.3, 0.3, 0.4)

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

beep(3)
