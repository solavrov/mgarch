options("max.print"=50)
library(mgarch)


t <- 'sber'
data_stock <- load_stock_data_by_ticker('sber')
tail_len <- 2000
calc_win <- c(250, 500, 1000, 1500, 2000, 3000)
conf_level_var <- c(0.99, 0.95, 0.90)
start_win <- 3000

cat(t, 'data len =', nrow(data_stock), '\n')

PV_simple <- c()
for (cl in conf_level_var) {
  pv_w <- c()
  for (i in 1:length(calc_win)) {
    w <- calc_win[i]
    cat('\n\nsimple windowed', cl, w)
    sd <- calc_sd_simple(data_stock, w)
    pv_w[i] <- analytics.p_value(sd, cl, tail_len = tail_len)
  }
  PV_simple <- rbind(PV_simple, pv_w)
}

PV_garch <- c()
for (cl in conf_level_var) {
  pv_w <- c()
  for (i in 1:length(calc_win)) {
    w <- calc_win[i]
    cat('\n\ngarch windowed', cl, w)
    sd <- calc_sd_garch(data_stock, w)
    pv_w[i] <- analytics.p_value(sd, cl, tail_len = tail_len)
  }
  PV_garch <- rbind(PV_garch, pv_w)
}

PV_garch_not <- c()
for (i in 1:length(conf_level_var)) {
  cl <- conf_level_var[i]
  cat('\n\ngarch not windowed', cl)
  sd <- calc_sd_garch(data_stock, start_win, FALSE)
  PV_garch_not[i] <- analytics.p_value(sd, cl)
}



