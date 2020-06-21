options("max.print"=50)
library(mgarch)

sber <- load_stock_data_by_ticker('sber')
gazp <- load_stock_data_by_ticker('gazp')
lkoh <- load_stock_data_by_ticker('lkoh')
chmf <- load_stock_data_by_ticker('chmf')
gmkn <- load_stock_data_by_ticker('gmkn')

l <- make_los(sber, gazp, lkoh, chmf, gmkn)

w <- c(0.3, 0.3, 0.2, 0.1, 0.1)

# d <- l[[1]]$Date[2000]

# cv_simple <- calc_cov_simple(l, d, calc_win=250)
# cv_gogarch <- gogarch.calc_cov(l, d, calc_win=1000, df=5)


# sd_simple <- calc_sd_port_simple(l, w, calc_win = 250)
sd_gogarch <- gogarch.calc_sd_port(l, w, calc_win = 1000)
cov_data <- gogarch.calc_cov_data(l, calc_win = 1000)
sd_gogarch_2 <- gogarch.calc_sd(cov_data, w)

# print_var_analytics(sd_simple, tail_len = 2200)
# print_var_analytics(sd_gogarch, tail_len = 2200)
#
# draw_var_analytics(sd_simple, tail_len = 2200)
# draw_var_analytics(sd_gogarch, tail_len = 2200)
