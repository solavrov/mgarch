options("max.print"=50)

sber <- load_stock_data_by_ticker('sber')
gazp <- load_stock_data_by_ticker('gazp')
lkoh <- load_stock_data_by_ticker('lkoh')

l <- make_los(sber, gazp, lkoh)

w <- c(0.4, 0.3, 0.3)

d <- l[[1]]$Date[2000]

# cv_simple <- calc_cov_simple(l, d, calc_win=250)
# cv_gogarch <- gogarch.calc_cov(l, d, calc_win=1000, df=5)



sd_simple <- calc_sd_port_simple(l, w, calc_win = 250)
sd_gogarch <- gogarch.calc_sd_port(l, w, calc_win = 1000)
