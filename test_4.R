options("max.print"=50)
library(mgarch)

sber <- load_stock_data_by_ticker('sber')
gazp <- load_stock_data_by_ticker('gazp')
lkoh <- load_stock_data_by_ticker('lkoh')
chmf <- load_stock_data_by_ticker('chmf')
gmkn <- load_stock_data_by_ticker('gmkn')


sd_sber_std_6 <- sgarch.calc_sd(sber, 6)
sd_gazp_std_6 <- sgarch.calc_sd(gazp, 6)
sd_lkoh_std_6 <- sgarch.calc_sd(lkoh, 6)
sd_chmf_std_6 <- sgarch.calc_sd(chmf, 6)
sd_gmkn_std_6 <- sgarch.calc_sd(gmkn, 6)


for (t in c('sber', 'gazp', 'lkoh', 'chmf', 'gmkn')) {
  sd <- eval(parse(text=paste0('sd_', t, '_std_6_1000')))
  print_var_analytics_student(sd, 6, conf_level_var = 0.99, type = 'student 6', ticker = t)
  print_var_analytics_student(sd, 6, conf_level_var = 0.95, type = 'student 6', ticker = t)
  print_var_analytics_student(sd, 6, conf_level_var = 0.90, type = 'student 6', ticker = t)
}


sd_sber_std_4_1000 <- sgarch.calc_sd(sber, 4, calc_win = 1000)
sd_gazp_std_4_1000 <- sgarch.calc_sd(gazp, 4, calc_win = 1000)
sd_lkoh_std_4_1000 <- sgarch.calc_sd(lkoh, 4, calc_win = 1000)
sd_chmf_std_4_1000 <- sgarch.calc_sd(chmf, 4, calc_win = 1000)
sd_gmkn_std_4_1000 <- sgarch.calc_sd(gmkn, 4, calc_win = 1000)

sd_sber_std_8_1000 <- sgarch.calc_sd(sber, 8, calc_win = 1000)
sd_gazp_std_8_1000 <- sgarch.calc_sd(gazp, 8, calc_win = 1000)
sd_lkoh_std_8_1000 <- sgarch.calc_sd(lkoh, 8, calc_win = 1000)
sd_chmf_std_8_1000 <- sgarch.calc_sd(chmf, 8, calc_win = 1000)
sd_gmkn_std_8_1000 <- sgarch.calc_sd(gmkn, 8, calc_win = 1000)


sd_sber <- calc_sd_garch(sber)
sd_gazp <- calc_sd_garch(gazp)
sd_lkoh <- calc_sd_garch(lkoh)
sd_chmf <- calc_sd_garch(chmf)
sd_gmkn <- calc_sd_garch(gmkn)

for (t in c('sber', 'gazp', 'lkoh', 'chmf', 'gmkn')) {
  sd <- eval(parse(text=paste0('sd_', t)))
  print_var_analytics(sd, conf_level_var = 0.99, type = 'garch', ticker = t)
  print_var_analytics(sd, conf_level_var = 0.95, type = 'garch', ticker = t)
  print_var_analytics(sd, conf_level_var = 0.90, type = 'garch', ticker = t)
}

