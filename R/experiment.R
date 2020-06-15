

#' run experiment
#'
#' @param stock_tickers
#' @param is_simple
#' @param is_garch
#' @param is_garch_windowed
#' @param is_garch_fixed
#' @param is_lambda
#' @param is_lambda_fixed
#' @param calc_win
#'
#' @return
#' @export
#'
#' @examples
run_experiment <- function(stock_tickers,
                           is_simple=FALSE,
                           is_garch=FALSE,
                           is_garch_windowed=FALSE,
                           is_garch_fixed=FALSE,
                           is_lambda=FALSE,
                           is_lambda_fixed=FALSE,
                           calc_win=250
) {

  for (t in stock_tickers) {

    dir_stock <- paste0(DIR,t,"/")
    data_stock <- load_stock_data_by_ticker(t)
    plot_div_gap_analytics(data_stock, ticker = t)

    save_data(data_stock, paste0(t, "_data_stock.csv"), dir=dir_stock)

    if(is_simple) {
      data_sd_simple <- calc_sd_simple(data_stock, calc_win = calc_win)
      print_var_analytics(data_sd_simple, ticker = t, type="simple")
      save_data(data_sd_simple,
                paste0(t, "_sd_simple_", calc_win, ".csv"), dir=dir_stock)
    }

    if(is_garch) {
      data_sd_garch <- calc_sd_garch(data_stock, ticker = t, calc_win = calc_win)
      print_var_analytics(data_sd_garch, ticker = t, type="garch")
      save_data(data_sd_garch,
                paste0(t, "_sd_garch_", calc_win, ".csv"), dir=dir_stock)
    }

    if(is_garch_windowed) {
      data_sd_garch <- calc_sd_garch(data_stock, ticker = t,
                                     is_windowed = TRUE, calc_win = calc_win)
      print_var_analytics(data_sd_garch, ticker = t, type="garch_windowed")
      save_data(data_sd_garch,
                paste0(t, "_sd_garch_windowed_", calc_win, ".csv"), dir=dir_stock)
    }

    if(is_garch_fixed) {
      g <- get_garch_params(data_stock$R)
      data_sd_garch_fixed <- calc_sd_garch_fixed(data_stock, g$a0, g$a1, g$b,
                                                 calc_win = calc_win)
      print_var_analytics(data_sd_garch_fixed, ticker = t, type="garch_fixed")
      save_data(data_sd_garch_fixed,
                paste0(t, "_sd_garch_fixed_", calc_win, ".csv"), dir=dir_stock)
    }

    if(is_lambda) {
      data_sd_lambda <- calc_sd_lambda(data_stock, ticker = t, calc_win = calc_win)
      print_var_analytics(data_sd_lambda, ticker = t, type="lambda")
      save_data(data_sd_lambda,
                paste0(t, "_sd_lambda_", calc_win, ".csv"), dir=dir_stock)
    }

    if(is_lambda_fixed) {
      l <- get_lambda(data_stock$R)
      data_sd_lambda_fixed <- calc_sd_lambda_fixed(data_stock, lambda = l,
                                                   calc_win = calc_win)
      print_var_analytics(data_sd_lambda_fixed, ticker = t, type="lambda_fixed")
      save_data(data_sd_lambda_fixed,
                paste0(t, "_sd_lambda_fixed_v", calc_win, ".csv"), dir=dir_stock)
    }

  }

  cat("\n\nDone!")

}
