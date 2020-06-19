
#' plot dividend gaps
#'
#' @param data_stock
#' @param ticker
#'
#' @return
#' @export
#'
#' @examples
plot_div_gap_analytics <- function(data_stock, ticker="noname") {
  j <- which(data_stock$DIV != 0)
  plot(data_stock$DIV[j], data_stock$DIV_GAP[j], main=ticker, xlab="div", ylab="gap")
}



#' print var analytics
#'
#' @param data_sd
#' @param conf_level_var
#' @param conf_level_backtest
#' @param tail_len
#' @param ticker
#' @param type
#'
#' @return
#' @export
#'
#' @examples
print_var_analytics <- function(data_sd,
                                conf_level_var=0.99 ,
                                conf_level_backtest=0.9,
                                tail_len=length(data_sd$R),
                                ticker="noname", type="none") {

  data_sd <- tail(data_sd, tail_len)
  n <- tail_len

  alpha <- qnorm(1 - conf_level_var)
  VAR <- alpha * data_sd$SD

  num_of_fails <- sum(as.numeric(data_sd$R < VAR))
  rate <- num_of_fails / n

  mean_fails <- (1 - conf_level_var) * n
  sd_fails <- sqrt(conf_level_var * (1 - conf_level_var) * n)
  alpha_backtest <- qnorm(1 - conf_level_backtest)
  conf_interval_backtest <- round(c(
    mean_fails + alpha_backtest * sd_fails,
    mean_fails - alpha_backtest * sd_fails
  ))

  cat("\n\nStock ticker:", ticker, "\n")
  cat("Analytics of", type, "method\n")
  cat("VAR at confidence", conf_level_var, "\n")
  cat("Cofidence interval at level", conf_level_backtest,"=",
      conf_interval_backtest[1],
      "...",
      round(mean_fails),
      "...",
      conf_interval_backtest[2], "\n")
  cat("Number of fails =", num_of_fails, ", Rate =", rate, "\n")

}



#' draw var analytics
#'
#' @param data_sd
#' @param conf_level_var
#' @param ticker
#' @param type
#'
#' @return
#' @export
#'
#' @examples
draw_var_analytics <- function(data_sd,
                               conf_level_var=0.99 ,
                               tail_len=length(data_sd$R),
                               ticker="noname", type="none") {

  data_sd <- tail(data_sd, tail_len)
  alpha <- qnorm(1 - conf_level_var)
  varisk <- data_sd$SD * alpha
  ymin <- min(data_sd$R, varisk)
  plot(data_sd$Date, varisk, type = 'l', col="green", ylim = c(ymin, 0),
       xlab="date", ylab="price change %")
  lines(data_sd$Date, data_sd$R, type = 'h', col="blue")
  j <- which(data_sd$R < varisk)
  lines(data_sd$Date[j], data_sd$R[j], col="red", lty=0, type='o')

}

