

#' plot dividend gaps
#'
#' @param data_stock
#' @param ticker
#'
#' @return
#' @export
#'
#' @examples
plot_div_gap_analytics <- function(data_stock, ticker = "noname") {
  j <- which(data_stock$DIV != 0)
  plot(
    data_stock$DIV[j],
    data_stock$DIV_GAP[j],
    main = ticker,
    xlab = "div",
    ylab = "gap"
  )
}



#' calculate p-value - prob to get given extrim if model correct
#'
#' @param data_sd
#' @param conf_level_var
#' @param tail_len
#'
#' @return
#' @export
#'
#' @examples
analytics.p_value <- function(data_sd,
                              conf_level_var = 0.99,
                              tail_len = length(data_sd$R)) {

  data_sd <- tail(data_sd, tail_len)
  n <- tail_len

  alpha <- qnorm(1 - conf_level_var)
  VAR <- alpha * data_sd$SD

  num_of_fails <- sum(as.numeric(data_sd$R < VAR))
  rate <- num_of_fails / n

  mean_fails <- (1 - conf_level_var) * n
  sd_fails <- sqrt(conf_level_var * (1 - conf_level_var) * n)
  p_left <- pnorm(num_of_fails, mean = mean_fails, sd = sd_fails)
  p_value <- 2 * min(p_left, 1 - p_left)

  return (p_value)

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
                                conf_level_var = 0.99,
                                conf_level_backtest = 0.9,
                                tail_len = length(data_sd$R),
                                ticker = "noname",
                                type = "none") {

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
  p_left <- pnorm(num_of_fails, mean=mean_fails, sd=sd_fails)
  p_value <- 2 * min(p_left, 1 - p_left)

  cat("\n\nStock ticker:", ticker, "\n")
  cat("Analytics of", type, "method\n")
  cat("VAR at confidence", conf_level_var, "\n")
  cat(
    "Cofidence interval at level",
    conf_level_backtest,
    "=",
    conf_interval_backtest[1],
    "...",
    round(mean_fails),
    "...",
    conf_interval_backtest[2],
    "\n"
  )
  cat("Number of fails =", num_of_fails, ", Rate =", rate, "\n")
  cat("p-value =", p_value, '\n')

}


#' print var analytics for student
#'
#' @param data_sd
#' @param nu
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
print_var_analytics_student <- function(data_sd,
                                        nu,
                                        conf_level_var = 0.99,
                                        conf_level_backtest = 0.9,
                                        tail_len = length(data_sd$R),
                                        ticker = "noname",
                                        type = "none") {

  data_sd <- tail(data_sd, tail_len)
  n <- tail_len

  alpha <- sgarch.qt(1 - conf_level_var, nu)
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
  cat(
    "Cofidence interval at level",
    conf_level_backtest,
    "=",
    conf_interval_backtest[1],
    "...",
    round(mean_fails),
    "...",
    conf_interval_backtest[2],
    "\n"
  )
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
                               conf_level_var = 0.99,
                               tail_len = length(data_sd$R),
                               ticker = "noname",
                               type = "none") {
  data_sd <- tail(data_sd, tail_len)
  alpha <- qnorm(1 - conf_level_var)
  varisk <- data_sd$SD * alpha
  ymin <- min(data_sd$R, varisk)
  plot(
    data_sd$Date,
    varisk,
    type = 'l',
    col = "green",
    ylim = c(ymin, 0),
    xlab = "date",
    ylab = "price change %"
  )
  lines(data_sd$Date, data_sd$R, type = 'h', col = "blue")
  j <- which(data_sd$R < varisk)
  lines(
    data_sd$Date[j],
    data_sd$R[j],
    col = "red",
    lty = 0,
    type = 'o'
  )

}


#' calculate mean correlation from covariance matrix
#'
#' @param cov_matrix
#'
#' @return
#' @export
#'
#' @examples
calc_mean_corr <- function(cov_matrix) {
  n <- nrow(cov_matrix)
  m <- sum(cov2cor(cov_matrix) - diag(n)) / (n ^ 2 - n)
  return (m)
}


#' draw mean correlation from covarience data
#'
#' @param cov_data
#'
#' @return
#' @export
#'
#' @examples
draw_mean_corr <- function(cov_data) {
  d <- as.Date(x = integer(0), origin = "1970-01-01")
  m <- c()
  for (i in 1:length(cov_data)) {
    d[i] <- cov_data[[i]]$Date
    m[i] <- calc_mean_corr(cov_data[[i]]$COV)
  }
  plot(
    d,
    m,
    type = 'l',
    col = 'blue',
    xlab = "date",
    ylab = "mean corr"
  )
}
