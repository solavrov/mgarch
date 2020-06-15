
#' calculate simple sd
#'
#' @param data_stock
#' @param calc_win
#'
#' @return
#' @export
#'
#' @examples
calc_sd_simple <- function(data_stock, calc_win=250) {

  n <- length(data_stock$Date) - calc_win

  SD <- c()
  for (i in 1:n) {
    j <- i:(calc_win + i - 1)
    SD[i] <- sd(data_stock$R[j])
  }

  j <- (calc_win + 1):length(data_stock$Date)

  df <- data.frame(data_stock$Date[j],
                   data_stock$R[j],
                   SD)
  colnames(df) <- c("Date", "R", "SD")

  return (df)

}
