
#' global params
#'
#' @export
DIR <- "C:/R/var/"
#options('stringsAsFactors'=FALSE)
#options("max.print"=50)


#' load stock data from file
#'
#' @param file_px
#' @param file_div
#' @param dir
#' @param div_factor
#'
#' @return
#' @export
#'
#' @examples
load_stock_data <- function(file_px, file_div, dir=DIR, div_factor=1) {

  # read and preprocess
  data_px <- read.csv(paste0(dir,file_px))
  data_div <- read.csv(paste0(dir,file_div))
  data_px$Date <- as.Date(data_px$Date, "%d.%m.%Y")
  data_px <- data_px[order(data_px$Date),]
  data_div$Ex.Date <- as.Date(data_div$Ex.Date, "%d.%m.%Y")
  data_div <- data_div[order(data_div$Ex.Date),]
  n <- length(data_px$Date)

  # build dividend vector
  DIV <- c()
  for (i in 1:n) {
    j <- which(data_div$Ex.Date == data_px$Date[i])
    if (length(j) == 0) DIV[i] <- 0 else DIV[i] <- sum(data_div$Dividend[j])
  }

  #build div gap vector
  DIV_GAP <- rep(0, n)
  j <- which(DIV != 0)
  for (i in j) DIV_GAP[i] <- data_px$PX_LAST[i] - data_px$PX_LAST[i-1]

  # build price change vector
  R <- c()
  for (i in 2:length(data_px$Date)) {
    R[i-1] <- log(
      (data_px$PX_LAST[i] + DIV[i] * div_factor) / data_px$PX_LAST[i-1]
    ) * 100
  }

  # build final data frame
  df <- data.frame(tail(data_px$Date, -1),
                   tail(data_px$PX_LAST, -1),
                   tail(DIV, -1),
                   tail(DIV_GAP, -1),
                   R)
  colnames(df) <- c("Date", "PX_LAST", "DIV", "DIV_GAP", "R")

  return (df)

}


#' load stock data by ticket
#'
#' @param ticker
#' @param dir
#' @param div_factor
#'
#' @return
#' @export
#'
#' @examples
load_stock_data_by_ticker <- function(ticker, dir=DIR, div_factor=1) {
  file_px <- paste0(ticker, '/', ticker, "_price.csv")
  file_div <- paste0(ticker, '/', ticker, "_div.csv")
  df <- load_stock_data(file_px, file_div, dir, div_factor)
  return (df)
}



#' select dividend events
#'
#' @param data_stock
#'
#' @return
#' @export
#'
#' @examples
select_div_events <- function(data_stock) {
  j <- which(data_stock$DIV != 0)
  df <- data_stock[sort(c(j,j-1)),]
  return (df)
}


#' save data
#'
#' @param data_sd
#' @param file_name
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
save_data <- function(data, file_name, dir=DIR) {
  write.csv(data, paste0(dir, file_name), row.names=FALSE)
}


#' load data
#'
#' @param file_name
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
load_df <- function(file_name, dir=DIR) {
  df <- read.csv(paste0(dir, file_name))
  df$Date <- as.Date(df$Date, "%Y-%m-%d")
  return (df)
}


#' load sd data by ticker
#'
#' @param ticker
#' @param method
#' @param calc_win
#' @param dir
#'
#' @return
#' @export
#'
#' @examples
load_sd <- function(ticker, method, calc_win, dir=DIR) {
  file <- paste0(ticker, '/', ticker, '_sd_', method, '_', calc_win, ".csv")
  df <-  load_df(file, dir)
  return (df)
}
