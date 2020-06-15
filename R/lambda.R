
#' calculate variances for given lambda
#'
#' @param r
#' @param lambda
#' @param h1
#'
#' @return
#' @export
#'
#' @examples
get_h_lambda <- function(r, lambda, h1=var(r)) {
  h <- c()
  h[1] <- h1
  for (i in 2:length(r)) h[i] <- (1 - lambda) * r[i-1] ^ 2 + lambda * h[i-1]
  return (h)
}


#' calculate reversed likelihood for given lambda
#'
#' @param x
#' @param r
#'
#' @return
#' @export
#'
#' @examples
get_rlh_lambda <- function(x, r) {
  h <- tail(get_h_lambda(r, x), -1)
  r <- tail(r, -1)
  return (sum(r^2/h + log(h)))
}


#' estimate lambda
#'
#' @param r
#' @param init_lambda
#'
#' @return
#' @export
#'
#' @examples
get_lambda <- function(r, init_lambda=0.95) {
  optim_result <- optimize(get_rlh_lambda, c(0,1), r)
  return (optim_result$minimum)
}


#' calculate sd using lambda method
#'
#' @param data_stock
#' @param calc_win
#' @param ticker
#'
#' @return
#' @export
#'
#' @examples
calc_sd_lambda <- function(data_stock, calc_win=250, ticker="noname") {

  n <- length(data_stock$Date) - calc_win

  SD <- c()
  lambda <- c()

  cat("\n\nDoing", ticker, "\n")
  cat("lambda progress: 0% ")
  progress <- 1
  for (i in 1:n) {
    r <- data_stock$R[1:(calc_win + i - 1)]
    lambda[i] <- get_lambda(r)
    r <- data_stock$R[1:(calc_win + i)]
    h <- get_h_lambda(r, lambda[i], h1=var(head(r, -1)))
    SD[i] <- sqrt(tail(h, 1))
    if (i/n *100 > progress) {
      cat(progress,"% ", sep = '')
      progress <- progress + 1
    }
  }

  j <- (calc_win + 1):length(data_stock$Date)

  df <- data.frame(data_stock$Date[j],
                   data_stock$R[j],
                   SD,
                   lambda
  )
  colnames(df) <- c("Date", "R", "SD", "lambda")

  return (df)

}


#' calculate sd using given fixed lambda
#'
#' @param data_stock
#' @param lambda
#' @param calc_win
#'
#' @return
#' @export
#'
#' @examples
calc_sd_lambda_fixed <- function(data_stock, lambda, calc_win=250) {

  n <- length(data_stock$Date) - calc_win
  h <- get_h_lambda(data_stock$R, lambda)
  SD <- sqrt(h)
  j <- (calc_win + 1):length(data_stock$Date)
  df <- data.frame(data_stock$Date[j],
                   data_stock$R[j],
                   SD[j],
                   rep(lambda, n)
  )
  colnames(df) <- c("Date", "R", "SD", "lambda")

  return (df)

}

