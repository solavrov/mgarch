

#' make list of stocks
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
make_los <- function(...) {

  x <- list(...)

  com_dates <- c()
  for (e in x) {
    if (is.null(com_dates))
      com_dates <- e$Date
    else
      com_dates <- intersect(com_dates, e$Date)
  }

  for (i in 1:length(x)) {
    j <- x[[i]]$Date %in% com_dates
    x[[i]] <- x[[i]][j, ]
    rownames(x[[i]]) <- c()
  }

  return (x)

}


#' calculate covarience matrix for given los
#'
#' @param los
#' @param date
#' @param calc_win
#'
#' @return
#' @export
#'
#' @examples
calc_cov_simple <- function(los, date, calc_win) {
  j <- which(los[[1]]$Date == date)
  k <- (j - calc_win):(j - 1)
  r_matrix <- NULL
  for (e in los) r_matrix <- cbind(r_matrix, e$R[k])
  cov_matrix <- cov(r_matrix)
  return (cov_matrix)
}


#' get vector of stocks returns by date
#'
#' @param los
#' @param date
#'
#' @return
#' @export
#'
#' @examples
get_r_vector <- function(los, date) {
  i <- which(los[[1]]$Date == date)
  r_vector <- c()
  for (j in 1:length(los)) {
    r_vector[j] <- los[[j]]$R[i]
  }
  return (cbind(r_vector))
}


#' calculate sd for portfolio using simple method
#'
#' @param los
#' @param w
#' @param calc_win
#'
#' @return
#' @export
#'
#' @examples
calc_sd_port_simple <- function(los, w, calc_win) {

  w <- cbind(w)

  n <- nrow(los[[1]]) - calc_win

  SD <- c()
  R <- c()
  for (i in 1:n) {
    d <- los[[1]]$Date[calc_win + i]
    cov_matrix <- calc_cov_simple(los, d, calc_win)
    SD[i] <- sqrt(t(w) %*% cov_matrix %*% w)
    R[i] <- t(w) %*% get_r_vector(los, d)
  }

  j <- (calc_win + 1):(calc_win + n)

  df <- data.frame(los[[1]]$Date[j],
                   R,
                   SD)
  colnames(df) <- c("Date", "R", "SD")

  return (df)

}

