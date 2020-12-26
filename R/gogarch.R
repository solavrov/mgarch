
#' make sample of factor processes for go-garch
#'
#' @param a
#' @param b
#' @param n
#' @param chisq_df
#'
#' @return
#' @export
#'
#' @examples
gogarch.make_factor_process_sample <- function(a, b, n, chisq_df=3) {

  X <- c()
  SD2 <- c()
  m <- length(a)

  for (j in 1:m) {
    sd2_j <- rchisq(1, chisq_df) / chisq_df
    x_j <- rnorm(1, sd = sqrt(sd2_j))
    a0 <- (1 - a[j] - b[j])
    for (i in 2:n) {
      sd2_j[i] <- a0 + a[j] * x_j[i-1]^2 + b[j] * sd2_j[i-1]
      x_j[i] <- rnorm(1, sd = sqrt(sd2_j[i]))
    }
    X <- rbind(X, x_j, deparse.level = 0)
    SD2 <- rbind(SD2, sd2_j, deparse.level = 0)
  }

  return (list(X=X, SD2=SD2))

}


#' make sample of list of stocks
#'
#' @param R
#'
#' @return
#' @export
#'
#' @examples
gogarch.make_los_sample <- function(R) {
  m <- nrow(R)
  n <- ncol(R)
  dates <- seq(Sys.Date(), Sys.Date() + n - 1, by=1)
  x <- list()
  for (i in 1:m)
    x[[i]] <- data.frame(Date=dates, R=R[i,])
  return(x)
}


#' find sigma estimate
#'
#' @param r
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_SIGMA_estim <- function(R) {
  return (R %*% t(R) / ncol(R))
}


#' decompose matrix
#'
#' @param M
#'
#' @return
#' @export
#'
#' @examples
gogarch.decompose_matrix <- function(M) {
  e <- eigen(M)
  return (list(U=e$vectors, L=diag(e$values)))
}


#' calculate H2 from H1
#'
#' @param H1
#' @param S
#' @param df
#'
#' @return
#' @export
#'
#' @examples
gogarch.H2 <- function(H1, S, df=5) {

  n <- ncol(S)
  m <- nrow(S)
  inv_H1 <- solve(H1)

  H2 <- matrix(0, m, m)
  for (i in 1:n) {
    g <- c(t(S[,i]) %*% inv_H1 %*% S[,i])
    w <- (m + df) / (df + g)
    H2 <- H2 + w * S[,i] %*% t(S[,i])
  }
  H2 <- H2 / n

  return (H2)

}


#' get square of distance between matrices
#'
#' @param A
#' @param B
#'
#' @return
#' @export
#'
#' @examples
gogarch.get_matrix_distance <- function(A,B) {
  return (sum((A-B)^2)/length(A))
}


#' find H for given accuracy
#'
#' @param S
#' @param H
#' @param df
#' @param acuracy
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_H <- function(S, df=5, accuracy=10^-2) {

  accuracy_square <- accuracy ^ 2

  H <- diag(nrow(S))
  H2 <- gogarch.H2(H, S, df)
  a <- gogarch.get_matrix_distance(H, H2)
  H <- H2

  while (a > accuracy_square) {
    H2 <- gogarch.H2(H, S, df)
    a2 <- gogarch.get_matrix_distance(H, H2)
    if (a2 >= a) print("find_H may not be converging")
    a <- a2
    H <- H2
  }

  return (H)

}

#' find OMEGA
#'
#' @param R
#' @param SIGMA_est
#' @param df
#' @param accuracy
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_OMEGA_estim <- function(R, SIGMA_est, df=5, accuracy=10^-2) {
  decomp <- gogarch.decompose_matrix(SIGMA_est)
  U <- decomp$U
  L <- decomp$L
  S <- sqrt(solve(L)) %*% t(U) %*% R
  H <- gogarch.find_H(S, df=df, accuracy=accuracy)
  decomp <- gogarch.decompose_matrix(H)
  V <- decomp$U
  OMEGA_est <- U %*% sqrt(L) %*% t(V)
  return (OMEGA_est)
}


#' find estimate of X process
#'
#' @param R
#' @param OMEGA_est
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_X_process_estim <- function(R, OMEGA_est) {
  return (solve(OMEGA_est) %*% R)
}


#' find estimates of garch params for factor processes
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_garch_param_estim <- function(X) {
  a <- c()
  b <- c()
  for (i in 1:nrow(X)) {
    gp <- garch.find_params_2(X[i,])
    a[i] <- gp$a
    b[i] <- gp$b
  }
  return (list(a=a, b=b))
}


#' find estimate for conditional covariance matrix for moment t
#'
#' @param X
#' @param gp garch params of X
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_SD2_estim <- function(X, gp) {
  SD2 <- c()
  for (i in 1:nrow(X)) {
    a0 <- 1 - gp$a[i] - gp$b[i]
    SD2 <- rbind(SD2, get_h_garch(X[i,], a0, gp$a[i], gp$b[i], extra_h=TRUE))
  }
  return (SD2)
}


#' find conditional covariance matrix for moment t
#'
#' @param SD2
#' @param OMEGA
#' @param t
#'
#' @return
#' @export
#'
#' @examples
gogarch.find_conditional_SIGMA <- function(SD2, OMEGA, t) {
  return (OMEGA %*% diag(SD2[,t]) %*% t(OMEGA))
}


#' permutate columns of M to get close to M0
#'
#' @param M given matrix for permutation
#' @param M0 target
#'
#' @return
#' @export
#'
#' @examples
gogarch.permutate_colunms <- function(M, M0) {
  S <- t(M) %*% M0
  P <- abs(S)
  I <- matrix(rep(1, length(P)), nrow = nrow(P), ncol = ncol(P))
  while (max(P) >= 0) {
    ind <- which(P==max(P), arr.ind = TRUE)[1,]
    P[ind['row'],] <- rep(-2, ncol(P))
    P[,ind['col']] <- rep(-2, nrow(P))
    P[ind['row'], ind['col']] <- -1
  }
  P <- (P + 2 * I) * sign(S)
  return (M %*% P)
}


#' calculate covariance matrix for given portfolio on given date
#'
#' @param los
#' @param date
#' @param calc_win
#' @param df
#'
#' @return
#' @export
#'
#' @examples
gogarch.calc_cov <- function(los, date, calc_win, df=5) {
  j <- which(los[[1]]$Date == date)
  k <- (j - calc_win):(j - 1)
  R <- c()
  for (e in los) R <- rbind(R, e$R[k])

  SIGMA_est <- gogarch.find_SIGMA_estim(R)
  OMEGA_est <- gogarch.find_OMEGA_estim(R, SIGMA_est, df)
  X_est <- gogarch.find_X_process_estim(R, OMEGA_est)
  gp <- gogarch.find_garch_param_estim(X_est)
  SD2_est <- gogarch.find_SD2_estim(X_est, gp)
  cov_matrix <- gogarch.find_conditional_SIGMA(SD2_est, OMEGA_est, ncol(SD2_est))

  return (cov_matrix)
}


#' calculate gogarch COV for portfolio for moving window for all dates
#'
#' @param los
#' @param w
#' @param calc_win
#'
#' @return
#' @export
#'
#' @examples
gogarch.calc_cov_data <- function(los, calc_win, df=5) {

  n <- nrow(los[[1]]) - calc_win

  cov_data <- list()

  cat("Doing go-garch, progress: 0% ")
  progress <- 1
  for (i in 1:n) {
    d <- los[[1]]$Date[calc_win + i]
    cov_matrix <- gogarch.calc_cov(los, d, calc_win, df)
    cov_data[[i]] <- list(Date=d, COV=cov_matrix, r_vector=get_r_vector(los, d))
    if (i / n * 100 > progress) {
      cat(progress,'% ', sep = '')
      progress <- progress + 1
    }
  }

  return (cov_data)

}


#' calculate gogarch SD for portfolio for moving window form cov_data
#'
#' @param cov_data
#' @param w
#'
#' @return
#' @export
#'
#' @examples
gogarch.calc_sd <- function(cov_data, w) {

  w <- cbind(w)

  Date <- as.Date(x = integer(0), origin = "1970-01-01")
  SD <- c()
  R <- c()

  for (i in 1:length(cov_data)) {
    Date[i] <- cov_data[[i]]$Date
    SD[i] <- sqrt(t(w) %*% cov_data[[i]]$COV %*% w)
    R[i] <- t(w) %*% cov_data[[i]]$r_vector
  }

  df <- data.frame(Date, R, SD)
  colnames(df) <- c("Date", "R", "SD")

  return (df)

}


