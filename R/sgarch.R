
#' sample of student distribution
#'
#' @param n
#' @param nu
#' @param sd
#'
#' @return
#' @export
#'
#' @examples
sgarch.rt <- function(n, nu, sd=1) {
  return (sd * sqrt(1 - 2 / nu) * rt(n, nu) )
}


#' inverted probability distribution for student
#'
#' @param p
#' @param nu
#' @param sd
#'
#' @return
#' @export
#'
#' @examples
sgarch.qt <- function(p, nu, sd=1) {
  return (sd * sqrt(1 - 2 / nu) * qt(p, nu))
}


#' make garch process sample for student
#'
#' @param a0
#' @param a1
#' @param b
#' @param nu
#' @param n
#' @param chisq_df
#'
#' @return
#' @export
#'
#' @examples
sgarch.make_process_sample <- function(a0, a1, b, nu, n, chisq_df=3) {
  h <- rchisq(1, chisq_df) / chisq_df
  r <- sgarch.rt(1, nu, sqrt(h))
  for (i in 2:n) {
    h[i] <- a0 + a1 * r[i-1]^2 + b * h[i-1]
    r[i] <- sgarch.rt(1, nu, sqrt(h[i]))
  }
  return (list(r=r, h=h))
}


#' calculate reversed likelihood for student distribution
#'
#' @param x
#' @param r
#' @param nu
#'
#' @return
#' @export
#'
#' @examples
sgarch.get_rlh <- function(x, r, nu) {
  a0 <- x[1]
  a1 <- x[2]
  b <- x[3]
  h <- tail(get_h_garch(r, a0, a1, b),-1)
  r <- tail(r,-1)
  rlh <- sum(log(h) + (nu + 1) * log(1 + r ^ 2 / h / (nu - 2)))
  return (rlh)
}



#' calculate reversed likelihood for student distribution for a0=1-a1-b case
#'
#' @param x
#' @param r
#' @param nu
#'
#' @return
#' @export
#'
#' @examples
sgarch.get_rlh_2 <- function(x, r, nu) {
  a <- x[1]
  b <- x[2]
  h <- tail(get_h_garch(r, 1-a-b, a, b), -1)
  r <- tail(r, -1)
  rlh <- sum(log(h) + (nu + 1) * log(1 + r ^ 2 / h / (nu - 2)))
  return (rlh)
}


#' rough estimate garch params student
#'
#' @param r
#' @param nu
#' @param a0
#' @param a1
#' @param b
#' @param n_desize
#' @param timer_on
#'
#' @return
#' @export
#'
#' @examples
sgarch.find_param_estim <- function(r,
                                    nu,
                                    a0 = c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 4, 5),
                                    a1 = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                                    b = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99),
                                    n_desize = 0,
                                    timer_on = FALSE) {
  if (timer_on)
    t1 <- Sys.time()
  a0 <- desize(a0, n_desize)
  a1 <- desize(a1, n_desize)
  b <- desize(b, n_desize)
  min_rlh <- sgarch.get_rlh(c(a0[1], a1[1], b[1]), r, nu)
  min_params <- c(a0[1], a1[1], b[1])
  for (k in 1:length(b)) {
    for (j in 1:length(a1)) {
      for (i in 1:length(a0)) {
        if ((i == 1 & j == 1 & k == 1) | a1[j] + b[k] >= 1)
          next
        rlh <- sgarch.get_rlh(c(a0[i], a1[j], b[k]), r, nu)
        if (rlh < min_rlh) {
          min_rlh <- rlh
          min_params <- c(a0[i], a1[j], b[k])
        }
      }
    }
  }
  if (timer_on) {
    t2 <- Sys.time()
    print(t2 - t1)
  }
  return (min_params)
}


#' rough estimate garch params student for a0=1-a1-b case
#'
#' @param r
#' @param nu
#' @param a
#' @param b
#' @param n_desize
#' @param timer_on
#'
#' @return
#' @export
#'
#' @examples
sgarch.find_param_estim_2 <- function(r,
                                      nu,
                                      a = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                                      b = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99),
                                      n_desize = 0,
                                      timer_on = FALSE) {
  if (timer_on)
    t1 <- Sys.time()
  a <- desize(a, n_desize)
  b <- desize(b, n_desize)
  min_rlh <- sgarch.get_rlh_2(c(a[1], b[1]), r, nu)
  min_params <- c(a[1], b[1])
  for (k in 1:length(b)) {
    for (j in 1:length(a)) {
      if ((j == 1 & k == 1) | a[j] + b[k] >= 1)
        next
      rlh <- sgarch.get_rlh_2(c(a[j], b[k]), r, nu)
      if (rlh < min_rlh) {
        min_rlh <- rlh
        min_params <- c(a[j], b[k])
      }
    }
  }
  if (timer_on) {
    t2 <- Sys.time()
    print(t2 - t1)
  }
  return (min_params)
}

#' find garch params estimate for student
#'
#' @param r
#' @param init_par
#' @param nu
#' @param a0_min
#' @param a1_min
#' @param b_min
#' @param a0_max
#' @param a1_max
#' @param b_max
#' @param max_persist
#'
#' @return
#' @export
#'
#' @examples
sgarch.find_params <- function(r,
                               nu,
                               init_par = sgarch.find_param_estim(r, nu),
                               a0_min = 0,
                               a1_min = 0,
                               b_min = 0,
                               a0_max = 5,
                               a1_max = 1,
                               b_max = 1,
                               max_persist = 0.999) {
  g <- list()
  optim_result <- constrOptim(
    init_par,
    sgarch.get_rlh,
    grad = NULL,
    ui = rbind(
      c(1, 0, 0),
      c(0, 1, 0),
      c(0, 0, 1),
      c(-1, 0, 0),
      c(0,-1, 0),
      c(0, 0,-1),
      c(0,-1,-1)
    ),
    ci = c(a0_min,
           a1_min,
           b_min, -a0_max, -a1_max, -b_max, -max_persist),
    mu = 1e-04,
    control = list(),
    method = "Nelder-Mead",
    outer.iterations = 100,
    outer.eps = 1e-05,
    r,
    nu,
    hessian = FALSE
  )
  g$a0 <- optim_result$par[1]
  g$a1 <- optim_result$par[2]
  g$b <- optim_result$par[3]
  return (g)
}


#' find garch params estimate for student for a0=1-a1-b case
#'
#' @param r
#' @param nu
#' @param init_par
#' @param a_min
#' @param b_min
#' @param a_max
#' @param b_max
#'
#' @return
#' @export
#'
#' @examples
sgarch.find_params_2 <- function(r,
                                 nu,
                                 init_par = sgarch.find_param_estim_2(r, nu),
                                 a_min = 0,
                                 b_min = 0,
                                 a_max = 1,
                                 b_max = 1) {
  g <- list()
  optim_result <- constrOptim(
    init_par,
    sgarch.get_rlh_2,
    grad = NULL,
    ui = rbind(c(1, 0),
               c(0, 1),
               c(-1, 0),
               c(0, -1),
               c(-1, -1)),
    ci = c(a_min,
           b_min,-a_max,-b_max,-1),
    mu = 1e-04,
    control = list(),
    method = "Nelder-Mead",
    outer.iterations = 100,
    outer.eps = 1e-05,
    r,
    nu,
    hessian = FALSE
  )
  g$a <- optim_result$par[1]
  g$b <- optim_result$par[2]
  return (g)
}


#' calculate sd using garch for student
#'
#' @param data_stock
#' @param nu
#' @param calc_win
#' @param is_windowed
#' @param ticker
#'
#' @return
#' @export
#'
#' @examples
sgarch.calc_sd <- function(data_stock, nu, calc_win=1000, is_windowed=TRUE, ticker="noname") {

  n <- length(data_stock$Date) - calc_win

  sd <- c()
  a0 <- c()
  a1 <- c()
  b <- c()

  if (is_windowed) comment <- "window" else comment <- ''
  cat("\n\nDoing", ticker, "\n")
  cat("student\'s garch", comment, "progress: 0% ")
  progress <- 1
  for (i in 1:n) {
    if (is_windowed) k <- i else k <- 1
    r <- data_stock$R[k:(calc_win + i - 1)]
    if (i==1) {
      g <- sgarch.find_params(r, nu)
    } else {
      g <- sgarch.find_params(r, nu, c(a0[i-1], a1[i-1], b[i-1]))
    }
    a0[i] <- g$a0
    a1[i] <- g$a1
    b[i] <- g$b
    r <- data_stock$R[k:(calc_win + i)]
    h <- get_h_garch(r, g$a0, g$a1, g$b, h1=var(head(r, -1)))
    sd[i] <- sqrt(tail(h, 1))
    if (i/n * 100 > progress) {
      cat(progress,"% ", sep = '')
      progress <- progress + 1
    }
  }

  j <- (calc_win + 1):length(data_stock$Date)

  df <- data.frame(data_stock$Date[j],
                   data_stock$R[j],
                   sd,
                   a0,
                   a1,
                   b
  )
  colnames(df) <- c("Date", "R", "SD", "a0", "a1", "b")

  return (df)

}
