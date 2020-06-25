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
  return (sum(log(1 + r ^ 2 / h / (nu - 2))))
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
                               init_par = sgarch.find_param_estim(r, nu),
                               nu,
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
