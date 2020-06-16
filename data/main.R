options('stringsAsFactors'=FALSE)
DIR <- "C:/R/var/"
options("max.print"=50)

#function
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

#function
load_stock_data_by_ticker <- function(ticker, dir=DIR, div_factor=1) {
  file_px <- paste0(ticker, '/', ticker, "_price.csv")
  file_div <- paste0(ticker, '/', ticker, "_div.csv")
  df <- load_stock_data(file_px, file_div, dir, div_factor)
  return (df)
}


# function
select_div_events <- function(data_stock) {
  j <- which(data_stock$DIV != 0)
  df <- data_stock[sort(c(j,j-1)),]
  return (df)
}

# function
plot_div_gap_analytics <- function(data_stock, ticker="noname") {
  j <- which(data_stock$DIV != 0)
  plot(data_stock$DIV[j], data_stock$DIV_GAP[j], main=ticker, xlab="div", ylab="gap")
}

# function
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

# function
get_h_garch <- function(r, a0, a1, b, h1=var(r)) {
  h <- c()
  h[1] <- h1
  for (i in 2:length(r)) h[i] <- a0 + a1 * r[i-1] ^ 2 + b * h[i-1]
  return (h)
}

# function
get_rlh_garch <- function(x, r) {
  a0 <- x[1]
  a1 <- x[2]
  b <- x[3]
  h <- tail(get_h_garch(r, a0, a1, b), -1)
  r <- tail(r, -1)
  return (sum(r^2/h + log(h)))
}

# function
get_h_lambda <- function(r, lambda, h1=var(r)) {
  h <- c()
  h[1] <- h1
  for (i in 2:length(r)) h[i] <- (1 - lambda) * r[i-1] ^ 2 + lambda * h[i-1]
  return (h)
}

# function
get_rlh_lambda <- function(x, r) {
  h <- tail(get_h_lambda(r, x), -1)
  r <- tail(r, -1)
  return (sum(r^2/h + log(h)))
}

#function
desize <- function(v, n=1) {
  if (n>=1) iter <- 1:n else iter <- c()
  for (j in iter) {
    v2 <- c()
    for (i in 1:(length(v)-1))
      v2[i] <- (v[i] + v[i+1]) / 2
    v <- sort(c(v, v2))
  }
  return (v)
}

# function
get_garch_estim <- function(r, 
                            a0=c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 4, 5), 
                            a1=c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
                            b=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99), 
                            n_desize=0,
                            timer_on=FALSE) {
  if (timer_on) t1 <- Sys.time()
  a0 <- desize(a0, n_desize)
  a1 <- desize(a1, n_desize)
  b <- desize(b, n_desize)
  min_rlh <- get_rlh_garch(c(a0[1], a1[1], b[1]), r)
  min_params <- c(a0[1], a1[1], b[1])
  for (k in 1:length(b)) {
    for (j in 1:length(a1)) {
      for (i in 1:length(a0)) {
        if ((i==1 & j==1 & k==1) | a1[j]+b[k]>=1) next
        rlh <- get_rlh_garch(c(a0[i], a1[j], b[k]), r)
        if (rlh < min_rlh) {
          min_rlh <- rlh
          min_params <- c(a0[i], a1[j], b[k])
        }
      }
    }
  }
  if (timer_on) {
    t2 <- Sys.time()
    print(t2-t1)
  }
  return (min_params)
}


# function
get_garch_params <- function(r, 
                             init_par=get_garch_estim(r), 
                             a0_min=0,
                             a1_min=0,
                             b_min=0,
                             a0_max=5,
                             a1_max=1,
                             b_max=1,
                             max_persist=0.999) {
  g <- list()
  optim_result <- constrOptim(init_par, get_rlh_garch, 
                              grad = NULL, 
                              ui = rbind(c(1,0,0), 
                                         c(0,1,0), 
                                         c(0,0,1), 
                                         c(-1,0,0), 
                                         c(0,-1,0), 
                                         c(0,0,-1),
                                         c(0,-1,-1)), 
                              ci = c(a0_min,
                                     a1_min,
                                     b_min,
                                     -a0_max,
                                     -a1_max,
                                     -b_max,
                                     -max_persist), 
                              mu = 1e-04, 
                              control = list(),
                              method = "Nelder-Mead", 
                              outer.iterations = 100, 
                              outer.eps = 1e-05, 
                              r,
                              hessian = FALSE)
  g$a0 <- optim_result$par[1]
  g$a1 <- optim_result$par[2]
  g$b <- optim_result$par[3]
  return (g)
}

# function
get_lambda <- function(r, init_lambda=0.95) {
  optim_result <- optimize(get_rlh_lambda, c(0,1), r)
  return (optim_result$minimum)
}

# function
calc_sd_garch <- function(data_stock, calc_win=250, is_windowed=FALSE, ticker="noname") {
  
  n <- length(data_stock$Date) - calc_win
  
  SD <- c()
  a0 <- c()
  a1 <- c()
  b <- c()
  
  if (is_windowed) comment <- "window" else comment <- ''
  cat("\n\nDoing", ticker, "\n")
  cat("garch", comment, "progress: 0% ")
  progress <- 1
  for (i in 1:n) {
    if (is_windowed) k <- i else k <- 1
    r <- data_stock$R[k:(calc_win + i - 1)]
    if (i==1) {
      g <- get_garch_params(r)
    } else {
      g <- get_garch_params(r, c(a0[i-1], a1[i-1], b[i-1]))
    }
    a0[i] <- g$a0
    a1[i] <- g$a1
    b[i] <- g$b
    r <- data_stock$R[k:(calc_win + i)]
    h <- get_h_garch(r, g$a0, g$a1, g$b, h1=var(head(r, -1)))
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
                   a0,
                   a1,
                   b
  )
  colnames(df) <- c("Date", "R", "SD", "a0", "a1", "b")
  
  return (df)
  
}

# fuction
calc_sd_garch_fixed <- function(data_stock, a0, a1, b, calc_win=250) {
  
  n <- length(data_stock$Date) - calc_win
  h <- get_h_garch(data_stock$R, a0, a1, b)
  SD <- sqrt(h)
  j <- (calc_win + 1):length(data_stock$Date)
  df <- data.frame(data_stock$Date[j], 
                   data_stock$R[j], 
                   SD[j],
                   rep(a0, n),
                   rep(a1, n),
                   rep(b, n)
  )
  colnames(df) <- c("Date", "R", "SD", "a0", "a1", "b")
  
  return (df)
  
}

# function
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

# function
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

# function
print_var_analytics <- function(data_sd, 
                                conf_level_var=0.99 , 
                                conf_level_backtest=0.9, 
                                ticker="noname", type="none") {
  
  n <- length(data_sd$Date)
  alpha <- qnorm(1 - conf_level_var)
  VAR <- alpha * data_sd$SD
  # VAR_ABS <- VAR + mean(data_sd$R)
  
  num_of_fails <- sum(as.numeric(data_sd$R < VAR))
  # num_of_fails_abs <- sum(as.numeric(data_sd$R < VAR_ABS))
  rate <- num_of_fails / n
  # rate_abs <- num_of_fails_abs / n
  
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
  cat("Cofidence interval at level",conf_level_backtest,"=",
      conf_interval_backtest[1], "...", conf_interval_backtest[2], "\n")
  cat("Number of fails =", num_of_fails, ", Rate =", rate, "\n")
  # cat("Number of fails absolute =", num_of_fails_abs, ", Rate =", rate_abs, "\n")
  
}

# function
save_data <- function(data_sd, file_name, dir=DIR) {
  write.csv(data_sd, paste0(dir, file_name), row.names=FALSE)
}

# function
load_df <- function(file_name, dir=DIR) {
  df <- read.csv(paste0(dir, file_name))
  df$Date <- as.Date(df$Date, "%Y-%m-%d")
  return (df)
}

# function
load_sd <- function(ticker, method, calc_win, dir=DIR) {
  file <- paste0(ticker, '/', ticker, '_sd_', method, '_', calc_win, ".csv")
  df <-  load_df(file, dir)
  return (df)
}

# function
draw_var_analytics <- function(data_sd, 
                               conf_level_var=0.99 , 
                               ticker="noname", type="none") {
  alpha <- qnorm(1 - conf_level_var)
  varisk <- data_sd$SD * alpha
  ymin <- min(data_sd$R, varisk)
  plot(data_sd$Date, varisk, type = 'l', col="green", ylim = c(ymin, 0), 
       xlab="date", ylab="price change %")
  lines(data_sd$Date, data_sd$R, type = 'h', col="blue")
  j <- which(data_sd$R < varisk)
  lines(data_sd$Date[j], data_sd$R[j], col="red", lty=0, type='o')
}


# function
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


# function los - list of stocks
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


# function
calc_cov_simple <- function(los, date, calc_win) {
  j <- which(los[[1]]$Date == date)
  k <- (j - calc_win):(j - 1)
  r_matrix <- NULL
  for (e in los) r_matrix <- cbind(r_matrix, e$R[k])
  cov_matrix <- cov(r_matrix)
  return (cov_matrix)
}

#function
get_r_vector <- function(los, date) {
  i <- which(los[[1]]$Date == date)
  r_vector <- c()
  for (j in 1:length(los)) {
    r_vector[j] <- los[[j]]$R[i]
  }
  return (cbind(r_vector))
}

# function
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



# prog
calc_win <- 1500

sber <- load_stock_data_by_ticker("sber")
lkoh <- load_stock_data_by_ticker("lkoh")
chmf <- load_stock_data_by_ticker("chmf")
gazp <- load_stock_data_by_ticker("gazp")
gmkn <- load_stock_data_by_ticker("gmkn")

los1 <- make_los(sber, lkoh)
los2 <- make_los(sber, lkoh, chmf)
los3 <- make_los(sber, lkoh, chmf, gazp)
los4 <- make_los(sber, lkoh, chmf, gazp, gmkn)

sd1 <- calc_sd_port_simple(los1, rep(1,2)/2, calc_win)
sd2 <- calc_sd_port_simple(los2, rep(1,3)/3, calc_win)
sd3 <- calc_sd_port_simple(los3, rep(1,4)/4, calc_win)
sd4 <- calc_sd_port_simple(los4, rep(1,5)/5, calc_win)

print_var_analytics(sd1, ticker = "sber, lkoh")
print_var_analytics(sd2, ticker = "sber, lkoh, chmf")
print_var_analytics(sd3, ticker = "sber, lkoh, chmf, gazp")
print_var_analytics(sd4, ticker = "sber, lkoh, chmf, gazp, gmkn")


# main program
# stock_tickers <- c("sber", "chmf", "gazp", "gmkn", "lkoh")
stock_tickers <- c("sber")
run_experiment(stock_tickers,
               is_simple = TRUE,
               is_garch = FALSE,
               is_garch_fixed = FALSE,
               is_garch_windowed = FALSE,
               calc_win = 2000)


# t <- c("sber")
# file_px <- paste0(t, "_price.csv")
# file_div <- paste0(t, "_div.csv")
# dir_stock <- paste0(DIR,t,"/")
# data_stock <- load_stock_data(file_px, file_div, dir=dir_stock)
# r1 <- data_stock$R[3504:4753]
# r2 <- data_stock$R[3507:4756]
# e1 <- get_garch_estim(r1)
# e2 <- get_garch_estim(r2)
# p1 <- unlist(get_garch_params(r1), use.names = FALSE)
# p2 <- unlist(get_garch_params(r2), use.names = FALSE)



# for (t in stock_tickers) {
#   d <- load_data(paste0(t, "_sd_garch_windowed_1250.csv"), t)
#   print_var_analytics(d, ticker = t)
# }
  
  





