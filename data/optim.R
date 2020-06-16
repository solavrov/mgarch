f <- function(v) {
  q <- 1:10
  y <- sum((v - q) ^ 4)
  return (y)
}

v0 <- optim(rep(0,10), f, method = "BFGS")


g <- function(x) {
  return ((x-5)^2)
}

v1 <- optim(0, g, method = "BFGS")

v2 <- optimize(g, c(0,3))


DIR <- "C:/R/var/"
FILE_CHNG <- "sber_chng.csv"

r <- read.csv(paste0(DIR,FILE_CHNG))
r250 <- head(r$PX_CHNG, 250)


h <- function(x,y,z) {
  return ((x-1/2)^2+(y-1/3)^2+(z-1/4)^2)
}


min3 <- function(f, x, y, z) {
  min_f <- f(x[1], y[1], z[1])
  min_xyz <- c(x[1], y[1], z[1])
  for (k in 1:length(z)) {
    for (j in 1:length(y)) {
      for (i in 1:length(x)) {
        if (f(x[i], y[j], z[k]) < min_f) {
          min_f <- f(x[i], y[j], z[k])
          min_xyz <- c(x[i], y[j], z[k])
        }
      }
    }
  }
  return (min_xyz)
}
