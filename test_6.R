options("max.print"=50)
library(mgarch)

n <- 10^8
t <- 10

a0 <- 0.1
a1 <- 0.8
b <- 0.1

cnt <- 1
x <- c()
for (i in 1:n) {
  x[i] <- sum(garch.make_process_sample(a0, a1, b, n=t, h0=1)$r)
  if (cnt < i / n * 100) {
    cat(paste0(cnt, '% '))
    cnt <- cnt + 1
  }
}


hist(x, breaks = 1000)
