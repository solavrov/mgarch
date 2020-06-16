library(fGarch)

DIR <- "C:/R/var/"
FILE_CHNG <- "sber_chng.csv"

r <- read.csv(paste0(DIR,FILE_CHNG))
r250 <- head(r$PX_CHNG, 250)
g <- garchFit(~garch(1,1), data = r250, cond.dist = 'std', include.mean = FALSE,
              include.shape=FALSE, include.skew=FALSE)

a0 <- g@fit$params$params$omega
a1 <- g@fit$params$params$alpha1
b <- g@fit$params$params$beta1
Vyt11)