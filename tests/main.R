
##' # Test code (should convert to RUnit):
##' dat <- data.frame(y=2, x1=rnorm(500,1), x2=rnorm(500,1))
##' dat <- rbind(dat, data.frame(y=1, x1=rnorm(500,-1), x2=rnorm(500,-1)))
##' plot(x2 ~ x1, dat, col=c("red","green")[dat$y])
##' mod <- LinearizedSVRTrain(X=as.matrix(dat[-1]), Y=dat$y, nump=6)
##' res <- predict(mod, newdata=as.matrix(dat[-1]))
##' table(res>0)
##' plot(x2 ~ x1, dat, col=c("red","green")[1+(res>0.5)])
##'
##' # An example that goes very badly with clusterY=FALSE, but often
##' # (not always) does well with clusterY=TRUE:
##' dat2 <- data.frame(y=0, x1=rnorm(500,0), x2=rnorm(500,3))
##' dat2 <- rbind(dat2, data.frame(y=0, x1=rnorm(500,0), x2=rnorm(500,-3)))
##' dat2$y <- 1+(dat2$x1 > 0)
##' plot(x2 ~ x1, dat2, col=c("red","green")[dat2$y])
##' mod2 <- LinearizedSVRTrain(X=as.matrix(dat2[-1]), Y=dat2$y, nump=4, clusterY=TRUE)
##' res2 <- predict(mod2, newdata=as.matrix(dat2[-1]))
##' table(res2>1.5)
##' plot(x2 ~ x1, dat2, col=c("red","green")[1+(res2>1.5)])
