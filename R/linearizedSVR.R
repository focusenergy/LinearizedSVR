library(kernlab)
library(LiblineaR)

LinearizedSVRTrain <- function(X, Y,
                C = 1, epsilon = 0.01, nump = floor(sqrt(N)),
                ktype=rbfdot, kpar, prototypes=c("kmeans","random"), clusterY=FALSE,
                epsilon.up=epsilon, epsilon.down=epsilon){

  N <- nrow(X); D <- ncol(X)
  tmp <- normalize(cbind(Y,X))
  Xn <- tmp$Xn[,-1]
  Yn <- tmp$Xn[,1]
  pars <- tmp$params

  prototypes <- switch(match.arg(prototypes),
                       kmeans = if(clusterY) {
                             suppressWarnings(kmeans(tmp$Xn, centers=nump))$centers[,-1]
                           } else {
                             suppressWarnings(kmeans(Xn, centers=nump))$centers
                           },
                       random = Xn[sample(nrow(Xn),nump),])
  rm(tmp, X, Y)  ## Free some memory


  if(!is.na(match("sigma", names(formals(ktype))))){
    if (missing(kpar)) {
      kpar <- list()
    }
    if(is.null(kpar$sigma)){
      kpar$sigma <- median(dist(Xn[sample(nrow(Xn),min(nrow(Xn),50)),]))
    }
  }
  if ('sigma' %in% names(kpar))
    message("Sigma: ", kpar$sigma)

  kernel <- do.call(ktype, kpar)

  Xt <- kernelMatrix(kernel, Xn, prototypes)
  message("Kernel dimensions: [", paste(dim(Xt), collapse=' x '), "]")

  Xt0 <- cbind(Yn-epsilon.down, Xt)
  Xt1 <- cbind(Yn+epsilon.up, Xt)
  data <- rbind(Xt0, Xt1)
  labels <- rep(c(0,1), each=N)

  svc <- LiblineaR(data, labels, type=2, cost=C, bias = TRUE)
  model <- list(W = svc$W, prototypes=prototypes, params=pars, kernel=kernel)
  class(model) <- 'LinearizedSVR'
  return(model)
}


predict.LinearizedSVR <- function(model, newdata){
  tmp <- normalize(cbind(0, newdata), model$params) #the zero column is because the params had the target also
  Xn <- tmp$Xn[, -1, drop=FALSE]
  Xt <- kernelMatrix(model$kernel, Xn, model$prototypes)
  Xt <- cbind(Xt, array(1, dim(Xt)[1]))
  wx.b <- Xt %*% model$W[-1] #all but the first weight
  Y.hat <- as.vector(-wx.b / model$W[1])
  Y.hat <- Y.hat *(model$params$MM[1]-model$params$mm[1]) + model$params$mm[1] #unnormalize predictions
  return(Y.hat)
}

## Old name, for backward compatibility
LinearizedSVRPredict <- predict.LinearizedSVR


normalize <- function(X, params){
  if (missing(params)){
    params <- list(MM = apply(X, 2, max), mm = apply(X, 2, min))
  }
  Xn <- sweep(X, 2, params$mm)
  Xn <- sweep(Xn, 2, ifelse(params$MM==params$mm, 1, params$MM-params$mm), FUN="/")
  return(list(Xn=Xn, params=params))
}

unnormalize <- function(X, params){
  Xn <- sweep(X, 2, ifelse(params$MM==params$mm, 1, params$MM-params$mm), FUN=`*`)
  Xn <- sweep(Xn, 2, params$mm, FUN=`+`)
  return(Xn)
}


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
