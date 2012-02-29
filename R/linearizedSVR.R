##' Train and predict using prototype-based Linearized Support-Vector Regression methods.
##'
##' .. content for \\details{} ..
##'
##' @name LinearizedSVR-package
##' @docType package
##' @title Linearized Support Vector Regression

library(kernlab)
library(LiblineaR)
library(expectreg)



##' .. content for \\description{} (no empty lines) ..
##'
##' .. content for \\details{} ..
##' @title LinearizedSVRTrain
##' @param X matrix of examples, one example per row.
##' @param Y vector of target values.  Must be the same length as the number of rows in \code{X}.
##' @param C cost of constraints violation
##' @param epsilon tolerance of termination criterion for optimization
##' @param nump number of prototypes by which to represent each example in \code{X}
##' @param ktype kernel-generating function, typically from the \pkg{kernlab} package
##' @param kpar a list of any parameters necessary for \code{ktype}.  See Details.
##' @param prototypes the method by which prototypes will be chosen
##' @param clusterY whether to cluster \code{X} and \code{Y} jointly
##' when using \code{prototypes="kmeans"}.  Otherwise \code{X} is
##' clustered without influence from \code{Y}.
##' @param epsilon.up allows you to use a different setting for
##' \code{epsilon} in the positive direction.
##' @param epsilon.down allows you to use a different setting for
##' \code{epsilon} in the negative direction.
##' @param quantile if non-null, do quantile regression using the
##' given quantile value.  Currently uses the \code{expectreg}
##' package.
##' @return a model object that can later be used as the first
##' argument for the \code{predict()} method.
##' @export
##' @examples
##' dat <- rbind(data.frame(y=2, x1=rnorm(500, 1), x2=rnorm(500, 1)),
##'              data.frame(y=1, x1=rnorm(500,-1), x2=rnorm(500,-1)))
##' mod <- LinearizedSVRTrain(X=as.matrix(dat[-1]), Y=dat$y, nump=6)
##' res <- predict(mod, newdata=as.matrix(dat[-1]))
##' plot(x2 ~ x1, dat, col=c("red","green")[1+(res>1.5)], pch=c(3,20)[dat$y])

LinearizedSVRTrain <- function(X, Y,
                C = 1, epsilon = 0.01, nump = floor(sqrt(N)),
                ktype=rbfdot, kpar, prototypes=c("kmeans","random"), clusterY=FALSE,
                epsilon.up=epsilon, epsilon.down=epsilon, quantile = NULL){

  N <- nrow(X); D <- ncol(X)
  tmp <- .normalize(cbind(Y,X))
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

  if(is.null(quantile)){
    svc <- LiblineaR(data, labels, type=2, cost=C, bias = TRUE)
    W <- svc$W
  }
  else{
    ex <- expectreg.ls(Yn~rb(Xt, type="special", B=Xt, P=diag(rep(1, nump))),
                      estimate="bundle", smooth="fixed", expectiles=quantile)
    W <- c(-1, unlist(ex$coefficients), ex$intercept)
  }
  model <- list(W = W, prototypes=prototypes, params=pars, kernel=kernel)
  class(model) <- 'LinearizedSVR'
  return(model)
}

##' Predict method for LinearizedSVR models
##'
##' This method produces predicted value, obtained by evaluating the
##' trained model \code{object} on the given data set \code{newdata}.
##' The columns of \code{newdata} must correspond exactly to the
##' columns of \code{X} when the model \code{object} was created.
##'
##' @title predict
##' @param object a model previously trained using \code{LinearizedSVRTrain()}
##' @param newdata a matrix of new data to run predictions on, with
##' the same columns as \code{X} had during training
##' @param ... further arguments passed to or from other methods
##' @return a vector of predicted regression values, with length equal
##' to the number of rows in \code{newdata}.
##' @method predict LinearizedSVR
##' @S3method predict LinearizedSVR
predict.LinearizedSVR <- function(object, newdata, ...){
  model <- object
  tmp <- .normalize(cbind(0, newdata), model$params) #the zero column is because the params had the target also
  Xn <- tmp$Xn[, -1, drop=FALSE]
  Xt <- kernelMatrix(model$kernel, Xn, model$prototypes)
  Xt <- cbind(Xt, array(1, dim(Xt)[1]))
  wx.b <- Xt %*% model$W[-1] #all but the first weight
  Y.hat <- as.vector(-wx.b / model$W[1])
  Y.hat <- Y.hat *(model$params$MM[1]-model$params$mm[1]) + model$params$mm[1] #unnormalize predictions
  return(Y.hat)
}

.normalize <- function(X, params){
  if (missing(params)){
    params <- list(MM = apply(X, 2, max), mm = apply(X, 2, min))
  }
  Xn <- sweep(X, 2, params$mm)
  Xn <- sweep(Xn, 2, ifelse(params$MM==params$mm, 1, params$MM-params$mm), FUN="/")
  return(list(Xn=Xn, params=params))
}

.unnormalize <- function(X, params){
  Xn <- sweep(X, 2, ifelse(params$MM==params$mm, 1, params$MM-params$mm), FUN=`*`)
  Xn <- sweep(Xn, 2, params$mm, FUN=`+`)
  return(Xn)
}
