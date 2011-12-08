library(kernlab)
library(LiblineaR)

LinearizedSVRTrain <- function(X, Y,
                C = 1, epsilon = 0.01, nump = floor(sqrt(N)),
                ktype=rbfdot, kpar=list(sigma=median(dist(X[sample(1:nrow(X),max(nrow(X),50)),])))){

  N <- dim(X)[1]; D <- dim(X)[2]
  tmp <- normalize(cbind(Y,X))
  Xn <- tmp$Xn[,-1]
  Yn <- tmp$Xn[,1]
  km <- suppressWarnings(kmeans(Xn, centers=nump))
  prototypes <- km$centers
  kernel <- do.call(ktype, kpar)

  print(kpar$sigma)

  Xt <- kernelMatrix(kernel, Xn, prototypes)
  print(dim(Xt))

  Xt0 <- cbind(Yn-epsilon, Xt)
  Xt1 <- cbind(Yn+epsilon, Xt)
  data <- rbind(Xt0, Xt1)
  labels <- append(array(0, N), array(1, N))

  svc <- LiblineaR(data, labels, type=2, cost=C, bias = TRUE)
  model <- list(W = svc$W, prototypes=prototypes, params=tmp$params, kernel=kernel)
  class(model) <- 'LinearizedSVR'
  return(model)
}


predict.LinearizedSVR <- function(model, newdata){
  tmp <- normalize(cbind(array(0, dim(newdata)[1]),newdata), model$params) #the zero array is because the params had the target also
  Xn <- tmp$Xn[,-1]
  Xt <- kernelMatrix(model$kernel, Xn, model$prototypes)
  Xt <- cbind(Xt, array(1, dim(Xt)[1]))
  wx.b <- Xt %*% model$W[-1] #all but the first weight
  Y.hat <- (-wx.b / model$W[1])
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
  Xn <- sweep(Xn, 2, (params$MM-params$mm), FUN="/")
  return(list(Xn=Xn, params=params))
}

