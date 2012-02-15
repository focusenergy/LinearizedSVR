\name{LinearizedSVRTrain}
\alias{LinearizedSVRTrain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LinearizedSVRTrain(X, Y, C = 1, epsilon = 0.01, nump = floor(sqrt(N)), ktype = rbfdot, kpar, prototypes = c("kmeans", "random"), clusterY = FALSE, epsilon.up = epsilon, epsilon.down = epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{nump}{
%%     ~~Describe \code{nump} here~~
}
  \item{ktype}{
%%     ~~Describe \code{ktype} here~~
}
  \item{kpar}{
%%     ~~Describe \code{kpar} here~~
}
  \item{prototypes}{
%%     ~~Describe \code{prototypes} here~~
}
  \item{clusterY}{
%%     ~~Describe \code{clusterY} here~~
}
  \item{epsilon.up}{
%%     ~~Describe \code{epsilon.up} here~~
}
  \item{epsilon.down}{
%%     ~~Describe \code{epsilon.down} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, C = 1, epsilon = 0.01, nump = floor(sqrt(N)), 
    ktype = rbfdot, kpar, prototypes = c("kmeans", "random"), 
    clusterY = FALSE, epsilon.up = epsilon, epsilon.down = epsilon) 
{
    N <- nrow(X)
    D <- ncol(X)
    tmp <- normalize(cbind(Y, X))
    Xn <- tmp$Xn[, -1]
    Yn <- tmp$Xn[, 1]
    pars <- tmp$params
    prototypes <- switch(match.arg(prototypes), kmeans = if (clusterY) {
        suppressWarnings(kmeans(tmp$Xn, centers = nump))$centers[, 
            -1]
    } else {
        suppressWarnings(kmeans(Xn, centers = nump))$centers
    }, random = Xn[sample(nrow(Xn), nump), ])
    rm(tmp, X, Y)
    if (!is.na(match("sigma", names(formals(ktype))))) {
        if (missing(kpar)) {
            kpar <- list()
        }
        if (is.null(kpar$sigma)) {
            kpar$sigma <- median(dist(Xn[sample(nrow(Xn), min(nrow(Xn), 
                50)), ]))
        }
    }
    if ("sigma" \%in\% names(kpar)) 
        message("Sigma: ", kpar$sigma)
    kernel <- do.call(ktype, kpar)
    Xt <- kernelMatrix(kernel, Xn, prototypes)
    message("Kernel dimensions: [", paste(dim(Xt), collapse = " x "), 
        "]")
    Xt0 <- cbind(Yn - epsilon.down, Xt)
    Xt1 <- cbind(Yn + epsilon.up, Xt)
    data <- rbind(Xt0, Xt1)
    labels <- rep(c(0, 1), each = N)
    svc <- LiblineaR(data, labels, type = 2, cost = C, bias = TRUE)
    model <- list(W = svc$W, prototypes = prototypes, params = pars, 
        kernel = kernel)
    class(model) <- "LinearizedSVR"
    return(model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
