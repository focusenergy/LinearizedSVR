\name{LinearizedSVRPredict}
\alias{LinearizedSVRPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LinearizedSVRPredict(model, newdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
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
function (model, newdata) 
{
    tmp <- normalize(cbind(0, newdata), model$params)
    Xn <- tmp$Xn[, -1, drop = FALSE]
    Xt <- kernelMatrix(model$kernel, Xn, model$prototypes)
    Xt <- cbind(Xt, array(1, dim(Xt)[1]))
    wx.b <- Xt \%*\% model$W[-1]
    Y.hat <- as.vector(-wx.b/model$W[1])
    Y.hat <- Y.hat * (model$params$MM[1] - model$params$mm[1]) + 
        model$params$mm[1]
    return(Y.hat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
