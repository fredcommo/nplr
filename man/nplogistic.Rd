\name{nplogistic}
\alias{nplogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to compute n-parameter logistic regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nplogistic(x, y, useLog = TRUE, LPweight = 0.25, npars = "all", method = c("res", "sdw", "gw", "Y2", "pw"), B = 10000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{useLog}{
%%     ~~Describe \code{useLog} here~~
}
  \item{LPweight}{
%%     ~~Describe \code{LPweight} here~~
}
  \item{npars}{
%%     ~~Describe \code{npars} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, y, useLog = TRUE, LPweight = 0.25, npars = "all", 
    method = c("res", "sdw", "gw", "Y2", "pw"), B = 10000, ...) 
{
    method <- match.arg(method)
    if (is.numeric(npars) & (npars < 2 | npars > 5)) 
        stop("\nThe number of parameters (npars) has to be in [2, 5], or 'all'!\n")
    if (any(is.na(x) | is.na(y))) {
        NAs <- union(which(is.na(x)), which(is.na(y)))
        x <- x[-NAs]
        y <- y[-NAs]
    }
    y <- y[order(x)]
    x <- sort(x)
    pp <- sum(y < 0 | y > 1)/length(y)
    if (pp > 0.2) {
        warning(paste(round(pp * 100, 2), "\% of your y values fall outside the range [0, 1] - any results output may not be representative.", 
            sep = ""))
    }
    if (useLog) 
        x <- log10(x)
    object <- new("nplm", x = x, y = y, useLog = useLog, LPweight = LPweight)
    weights <- rep(1, length(y))
    .sce <- .chooseSCE(method)
    if (npars == "all") {
        npars <- .testAll(.sce, x, y, weights, LPweight)
        cat(sprintf("\%s-Parameters model seems to have better performance.\n", 
            npars))
    }
    nPL <- .chooseModel(npars)
    inits <- .initPars(x, y, npars)
    best <- nlm(f = .sce, p = inits, x = x, yobs = y, Weights = weights, 
        wcoef = LPweight, nPL)
    bottom <- best$estimate[1]
    top <- best$estimate[2]
    xmid <- best$estimate[3]
    scal <- best$estimate[4]
    s <- best$estimate[5]
    newX <- seq(min(x), max(x), length = 200)
    newY <- nPL(bottom, top, xmid, scal, s, newX)
    yFit <- nPL(bottom, top, xmid, scal, s, x)
    perf <- .getPerf(y, yFit)
    pars <- cbind.data.frame(bottom = bottom, top = top, xmid = xmid, 
        scal = scal, s = s)
    targets <- seq(0.9, 0.1, by = -0.1)
    estimates <- lapply(targets, function(target) {
        .estimateRange(target, perf$stdErr, pars, B, object@useLog)
    })
    estimates <- cbind.data.frame(Resp = targets, do.call(rbind, 
        estimates))
    colnames(estimates) <- c("Prop", "xmin", "x", "xmax")
    infl <- .inflPoint(pars)
    object@npars <- npars
    object@pars <- pars
    object@yFit <- yFit
    object@xCurve <- newX
    object@yCurve <- newY
    object@inflPoint <- infl
    object@goodness <- perf$goodness
    object@stdErr <- perf$stdErr
    object@estimates <- estimates
    object@AUC <- data.frame(trapezoide = .AUC(newX, newY), Simpson = .Simpson(newX, 
        newY))
    object@nPL <- nPL
    object@SCE <- .sce
    return(object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
