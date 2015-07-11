\name{getEstimates}
\alias{getEstimates}
\alias{getEstimates-methods}
\alias{getEstimates,nplr-method}
\title{
  Function to Estimate x Given y.
}
\description{
  This function takes as its first argument a model returned by \code{nplr()}. By inverting the logistic model, it estimates the x values corresponding to one (or a vector of) y target(s) provided. The standard error of the model, defined as the mean squared error on the fitted values, is used to estimate a confidence interval on the predicted x values, according to the specified \code{conf.level}. see \code{Details}.
}
\usage{
  \S4method{getEstimates}{nplr}(object, targets = seq(.9, .1, by = -.1), B = 1e4, conf.level = .95)
}
\arguments{
  \item{object}{ : an object of class \code{nplr}.}
  \item{targets}{ : one, of a vector of, numerical value(s) for which the corresponding x has to be estimated. Default are target values from .9 to .1.}
  \item{B}{ : the length of the y distribution from which the x confidence interval is estimated.}
  \item{conf.level}{ : the estimated x confidence interval, bounded by (1-conf.level)/2 and 1 - (1-conf.level)/2 (by default .95, which gives x.025 and x.975).}
}
\details{
  In n-parameter logistic regressions, none of the parameters follow any particular distribution from which confidence intervals can be estimated. To overcome this issue, the standard error is used to generate a normal distribution of the target(s) passed to the function. The quantiles of that distribution are used in order to provide estimated bounds for the corresponding x value, with respect to \code{conf.level}.
See also \code{Warning}.
}
\value{
A data set containing:
  \item{y}{ : the target value.}
  \item{x.05}{ : the lower bound of the estimated 95\% confidence interval (default). If another value is passed to conf.level, x will be labelled as x.(1-conf.level)/2.}
  \item{x}{ : the estimated value.}
  \item{x.95}{ : the upper bound of the estimated 95\% confidence interval (default). If another value is passed to conf.level, x will be labelled as x.1-(1-conf.level)/2.}

}
\author{
Frederic Commo, Brian M. Bot
}
\section{Warning}{
  Notice that, if any \eqn{target<=B} or \eqn{target>=T}, in other words outside the 2 asymptotes, the maximal (or minimal) possible value the model can estimates is returned.
}
\seealso{
  \code{\link{nplr}}, \code{\link{plot.nplr}}, , \code{\link{nplrAccessors}}
}
\note{
The data used in the examples are samples from the NCI-60 Growth Inhibition Data: \url{https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data}
}

\examples{
# Using the PC-3 data
  require(nplr)
  path <- system.file("extdata", "pc3.txt", package="nplr")
  pc3 <- read.delim(path)
  model <- nplr(x = pc3$CONC, y = pc3$GIPROP)
  getEstimates(model)
  getEstimates(model, c(.3, .6), conf.level = .9)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
