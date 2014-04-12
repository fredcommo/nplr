\name{getEstimates}
\alias{getEstimates}
\title{
  Function to estimate x given y.
}
\description{
  This function takes as its first argument a model returned by \code{nplr}. By inverting the logistic model, it returns the x values estimates corresponding to one (or a vector of) y target(s) provided. The standard error of the model, defined as the mean squared error on the fitted values, is used to estimate a 95\% confidence interval on the predicted x values. see \code{Details}.
}
\usage{
  getEstimates(object, target=seq(.9, .1, by=-.1), B=1e4)
}
\arguments{
  \item{object}{
  : an object of class \code{nplr}.
}
  \item{targets}{
  : one, of a vector of, numerical value(s) for which the corresponding x has to be estimated. Default are target values from .9 to .1.
}
  \item{B}{
  : the lenght of the y distribution from which the x 95\% confidence interval is estimated.
}
}
\details{
  In n-parameter logistic regressions, none of the parameters follow any particular distribution from which confidence intervals can be estimated. To overcome this issue, the standard error is used to generate a normal distribution of the target passed to the function. The q5 and q95 quantiles of that distribution are used in order to provide estimates of the x5 and x95 bounds for the corresponding x value.
  
  See also \code{Warning}.
}
\value{
A data set containing:
  \item{y}{
  : the target value provided.
}
  \item{x05}{
  : the lower bound of the estimated 95\% confidence interval.
}
  \item{x50}{
  : the estimated value.
}
  \item{x95}{
  : the upper bound of the estimated 95\% confidence interval.
}

}
\references{
}
\author{
Frederic Commo, Brian M. Bot
}
\note{
}
\section{Warning}{
  Note that, if any \eqn{target<=B} or \eqn{target>=T}, in other words outside the 2 asymptotes, the maximal (or minimal) possible value the model can estimates is returned.
}
\seealso{
  \code{\link{nplr}}
}
\examples{
# Using the PC-3 data
  require(nplr)
  path <- system.file("extdata", "pc3.txt", package="nplr")
  pc3 <- read.table(path, header=TRUE)
  model <- nplr(pc3$CONC, pc3$GIPROP)
  getEstimates(model)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
