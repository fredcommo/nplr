\name{getEstimates}
\alias{getEstimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to estimate x given y.
}
\description{
  This function takes as its first argument a model returned by \code{nplr}. By inverting the logistic model, it returns the x values estimates corresponding to one (or a vector of) y target(s) provided. The standard error of the model, defined as the mean squared error on the fitted values, is used to compute a 95\% confidence interval on the predicted x values.
}
\usage{
  getEstimates(object, target=seq(.9, .1, by=-.1), B=1e4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class \code{nplr}.
}
  \item{target}{
  one, of a vector of, numerical value(s) for which the corresponding x has tobe estimated. Default are target values from .9 to .1.
}
  \item{B}{
  the lenght of the y distribution from which the x 95\% confidence interval is estimated.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data set containing:
  \item{y}{
  the target value provided.
}
  \item{x05}{
  the lower bound of the estimated 95\% confidence interval.
}
  \item{x50}{
  the estimated value.
}
  \item{x95}{
  the upper bound of the estimated 95\% confidence interval.
}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Frederic Commo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{nplr}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
