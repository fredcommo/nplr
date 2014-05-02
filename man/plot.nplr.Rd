\name{plot.nplr}
\alias{plot.nplr}
\title{
  Plotting \code{nplr} objects
}
\description{
  This function allows to visualize logistic models fitted using \code{\link{nplr}}.
}
\usage{
  \method{plot}{nplr}(object, pcol="aquamarine1", lcol="red3",
  showTarget=FALSE, showCI=TRUE, showGOF=TRUE, showInfl=FALSE, B=1e4, conf.level=.95, unit=NA, ...)
}
\arguments{
  \item{object}{ : an object of class \code{\link{nplr}}}
  \item{pcol}{ : the points color.}
  \item{lcol}{ : the line color.}
  \item{showTarget}{ : logical. show the estimated value.}
  \item{showCI}{ : logical. show the estimated confidence interval}
  \item{showGOF}{ : logical. show the estimated goodness-of-fit.}
  \item{showInfl}{ : logical. add the inflexion point on the plot.}
  \item{B}{ : the length of simulated y values. Used to estimate the confidence interval}
  \item{conf.level}{ : the confidence level. See \code{\link{getEstimates}}}
  \item{unit}{ : the unit to specify when \code{showEstim} is \code{TRUE}}
  \item{...}{ : other graphical parameters.}
}
\details{
None
}
\source{
None
}
\references{
None
}
\examples{
# Using the PC-3 data
  require(nplr)
  path <- system.file("extdata", "pc3.txt", package="nplr")
  pc3 <- read.table(path, header=TRUE)
  model <- nplr(x=pc3$CONC, y=pc3$GIPROP)
  plot(model)
}
\seealso{
    \code{\link{nplr}}, \code{\link{convertToProp}}, \code{\link{getEstimates}}
}
\keyword{datasets}
