\name{plot.nplr}
\alias{plot.nplr}
\title{
  Plotting \code{nplr} objects
}
\description{
  This function allows to visualize logistic models fitted using \code{\link{nplr}}.
}
\usage{
  \method{plot}{nplr}(x, pcol = "aquamarine1", lcol = "red3",
  showEstim = FALSE, showCI = TRUE, showGOF = TRUE, showInfl = FALSE,
  B = 1e4, conf.level = .95, unit = "", ...)
}
\arguments{
  \item{x}{ : an object of class \code{\link{nplr}}}
  \item{pcol}{ : the points color.}
  \item{lcol}{ : the line color.}
  \item{showEstim}{ : logical/numeric. If a numerical value is passed (a y value to reach), the estimated x value , and interval, is displayed on the plot. Default is \code{FALSE} }
  \item{showCI}{ : logical. show the estimated confidence interval}
  \item{showGOF}{ : logical. show the estimated goodness-of-fit.}
  \item{showInfl}{ : logical. add the inflexion point on the plot.}
  \item{B}{ : the length of simulated y values. Used to estimate the confidence interval}
  \item{conf.level}{ : the confidence level. See \code{\link{getEstimates}}}
  \item{unit}{ : the unit to specify when \code{showEstim} is \code{TRUE}. Default is an empty string.}
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
\note{
The data used in the examples are samples from the NCI-60 Growth Inhibition Data: \url{https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data}
}
\examples{
# Using the PC-3 data
  require(nplr)
  path <- system.file("extdata", "pc3.txt", package = "nplr")
  pc3 <- read.delim(path)
  model <- nplr(x = pc3$CONC, y = pc3$GIPROP)
  plot(model, showEstim = 0.5, unit = "nM")
}
\seealso{
    \code{\link{nplr}}, \code{\link{convertToProp}}, \code{\link{getEstimates}}
}
\keyword{datasets}
