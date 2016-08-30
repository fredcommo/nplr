\name{convertToProp}
\alias{convertToProp}
\title{
  Function to Convert a Vector Into Proportions.
}
\description{
Convert a vector of values to proportions, given a minimun and a maximun value (optional).
See \code{Details} and \code{Examples}.
}
\usage{
convertToProp(y, T0 = NULL, Ctrl = NULL)
}
\arguments{
  \item{y}{ : a vector of values (responses to x).}
  \item{T0}{ : the minimal value to consider. If NULL (default), \code{min(y, na.rm=TRUE)} will be used. See \code{Details} and \code{Warning}.}
  \item{Ctrl}{ : the maximal value to consider. If NULL (default), \code{max(y, na.rm=TRUE)} will be used. See \code{Details} and \code{Warning}.}
}
\details{
In typical cell viability experiments, responses to drug concentrations (inhibition rate) may be estimated with respect to a \code{time zero (T0)} and an untreated condition values \code{(Ctrl)}, as described in [1]:

\deqn{}{ prop(y) = (y - T0)/(Ctrl - T0) }

If none of the T0 and Ctrl values are provided, \code{min(y, na.rm=TRUE)} and \code{max(y, na.rm=TRUE)} will be used, respectively. See \code{Warning}.
}
\value{
  a vector of values.
}
\references{
1 - https://dtp.nci.nih.gov/branches/btb/ivclsp.html
}
\author{
Frederic Commo, Brian M. Bot
}
\section{Warning}{
Note that, for drug response analyses, rescaling the responses between 0 to 1 using to the min and max of y, would lead to estimate a EC50 (the half effect between the maximum and the minimum of the observed effects), rather than a IC50.
}
\seealso{
  \code{\link{nplr}}
}
\note{
The data used as examples come from the NCI-60 Growth Inhibition Data:
\url{https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data},
except for multicell.tsv which are simulated data.
}
\examples{
## Using the MDA-N data
  op <- par(no.readonly=TRUE)         # save default parameters

  require(nplr)
  path <- system.file("extdata", "mdan.txt", package = "nplr")
  mdan <- read.delim(path)
  
# fit a model on the original responses (proportions of control):
  conc <- mdan$CONC
  y0 <- mdan$GIPROP
  model0 <- nplr(conc, y0)

# Adjust the data between 0 to 1, then fit a new model:
  y1 <- convertToProp(y0)
  model1 <- nplr(conc, y1)
  
  par(mfrow=c(1, 2))
  plot(model0, ylim = range(0, 1), main = "Original y values")
  plot(model1, ylim = range(0, 1), main = "Rescaled y values")
  par(op)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
