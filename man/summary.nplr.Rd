\name{summary.nplr}
\alias{summary.nplr}
\title{
  summaryzing \code{nplr} Objects
}
\description{
  A S3 method to visualize a model summary as a table.
}
\usage{
  \method{summary}{nplr}(object, ...)
}
\arguments{
  \item{object}{ : an object of class \code{\link{nplr}}}
  \item{...}{ : other optional parameters (not used).}
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
The data used in the examples are samples from the NCI-60 Growth Inhibition Data:\cr
\url{https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data},
except for multicell.tsv which are simulated data.
}
\examples{
# Using the PC-3 data
require(nplr)
path <- system.file("extdata", "pc3.txt", package = "nplr")
pc3 <- read.delim(path)
model <- nplr(x = pc3$CONC, y = pc3$GIPROP)
summary(model)
}
\seealso{
    \code{\link{plot.nplr}}
}
\keyword{datasets}
