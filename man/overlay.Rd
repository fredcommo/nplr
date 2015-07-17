\name{overlay}
\alias{overlay}
\title{
  Plotting Multiple \code{nplr} Objects
}
\description{
  To superimpose multiple logistic models fitted using \code{\link{nplr}}.
}
\usage{
  overlay(modelList = NULL, showLegend = TRUE, Cols = NULL, ...)
}

\arguments{
  \item{modelList}{ : list. A list of objects of class \code{\link{nplr}}.}
  \item{showLegend}{ : logical. Whether the legend has to be displayed.}
  \item{Cols}{ : character. A vector of colors to use. If \code{NULL} (default), greys will be used.}
  \item{...}{ : Other graphical parameters. See \code{\link{par}}.}
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

path <- system.file("extdata", "multicell.tsv", package="nplr")
multicell <- read.delim(path)

# Computing models (to store in a list)
cellsList <- split(multicell, multicell$cell)
Models <- lapply(cellsList, function(tmp){
  nplr(tmp$conc, tmp$resp, silent = TRUE)
  })

# Visualizing
overlay(Models, xlab = expression(Log[10](Conc.)), ylab = "Resp.",
  main="Superimposing multiple curves", cex.main=1.5)

}
\seealso{
    \code{\link{plot.nplr}}
}
\keyword{datasets}
