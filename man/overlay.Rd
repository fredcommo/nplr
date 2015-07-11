\name{overlay}
\alias{overlay}
\title{
  Plotting Multiple \code{nplr} Objects
}
\description{
  To superimpose multiple logistic models fitted using \code{\link{nplr}}.
}
\usage{
  overlay(modelList = NULL, ...)
}

\arguments{
  \item{modelList}{ : A list of objects of class \code{\link{nplr}}}
  \item{...}{ : other graphical parameters. See \code{\link{par}}.}
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

# Simulating responses
set.seed(123)
drug <- rep(seq(-8, -1), 3)
Resp <- lapply(c(-6, -4, -2), function(xmid){
    scal <- rnorm(1, -1, .1)
    s <- rnorm(1, 1, .1)
    err <- rnorm(length(drug), 0, .15)
    nplr:::.nPL5(0, 1, xmid, scal, s, drug) + err
    }
)

# Computing models (to store in a list)
Models <- lapply(Resp, function(resp){
  nplr(10^drug, resp, silent = TRUE)
  })

# Visualizing
overlay(Models, xlab = "Conc.", ylab = "Resp.",
  main="Superimposing multiple curves", cex.main=1.5)

}
\seealso{
    \code{\link{plot.nplr}}
}
\keyword{datasets}
