\name{convertToProp}
\alias{convertToProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to convert a y vector to proportions
}
\description{
Convert a y-vector of values to proportions, given a minimun and a maximun value (optional).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
convertToProp(y, T0 = NULL, Ctrl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ a vector of values (responses to x).
%%     ~~Describe \code{y} here~~
}
  \item{T0}{ the minimal value to consider. It should be the \code{time zero} measurement value in case of cell viability. If NULL (default), \code{min(y)} wil be used. See \code{Details}.
%%     ~~Describe \code{T0} here~~
}
  \item{Ctrl}{ the maximal value to consider. It should be the untreated condition measurement value in case of cell viability. If NULL (default), \code{max(y)} wil be used. See \code{Details}.
%%     ~~Describe \code{Ctrl} here~~
}
}
\details{
In typical cell viability experiments, responses to drug concentrations should be estimated with respect to a \code{time zero (T0)} value and a untreated condition value \code{(Ctrl)}, as described in [1].

\deqn{}{ prop(y) = (y - T0)/(Ctrl - T0) }

%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a vector of values.
}
\references{
1 - \url{http://dtp.nci.nih.gov/branches/btb/ivclsp.html}
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
