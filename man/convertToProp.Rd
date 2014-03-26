\name{convertToProp}
\alias{convertToProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to convert a y vector to proportions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
convertToProp(y, T0 = NA, Ctrl = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{T0}{
%%     ~~Describe \code{T0} here~~
}
  \item{Ctrl}{
%%     ~~Describe \code{Ctrl} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, T0 = NA, Ctrl = NA) 
{
    if (is.na(Ctrl)) 
        Ctrl <- max(y, na.rm = TRUE)
    if (is.na(T0)) 
        return(y/Ctrl)
    else return((y - T0)/(Ctrl - T0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
