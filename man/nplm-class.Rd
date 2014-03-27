\name{nplm-class}
\Rdversion{1.1}
\docType{class}
\alias{nplm-class}
\alias{getAUC,nplm-method}
\alias{getEstimates,nplm-method}
\alias{getFitValues,nplm-method}
\alias{getGoodness,nplm-method}
\alias{getInflexion,nplm-method}
\alias{getPar,nplm-method}
\alias{getStdErr,nplm-method}
\alias{getX,nplm-method}
\alias{getXcurve,nplm-method}
\alias{getY,nplm-method}
\alias{getYcurve,nplm-method}
\alias{getYProp,nplm-method}
\alias{plot,nplm-method}
\alias{predict,nplm-method}

\title{Class \code{"nplm"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("nplm", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"vector"} ~~ }
    \item{\code{y}:}{Object of class \code{"vector"} ~~ }
    \item{\code{useLog}:}{Object of class \code{"logical"} ~~ }
    \item{\code{yProp}:}{Object of class \code{"vector"} ~~ }
    \item{\code{npars}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{LPweight}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{yFit}:}{Object of class \code{"vector"} ~~ }
    \item{\code{xCurve}:}{Object of class \code{"vector"} ~~ }
    \item{\code{yCurve}:}{Object of class \code{"vector"} ~~ }
    \item{\code{inflPoint}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{goodness}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{stdErr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pars}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{estimates}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{AUC}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{nPL}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{SCE}:}{Object of class \code{"ANY"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getAUC}{\code{signature(object = "nplm")}: ... }
    \item{getEstimates}{\code{signature(object = "nplm")}: ... }
    \item{getFitValues}{\code{signature(object = "nplm")}: ... }
    \item{getGoodness}{\code{signature(object = "nplm")}: ... }
    \item{getInflexion}{\code{signature(object = "nplm")}: ... }
    \item{getPar}{\code{signature(object = "nplm")}: ... }
    \item{getStdErr}{\code{signature(object = "nplm")}: ... }
    \item{getX}{\code{signature(object = "nplm")}: ... }
    \item{getXcurve}{\code{signature(object = "nplm")}: ... }
    \item{getY}{\code{signature(object = "nplm")}: ... }
    \item{getYcurve}{\code{signature(object = "nplm")}: ... }
    \item{getYProp}{\code{signature(object = "nplm")}: ... }
    \item{plot}{\code{signature(x = "nplm")}: ... }
    \item{predict}{\code{signature(object = "nplm")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("nplm")
}
\keyword{classes}
