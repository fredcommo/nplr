\name{nplmAccessors}
\docType{methods}
\alias{getAUC}
\alias{getEstimates}
\alias{getFitValues}
\alias{getGoodness}
\alias{getInflexion}
\alias{getPar}
\alias{getStdErr}
\alias{getX}
\alias{getXcurve}
\alias{getY}
\alias{getYcurve}
\alias{getYProp}

\title{ \code{\link{nplm}} accessor functions }

\description{
  Methods for extracting information from objects of class \code{\link{nplm}}. Each of the below methods are simply convenience functions which extract the corresponding slots (as the name of each method suggests) from the object of class \code{\link{nplm}}.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "nplm")}}{
      getAUC(object)
      getEstimates(object)
      getFitValues(object)
      getGoodness(object)
      getInflexion(object)
      getPar(object)
      getStdErr(object)
      getX(object)
      getXcurve(object)
      getY(object)
      getYcurve(object)
      getYProp(object)
    }
  }
}

\seealso{
  \code{\link{nplogistic}}, \code{\link{nplm}}
}
