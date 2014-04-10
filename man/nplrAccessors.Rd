\name{nplrAccessors}
\docType{methods}
\alias{getAUC}
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

\title{ \code{\link{nplr}} accessor functions }

\description{
  Methods for extracting information from objects of class \code{\link{nplr}}. Each of the below methods are simply convenience functions which extract the corresponding slots (as the name of each method suggests) from the object of class \code{\link{nplr}}.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "nplr")}}{
      getAUC(object)
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
  \code{\link{nplr}}
}
