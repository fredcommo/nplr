\name{nplrAccessors}
\docType{methods}
\alias{getX}
\alias{getX-methods}
\alias{getX,nplr-method}
\alias{getY}
\alias{getY-methods}
\alias{getY,nplr-method}
\alias{getXcurve}
\alias{getXcurve-methods}
\alias{getXcurve,nplr-method}
\alias{getYcurve}
\alias{getYcurve-methods}
\alias{getYcurve,nplr-method}
\alias{getFitValues}
\alias{getFitValues-methods}
\alias{getFitValues,nplr-method}
\alias{getInflexion}
\alias{getInflexion-methods}
\alias{getInflexion,nplr-method}
\alias{getPar}
\alias{getPar-methods}
\alias{getPar,nplr-method}
\alias{getAUC}
\alias{getAUC-methods}
\alias{getAUC,nplr-method}
\alias{getGoodness}
\alias{getGoodness-methods}
\alias{getGoodness,nplr-method}
\alias{getStdErr}
\alias{getStdErr-methods}
\alias{getStdErr,nplr-method}
\alias{nplrAccessors}
\alias{nplrAccessors-methods}
\alias{nplrAccessors,nplr-method}

\title{ \code{\link{nplr}} accessor functions }

\description{
  Methods for extracting information from an object of class \code{nplr}. Each of the below methods are simply convenience functions which extract the corresponding slots (as the name of each method suggests) from the object of class \code{\link{nplr}}.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "nplr")}}{
    \itemize{
      \item getX(object)
      \item getY(object)
      \item getXcurve(object)
      \item getYcurve(object)
      \item getFitValues(object)  
      \item getInflexion(object)  
      \item getPar(object)
      \item getAUC(object)
      \item getGoodness(object)  
      \item getStdErr(object)
      }
    }
  }
}

\seealso{
  \code{\link{nplr}}, \code{\link{getEstimates}}
}
