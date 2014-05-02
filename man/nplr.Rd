\name{nplr}
\alias{nplr}
\alias{nplr-class}

\title{
  Function to compute n-parameters logistic regression.
}

\description{
  This function computes a weighted n-parameters logistic regression, given x and y values. See \code{Details}.
}

\usage{
nplr(x, y, useLog = TRUE, LPweight = 0.25, npars = "all", method = c("res", "sdw", "gw"), silent = FALSE)
}

\arguments{
  \item{x}{ : a vector of numeric values, e.g. a vector of drug concentrations.}
  \item{y}{ : a vector of numeric values, e.g. a vector of responses, provided as proportions of control.}
  \item{useLog}{ : Logical. Should x-values be transformed to \eqn{log10(x)}. Default is \code{TRUE}.}
  \item{LPweight}{ : a coefficient to adjust the weights. \eqn{LPweight = 0} will compute a non-weighted np-logistic regression.}
  \item{npars}{ : a numeric value (or \code{"all"}) to specify the number of parameters to use in the model. If \code{"all"} the logistic model will be tested with 2 to 5 parameters, and the best option will be returned. See \code{Details}}
  \item{method}{ : a character string to specify what weighted method to use. Options are \code{"res"(Default), "sdw", "gw"}. See \code{Details}}
  \item{silent}{ : Logical. Specify whether \code{warnings} ad/or \code{messages} has to be silenced. Default is \code{FALSE}.}
}

\details{
The n-parameter logistic regression is of the form:
\deqn{
y = B + (T - B)/[1 + 10^(b*(xmid - x))]^s
}
Where \code{B} and \code{T} are the bottom and top asymptotes, respectively, \code{b} and \code{xmid} are the Hill slope and the x-coordinate at the inflexion point, respectively, and s is an asymetric coefficient.
This equation is sometimes refered to as a 5-parameter logistic regression, or the Richards equation.

When specifying \code{npars = 4}, the \code{s} parameter is forced to be \code{1}, and the corresponding model is a 4-parameter logistic regression, symetrical around its inflexion point.
when specifying \code{npars = 3 or 2}, add 2 more constraints and force \code{B} and \code{T} to be \code{0} and \code{1}, respectively.

Weighted methods:

The model parameters are optimized, simultaneously, using \href{http://stat.ethz.ch/R-manual/R-devel/library/stats/html/nlm.html}{nlm}, given a sum of squared errors function, sse(Y), to minimize:

\deqn{
  sse(Y) = \Sigma [W.(Yobs - Yfit)^2 ]
}

where Yobs, Yfit and W are the vectors of observed values, fitted values and weights, respectively.

In order to reduce the effect of possible outliers, the weights passed to sse(Y) can be computed in different ways, specified in \code{nplr}.

The residuals-weights (method \code{"res"}) computes the weights as:
\deqn{
  W = (1/residuals)^LPweight
}
where \code{residuals} and \code{LPweight} are the quadratic difference between the observed and fitted values, and a weight-adjustment coefficient, respectively. Best results are generally obtained by setting \eqn{LPweight = 0.25} (default value), while setting \eqn{LPweight = 0} results in computing a non-weighted sum of squared errors.

The standard-weights method (\code{"sdw"}) computes the weights as:
\deqn{
  W = 1/Var(Yobs_r)
}
where \code{Var(Yobs_r)} is the vector of the within-replicates variances.

The general-weights method (\code{"gw"}) computes weights as:
\deqn{
  W = 1/yfit^LPweight
}
where \code{yfit} and \code{LPweight} are the fitted values and a weight-adjustment coefficient, respectively. As for the residuals-weights method, setting \eqn{LPweight = 0} results in computing a non-weighted sum of squared errors.
}

\value{
An object of class \code{nplr}.
}

\section{slots}{
  \itemize{
    \item{x}{ : the x values as they are used in the model. It can be \code{Log10(x)} if \code{useLog} was set to \code{TRUE}.}
    \item{y}{ : the y values.}
    \item{useLog}{ : logical.}
    \item{npars}{ : the best number of parameters if \code{npars="all", the specified number of parameters, otherwise.}}
    \item{LPweight}{ : the weights tuning parameter.}
    \item{yFit}{ : the y fitted values.}
    \item{xCurve}{ : the x values generated to draw the curve. 200 points between the \code{min} and \code{max} of x.}
    \item{yCurve}{ : the fitted values used to draw the curve. the fitted values corresponding to \code{xCurve}.}
    \item{inflPoint}{ : the inflexion point x and y coordinates.}
    \item{goodness}{ : the goodness-of-fit. The correlation between the fitted and the observed y values}
    \item{stdErr}{ : the mean squared error between the fitted and the observed y values}
    \item{pars}{ : the model parameters.}
    \item{AUC}{ : the area under the curve estimated using both the trapezoid method and the Simpson's rule.}
    }
}

\references{
Richards, F. J. (1959). A flexible growth function for empirical use. J Exp Bot 10, 290-300.

Giraldo J, Vivas NM, Vila E, Badia A. Assessing the (a)symmetry of concentration-effect curves: empirical versus mechanistic models. Pharmacol Ther. 2002 Jul;95(1):21-45.

Motulsky HJ, Brown RE. Detecting outliers when fitting data with nonlinear regression - a new method based on robust nonlinear regression and the false discovery rate. BMC Bioinformatics. 2006 Mar 9;7:123.

NCI-60 Growth Inhibition Data: \url{https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Growth+Inhibition+Data}
}
\author{
Frederic Commo, Brian M. Bot
}
\seealso{
    \code{\link{convertToProp}}, \code{\link{getEstimates}}, \code{\link{plot.nplr}}, \code{\link{nplrAccessors}}
}
\examples{
# Using the PC-3 data
  require(nplr)
  path <- system.file("extdata", "pc3.txt", package="nplr")
  pc3 <- read.table(path, header=TRUE)
  model <- nplr(x=pc3$CONC, y=pc3$GIPROP)
  plot(model)
}
