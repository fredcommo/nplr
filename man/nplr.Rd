\name{nplr}
\alias{nplr}
\alias{convertToProp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to compute n-parameters logistic regression.
}
\description{
  This function computes a weighted n-parameters logistic regression, given x and y values. See \code{Details}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nplr(x, y, useLog = TRUE, LPweight = 0.25, npars = "all", method = c("res", "sdw", "gw", "Y2", "pw"), B = 10000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of numeric values, e.g. a vector of drug concentrations.
%%     ~~Describe \code{x} here~~
}
  \item{y}{ a vector of numeric values, e.g. a vector of responses, provided as proportions of control.
%%     ~~Describe \code{y} here~~
}
  \item{useLog}{ Logical. Should x-values be transformed to \eqn{log10(x)}. Default is \code{TRUE}.
%%     ~~Describe \code{useLog} here~~
}
  \item{LPweight}{a coefficient to adjust the weights. \eqn{LPweight = 0} will compute a non-weighted np-logistic regression.
%%     ~~Describe \code{LPweight} here~~
}
  \item{npars}{ a numeric value (or \code{"all"}) to specify the number of parameters to use in the model. If \code{"all"} the logistic model will be tested with 2 to 5 parameters, and the best option will be returned. See \code{Details}
%%     ~~Describe \code{npars} here~~
}
  \item{method}{ a character string to specify what weighted method to use. Options are \code{"res"(\code{Default}), "sdw", "gw", "Y2", "pw"}. See \code{Details}
%%     ~~Describe \code{method} here~~
}
  \item{B}{a numeric value indicating the number of simulations to run, in order to estimate the x-values corresponding to responses from 0.1 to 0.9, and given the model. See the \code{\link{predict()}} method for more flexibility.
%%     ~~Describe \code{B} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
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

The model parameters are optimized, simultaneously, using \code{\href{http://stat.ethz.ch/R-manual/R-devel/library/stats/html/nlm.html}{nlm}}, given an objective function to minimize.

The following methods can be called in the function in order to reduce the effect of possible outliers.

The residuals-weights (method \code{"res"}) computes the sum of squared errors as:
\deqn{
  err(Y) = \Sigma (Yobs - Yfit)^[ 2.(1 - LPweight) ]
}
where \code{yobs}, \code{yfit}, and \code{LPweight} are the observed values, the fitted values and a weight-adjustment coefficient such that \eqn{0 \le LPweight < 1}, respectively. Best results are generally obtained by setting \eqn{LPweight = 0.25} (default value), while setting \eqn{LPweight = 0} results in computing a non-weighted sum of squared errors.

The standard-weights method (\code{"sdw"}) computes the sum of squared errors as:
\deqn{
  err(Y) = \Sigma [ (Yobs - Yfit)^2/\sigma(Yobs)^2 ]
}
where \eqn{\sigma(Yobs)} is the vector of the within-replicates standard deviations.

The general-weights method (\code{"gw"}) computes the sum of squared errors as:
\deqn{
  err(Y) = \Sigma [ (Yobs - Yfit)^2/Yobs^LPweight ]
}
where \code{yobs}, \code{yfit}, and \code{LPweight} are the observed values, the fitted values and a weight-adjustment coefficient, respectively. As for the residuals-weights method, setting \eqn{LPweight = 0} results in computing a non-weighted sum of squared errors.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{nplr}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Richards, F. J. (1959). A flexible growth function for empirical use. J Exp
Bot 10, 290–300.

Giraldo J, Vivas NM, Vila E, Badia A. Assessing the (a)symmetry of
concentration-effect curves: empirical versus mechanistic models. Pharmacol Ther.
2002 Jul;95(1):21-45.

Motulsky HJ, Brown RE. Detecting outliers when fitting data with nonlinear
regression - a new method based on robust nonlinear regression and the false
discovery rate. BMC Bioinformatics. 2006 Mar 9;7:123.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
