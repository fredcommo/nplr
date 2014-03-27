## METHODS FOR EXTRACTING INFORMATION FROM THE nplm CLASS
setMethod("getX", "nplm", function(object) return(object@x))
setMethod("getY", "nplm", function(object) return(object@y))
setMethod("getFitValues", "nplm", function(object) return(object@yFit))
setMethod("getXcurve", "nplm", function(object) return(object@xCurve))
setMethod("getYcurve", "nplm", function(object) return(object@yCurve))
setMethod("getInflexion", "nplm", function(object) return(object@inflPoint))
setMethod("getPar", "nplm", function(object){return(list(npar=object@npars, params=object@pars))})
setMethod('getGoodness', 'nplm', function(object) return(object@goodness))
setMethod('getStdErr', 'nplm', function(object) return(object@stdErr))
setMethod("getAUC", "nplm", function(object) return(object@AUC))
setMethod('getEstimates', 'nplm', function(object) return(object@estimates))


## MAIN nplm FUNCION
nplogistic <- function(x, y, useLog=TRUE, LPweight=0.25,
                       npars="all", method=c("res", "sdw", "gw", "Y2", "pw"), B=1e4,...){
  
  method <- match.arg(method)
  
  if(is.numeric(npars) & (npars<2 | npars>5))
    stop("\nThe number of parameters (npars) has to be in [2, 5], or 'all'!\n")
  
  if(any(is.na(x) | is.na(y))){
    NAs <- union(which(is.na(x)), which(is.na(y)))
    x <- x[-NAs]
    y <- y[-NAs]
  }
  y <- y[order(x)]
  x <- sort(x)
  
#   pp <- sum(y<0 | y>1)/length(y)
#   if(pp > .2){
#     warning(paste(round(pp*100, 2), "% of your y values fall outside the range [0, 1] - any results output may not be representative.", sep=""))
#   }
  
  if(useLog) x <- log10(x)
  object <- new("nplm", x=x, y=y, useLog=useLog, LPweight=LPweight)
  
  weights <- rep(1, length(y))
  .sce <- .chooseSCE(method)
  
  if(npars=="all"){
    npars <- .testAll(.sce, x, y, weights, LPweight)
    cat(sprintf("%s-Parameters model seems to have better performance.\n", npars))
  }
  
  nPL <- .chooseModel(npars)
  inits <- .initPars(x, y, npars)
  best <- nlm(f=.sce, p=inits, x=x, yobs=y, Weights=weights, wcoef=LPweight, nPL)
  
  # Best estimates
  bottom <- best$estimate[1]
  top <- best$estimate[2]
  xmid<-best$estimate[3]
  scal <- best$estimate[4]
  s <- best$estimate[5]
  
  # Estimating values
  newX <- seq(min(x), max(x), length=200)
  newY <- nPL(bottom, top, xmid, scal, s, newX)
  yFit <- nPL(bottom, top, xmid, scal, s, x)
  perf <- .getPerf(y, yFit)
  
  # Compute simulations to estimate the IC50 conf. interval
  pars <- cbind.data.frame(bottom=bottom, top=top, xmid=xmid, scal=scal, s=s)
  targets <- unique(yFit)
  estimates <- lapply(targets, function(target){.estimateRange(target, perf$stdErr, pars, B, object@useLog)})
  estimates <- cbind.data.frame(Resp = targets, do.call(rbind, estimates))
  colnames(estimates) <- c('y', 'xmin', 'x', 'xmax')
  
  # Inflexion point coordinates
  infl <- .inflPoint(pars)
  
  object@npars <- npars
  object@pars <- pars
  object@yFit <- yFit
  object@xCurve <- newX
  object@yCurve <- newY
  object@inflPoint <- infl
  object@goodness <- perf$goodness
  object@stdErr <- perf$stdErr
  object@estimates <- estimates
  object@AUC <- data.frame(trapezoide = .AUC(newX, newY), Simpson = .Simpson(newX, newY))
  object@nPL <- nPL
  object@SCE <- .sce
  
  return(object)
}
