## METHODS FOR EXTRACTING INFORMATION FROM THE nplr CLASS
setMethod("getX", "nplr", function(object) return(object@x))
setMethod("getY", "nplr", function(object) return(object@y))
setMethod("getFitValues", "nplr", function(object) return(object@yFit))
setMethod("getXcurve", "nplr", function(object) return(object@xCurve))
setMethod("getYcurve", "nplr", function(object) return(object@yCurve))
setMethod("getInflexion", "nplr", function(object) return(object@inflPoint))
setMethod("getPar", "nplr", function(object){return(list(npar=object@npars, params=object@pars))})
setMethod('getGoodness', 'nplr', function(object) return(object@goodness))
setMethod('getStdErr', 'nplr', function(object) return(object@stdErr))
setMethod("getAUC", "nplr", function(object) return(object@AUC))
setMethod('getEstimates', 'nplr', function(object) return(object@estimates))


## MAIN nplr FUNCION
nplr <- function(x, y, useLog=TRUE, LPweight=0.25,
                 npars="all", method=c("res", "sdw", "gw", "Y2", "pw"),
                 B=1e4, silent=FALSE,...){
  
  if(length(x)!=length(y))
    stop("x and y lengths differ.\n")
  
  if(is.numeric(npars) & (npars<2 | npars>5))
    stop("\nThe number of parameters (npars) has to be in [2, 5], or 'all'!\n")
  
  method <- match.arg(method)

  repTable <- table(x)
  maxrep <- max(repTable, na.rm=TRUE)
  minrep <- min(repTable, na.rm=TRUE)
  if(method=="sdw"){
    if(maxrep<2){
      method <- "res"
      warning("\nNone of the x-values seem to be replicated. The 'sdw' method has been replaced by 'res'.\n",
            call.=FALSE, immediate.=TRUE)
    } else if(minrep<2){
      warning("\nOne (or more) points have no replicates. The 'sdw' method may not be appropriate.\n",
              call.=FALSE, immediate.=TRUE)
    }
  }
  
  if(method=="gw" & any(y<0))
    warning("\nBecause of one (or more) y negative values, the 'gw' method may not be appropriate.\n",
            call.=FALSE, immediate.=TRUE)

  if(any(is.na(x) | is.na(y))){
    NAs <- union(which(is.na(x)), which(is.na(y)))
    x <- x[-NAs]
    y <- y[-NAs]
    warning(call.=FALSE,
            sprintf("%s point(s) has(ve) been removed for missingness.\n", length(NAs)),
            immediate.=TRUE)
  }
  y <- y[order(x)]
  x <- sort(x)
  
  pp <- sum(y<0 | y>1)/length(y)
  if(pp > .2 & !silent){
    warningtext <- "% of your y values fall outside the range [0, 1]"
    warning(call.=FALSE, sprintf("%s%s", round(pp*100, 2), warningtext), immediate.=TRUE)
    message("\t- any results output may not be representative.")
    message("\t- be sure you are using y-values as proportions.\n")
  }
  
  if(useLog) x <- log10(x)
  object <- new("nplr", x=x, y=y, useLog=useLog, LPweight=LPweight)
  
#  weights <- rep(1, length(y))
  .weight <- .chooseWeight(method)
  
  if(npars=="all"){
    npars <- .testAll(.sce, x, y, .weight, LPweight, silent)
    if(!silent)
      cat(sprintf("%s-Parameters model seems to have better performance.\n", npars))
  }
  
  nPL <- .chooseModel(npars)
  inits <- .initPars(x, y, npars)
  best <- nlm(f=.sce, p=inits, x=x, yobs=y, .weight, LPweight, nPL)
  
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
  if(length(unique(signif(yFit, 5)))==1)
    stop("nplr failed and returned constant fitted values. Your data may not be appropriate for such model.")
  w <- .weight(x, y, yFit, LPweight)
  perf <- .getPerf(y, yFit, w)
  
  # Compute simulations to estimate the IC50 conf. interval
  pars <- cbind.data.frame(bottom=bottom, top=top, xmid=xmid, scal=scal, s=s)
  targets <- seq(.9, .1, by=(-.1))
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
