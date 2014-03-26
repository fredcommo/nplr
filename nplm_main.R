nplm <- function(x, y, T0=NA, Ctrl=NA, isProp=TRUE, useLog=TRUE, LPweight=0.25,
                    npars="all", method=c("res", "sdw", "gw", "Y2", "pw"), B=1e4,...){
  
  method <- match.arg(method)
  
  if(is.numeric(npars) & (npars<2 | npars>5))
    stop("\nThe number of parameters (npars) has to be in {2, 5}, or 'all'!\n")
    
  if(any(is.na(x) | is.na(y))){
    NAs <- union(which(is.na(x)), which(is.na(y)))
    x <- x[-NAs]
    y <- y[-NAs]
  }
  y <- y[order(x)]
  x <- sort(x)
  
  if(useLog) x <- log10(x)
  object <- .nplmObj(x=x, y=y, useLog=useLog, LPweight=LPweight)
  
  if(!isProp){
    object@yProp <- .survProp(y, T0, Ctrl)
  } else {
    object@yProp <- y
  }
  
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
  targets <- seq(.1, .9, by = .1)
  estimates <- lapply(targets, function(target){.estimateRange(target, perf$stdErr, pars, B, object@useLog)})
  estimates <- cbind.data.frame(Resp = targets, do.call(rbind, estimates))
  colnames(estimates) <- c('Prop', 'xmin', 'x', 'xmax')
  
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
