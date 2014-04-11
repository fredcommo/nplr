setMethod(
  f = "getEstimates", 
  signature = "nplr", 
  definition = function(object, targets=NULL, B=1e4){
    
    if(is.null(targets)){
      targets <- seq(.9, .1, by=-.1)
    } else if(any(!is.numeric(targets)))
        stop("Target values have to be numeric")

    pars <- getPar(object)
    if(any(targets<=pars$params$bottom)){
      targets[which(targets<=pars$params$bottom)] <- min(getFitValues(object))
      warning("One (or more) of the provided values were below the estimated bottom asymptote.",
              call.=FALSE, immediate.=TRUE)
      message("These values have been replaced by the minimal possible value the model can estimate.")
    }
    if(any(targets>=pars$params$top)){
      targets[which(targets>=pars$params$top)] <- max(getFitValues(object))
      warning("One (or more) of the values were above the estimated top asymptote.",
              call.=FALSE, immediate.=TRUE)
      message("These values have been replaced by the maximal possible value the model can estimate.")
    }
    estim <- lapply(targets, function(target)
      .estimateRange(target, getStdErr(object), pars$params, B, object@useLog)
    )
    estim <- cbind.data.frame(y=targets, do.call(rbind, estim))
    colnames(estim)[-1] <- c('x05', 'x', 'x95')
    return(estim)
  }
)

.invModel <- function(pars, target){
  return(pars$xmid - 1/pars$scal*log10(((pars$top - pars$bottom)/(target - pars$bottom))^(1/pars$s)-1))
}
.estimateRange <- function(target, stdErr, pars, B, useLog){
  Xtarget = .invModel(pars, target)
  if(is.na(Xtarget)) Dmin <- D <- Dmax <- NA
  else{
    Ytmp <- target + rnorm(B, 0, stdErr)
    if(any(Ytmp<pars$bottom)) Ytmp <- Ytmp[-which(Ytmp<pars$bottom)]
    if(any(Ytmp>pars$top)) Ytmp <- Ytmp[-which(Ytmp>pars$top)]
    Q <- quantile(Ytmp, probs=c(.05, .95), na.rm=T)
    estimates <- .invModel(pars, c(Q[1], target, Q[2]))
    if(useLog) estimates <- 10^estimates
    x05 <- signif(min(estimates), 2)
    x50 <- signif(estimates[2], 2)
    x95 <- signif(max(estimates), 2)
    
#     estimates <- .invModel(pars, Ytmp)
#     if(useLog) estimates <- 10^estimates
#     Q <- quantile(estimates, c(.05, .5, .95), na.rm=TRUE)
#     x05 <- signif(Q[1], 2)
#     x <- signif(Q[2], 2)
#     x95 <- signif(Q[3], 2)
  }
  return(as.numeric(c(x05, x, x95)))
}
