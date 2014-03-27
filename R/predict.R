setMethod(
  f = "predict", 
  signature = "nplm", 
  definition = function(object, targets, B=1e4){
    pars <- getPar(object)
    if(any(targets<=pars$params$bottom)){
      targets[which(targets<=pars$params$bottom)] <- min(getFitValues(object))
      warning("Warning: One (or more) of the provided values were below the estimated bottom asymptote. These values have been replaced by the minimal possible value the model can estimates",
              call.=FALSE)
    }
    if(any(targets>=pars$params$top)){
      targets[which(targets>=pars$params$top)] <- max(getFitValues(object))
      warning("Warning: One (or more) of the values were above the estimated top asymptote. These values have been replaced by the maximal possible value the model can estimates",
              call.=FALSE)
    }
    estim <- lapply(targets, function(target)
      .estimateRange(target, getStdErr(object), pars$params, B, object@useLog)
    )
    estim <- cbind.data.frame(y=targets, do.call(rbind, estim))
    colnames(estim)[-1] <- c('xmin', 'x', 'xmax')
    return(estim)
  }
)
