setMethod(
  f = "getEstimates", 
  signature = "nplr", 
  definition = function(object, targets, B){
    
    if(any(!is.numeric(targets)))
        stop("Target values have to be numeric.")

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
