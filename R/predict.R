setMethod(
  f = "predict", 
  signature = "nplm", 
  definition = function(object, targets, B=1e4){
    if(any(targets<0 | targets>1))
      stop("Target values have to be between 0 and 1 (fraction of y)")
    pars <- getPar(object)
    estim <- lapply(targets, function(target)
      .estimateRange(target, getStdErr(object), pars$params, B, object@useLog)
    )
    estim <- cbind.data.frame(prop=targets, do.call(rbind, estim))
    colnames(estim)[-1] <- c('xmin', 'x', 'xmax')
    return(estim)
  }
)
