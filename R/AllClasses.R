## DEFINE NPLM CLASS
setClass(
  Class='nplm', 
  representation(
    x='vector', 
    y='vector', 
    useLog='logical',
    yProp='vector', 
    npars='numeric', 
    LPweight='numeric',
    yFit='vector', 
    xCurve='vector', 
    yCurve='vector',
    inflPoint='data.frame', 
    goodness='numeric', 
    stdErr='numeric',
    pars='data.frame', 
    estimates='data.frame', 
    AUC='data.frame',
    nPL='ANY', 
    SCE='ANY'),
  
  prototype = prototype(
    useLog = TRUE,
    yProp = NA,
    npars = 0,
    LPweight = 0,
    yFit = NA,
    xCurve=NA, 
    yCurve=NA, 
    inflPoint=data.frame(), 
    goodness=0, 
    stdErr=0, 
    pars=data.frame(),
    estimates=data.frame(), 
    AUC=data.frame(), 
    nPL=NULL, 
    SCE=NULL)
)

## SHOW METHOD FOR THIS CLASS
setMethod(
  f = 'show', 
  signature = 'nplm',
  definition = function(object){
    cat("Instance of class nplm\n")
    cat("\n")
    cat(sprintf("%s-P logistic model\n", object@npars))
    cat("Goodness of fit:", getGoodness(object), "\n")
    cat("Standard error:", getStdErr(object), "\n")
    cat("\n")
    cat("Estimated values:\n")
    show(getEstimates(object))
  }
)
