## DEFINE nplr CLASS
setClass(
  Class='nplr', 
  representation(
    x='numeric', 
    y='numeric', 
    useLog='logical',
    npars='numeric', 
    LPweight='numeric',
    yFit='numeric', 
    xCurve='numeric', 
    yCurve='numeric',
    inflPoint='data.frame', 
    goodness='numeric', 
    stdErr='numeric',
    pars='data.frame',
    AUC='data.frame'),
  
  prototype = prototype(
    useLog = TRUE,
    npars = 0,
    LPweight = 0,
    yFit = numeric(),
    xCurve = numeric(), 
    yCurve = numeric(), 
    inflPoint = data.frame(), 
    goodness = 0, 
    stdErr = 0, 
    pars = data.frame(),
    AUC = data.frame()
    )
)

## SHOW METHOD FOR THIS CLASS
setMethod(
  f = 'show', 
  signature = 'nplr',
  definition = function(object){
    cat("Instance of class nplr\n")
    cat("\n")
    cat(sprintf("%s-P logistic model\n", object@npars))
    cat("Bottom asymptote:", getPar(object)$params$bottom, "\n")
    cat("Top asymptote:", getPar(object)$params$top, "\n")
    cat("Inflexion point at (x, y):", as.numeric(getInflexion(object)), "\n")
    cat("Goodness of fit:", getGoodness(object), "\n")
    cat("Standard error:", getStdErr(object), "\n")
    cat("\n")
  }
)
