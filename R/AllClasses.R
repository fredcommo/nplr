# -------------------------------------------------------
## DEFINE nplr CLASS
# -------------------------------------------------------
setClass(
  Class='nplr', 
  representation(
    weightMethod='ANY',
    x='numeric', 
    y='numeric',
    w='numeric',
    useLog='logical',
    npars='numeric', 
    LPweight='numeric',
    yFit='numeric', 
    xCurve='numeric', 
    yCurve='numeric',
    inflPoint='data.frame', 
    goodness='ANY', 
    stdErr='numeric',
#    nlmErr='vector',
    pars='data.frame',
    AUC='data.frame',
    call='ANY'),
  
  prototype = prototype(
    weightMethod = NULL,
    useLog = TRUE,
    npars = 0,
    LPweight = 0,
    yFit = numeric(),
    xCurve = numeric(), 
    yCurve = numeric(), 
    inflPoint = data.frame(), 
#    goodness = list(),
    stdErr = 0,
#    nlmErr = 0,
    pars = data.frame(),
    AUC = data.frame(),
    call = NULL
    )
)

# -------------------------------------------------------
## SHOW METHOD FOR THIS CLASS
# -------------------------------------------------------
setMethod(
  f = 'show', 
  signature = 'nplr',
  definition = function(object){
    cat("Instance of class", class(object), "\n")
    cat("\n")
    cat("Call:\n")
    print(object@call)
    weightMethod <- NULL
    goodness <- getGoodness(object)
    if(object@weightMethod == "res")
        weightMethod <- "residuals"
        else if(object@weightMethod == "sdw")
            weightMethod <- "standard weights"
            else weightMethod <- "general weights"
    message("weights method: ", weightMethod)
    cat("\n")
    cat(sprintf("%s-P logistic model\n", object@npars))
    cat("Bottom asymptote:", getPar(object)$params$bottom, "\n")
    cat("Top asymptote:", getPar(object)$params$top, "\n")
    cat("Inflexion point at (x, y):", as.numeric(getInflexion(object)), "\n")
    cat("Goodness of fit:", goodness$gof, "\n")
    cat("Weighted Goodness of fit:", goodness$wgof, "\n")
    cat("Standard error:", getStdErr(object), "\n")
#    cat("nlm error(s):", object@nlmErr, "\n")
    cat("\n")
  }
)
