# -------------------------------------------------------
## SET GENERICS
# -------------------------------------------------------
setGeneric("getX", function(object) standardGeneric("getX"))
setGeneric("getY", function(object) standardGeneric("getY"))
setGeneric("getWeights", function(object) standardGeneric("getWeights"))
setGeneric("getFitValues", function(object) standardGeneric("getFitValues"))
setGeneric("getXcurve", function(object) standardGeneric("getXcurve"))
setGeneric("getYcurve", function(object) standardGeneric("getYcurve"))
setGeneric("getPar", function(object) standardGeneric("getPar"))
setGeneric("getInflexion", function(object) standardGeneric("getInflexion"))
setGeneric("getGoodness", function(object) standardGeneric("getGoodness"))
setGeneric("getStdErr", function(object) standardGeneric("getStdErr"))
setGeneric("getAUC", function(object) standardGeneric("getAUC"))
setGeneric("getEstimates",
	function(object, targets=seq(.9, .1, by=-.1), B=1e4, conf.level=.95)
		standardGeneric("getEstimates")
	)

# -------------------------------------------------------
## METHODS FOR EXTRACTING INFORMATION FROM THE nplr CLASS
# -------------------------------------------------------
setMethod("getX", "nplr", function(object) return(object@x))
setMethod("getY", "nplr", function(object) return(object@y))
setMethod("getWeights", "nplr", function(object) return(object@w))
setMethod("getFitValues", "nplr", function(object) return(object@yFit))
setMethod("getXcurve", "nplr", function(object) return(object@xCurve))
setMethod("getYcurve", "nplr", function(object) return(object@yCurve))
setMethod("getInflexion", "nplr", function(object) return(object@inflPoint))
setMethod("getPar", "nplr", function(object){
    return(list(npar=object@npars, params=object@pars))
    })
setMethod('getGoodness', 'nplr', function(object) return(object@goodness))
setMethod('getStdErr', 'nplr', function(object) return(object@stdErr))
setMethod("getAUC", "nplr", function(object) return(object@AUC))
