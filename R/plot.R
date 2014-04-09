setMethod(
  f = "plot", 
  signature = "nplr",
  definition = function(object, x=NA, y=NA, pcol="aquamarine1", lcol="red3", cex=1.5,
                        showTarget=.5, showGOF=TRUE, showCI=TRUE, showInfl=FALSE, B=1e4, unit='',
                        Title=NA, xlab='Log10(Drug[c])', ylab='Survival', ...){

    x <- getX(object)
    y <- getY(object)
    newx <- getXcurve(object)
    newy <- getYcurve(object)
    my <- as.numeric(by(y, x, mean, na.rm=TRUE))
    mx <- unique(x)
    gof <- round(getGoodness(object), 3)
    plot(x, y, col=pcol, cex=cex, pch=19, xlab=xlab, ylab=ylab,
         las = 1, cex.axis = 1.25, cex.lab = 1.5,...)
    points(x, y, pch = 1, cex = cex)
    
    if(showGOF)
      legend(ifelse(newy[length(newy)]<newy[1], 'topright', 'bottomright'),
             legend = paste('Goodness of fit:', gof), bty = 'n', cex = 1.5)
    
    if(!(!showTarget)){
      stdErr <- getStdErr(object)
      estim <- .estimateRange(showTarget, stdErr, getPar(object)$params, B, object@useLog)
      legend1 <- sprintf("IC%d : %s%s", showTarget*100, format(estim[2], scientific=TRUE), unit)
      legend2 <- sprintf("[%s, %s]", format(estim[1], scientific=TRUE), format(estim[3], scientific=TRUE))
      legend(ifelse(newy[length(newy)]<newy[1], 'bottomleft', 'topleft'),
             legend = c(legend1, legend2), cex = 1.5, text.col = 'steelblue4', bty = 'n')
    }
    
    if(showCI){
      bounds <- .confInt(getStdErr(object), getY(object), getFitValues(object), newy)
      xx <- c(newx, rev(newx))
      yy <- c(bounds$lo, rev(bounds$hi))
      polygon(xx, yy, border = NA, col = rgb(.8,.8,.8,.4))
    }
    
    if(showInfl)
      points(getInflexion(object), pch=19, cex=2, col="blue")
    
    lines(newy ~ newx, col=lcol, lwd=4)
    if(object@LPweight != 0){
      Sub = sprintf("Weighted %s-P logistic regr. (nplr package, version: %s)", object@npars, packageVersion("nplr"))
    } else{ 
      Sub = sprintf("Non-weighted %s-P logistic regr. (nplr package, version: %s)", object@npars, packageVersion("nplr"))
    }
    title (main = Title, sub = Sub, cex.sub = .75)
  }
)
