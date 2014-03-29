## HELPER FUNCTIONS (ALL INTERNAL TO THE PACKAGE)
.nPL5 <- function(bottom, top, xmid, scal, s,  X){
  yfit <- bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
  return(yfit)
}
.nPL4 <- function(bottom, top, xmid, scal, s,  X){
  yfit <- bottom+(top-bottom)/(1+10^((xmid-X)*scal))
  return(yfit)
}
.nPL3 <- function(bottom, top, xmid, scal, s,  X){
  yfit <- (top)/(1+10^((xmid-X)*scal))
  return(yfit)
}
.nPL2 <- function(bottom, top, xmid, scal, s,  X){
  yfit <- 1/(1+10^((xmid-X)*scal))
  return(yfit)
}
.wsqRes <- function(pars, x, yobs, Weights, wcoef, nPL) {
  bottom <- pars[1]
  top <- pars[2]
  xmid <- pars[3]
  scal <- pars[4]
  s <- pars[5]
  ytheo <- nPL(bottom, top, xmid, scal, s, x)
  residuals <- (yobs - ytheo)^2
  Weights <- (1/residuals)^(wcoef)
  return(sum(Weights*(yobs - ytheo)^2))
}
.sdWeight <- function(pars, x, yobs, Weights, wcoef, nPL){
  bottom <- pars[1]
  top <- pars[2]
  xmid <- pars[3]
  scal <- pars[4]
  s <- pars[5]
  ytheo <- nPL(bottom, top, xmid, scal, s, x)
  residuals <- (yobs - ytheo)^2
  Weights <- as.numeric(by(yobs, x, sd))
  Weights <- rep(Weights, length(x)/length(unique(x)))
  return(sum(residuals/Weights^2))
}
.generalWeight <- function(pars, x, yobs, Weights, wcoef, nPL){
  bottom <- pars[1]
  top <- pars[2]
  xmid <- pars[3]
  scal <- pars[4]
  s <- pars[5]
  ytheo <- nPL(bottom, top, xmid, scal, s, x)
  residuals <- (yobs - ytheo)^2
  return(sum(residuals/(ytheo^wcoef)))
}
.Y2 <- function(pars, x, yobs, Weights, wcoef, nPL){
  bottom <- pars[1]
  top <- pars[2]
  xmid <- pars[3]
  scal <- pars[4]
  s <- pars[5]
  ytheo <- nPL(bottom, top, xmid, scal, s, x)
  residuals <- (yobs - ytheo)^2
  return(sum(residuals/ytheo^2))
}
.poissonWeight <- function(pars, x, yobs, Weights, wcoef, nPL){
  bottom <- pars[1]
  top <- pars[2]
  xmid <- pars[3]
  scal <- pars[4]
  s <- pars[5]
  ytheo <- nPL(bottom, top, xmid, scal, s, x)
  residuals <- (yobs - ytheo)^2
  return(sum(residuals/abs(ytheo)))
}
.chooseSCE <- function(method){
  switch(method,
         sdw = {.sce <- .sdWeight},
         res = {.sce <- .wsqRes},
         Y2 = {.sce <- .Y2},
         pw = {.sce <- .poissonWeight},
         gw = {.sce <- .generalWeight}
  )
  return(.sce)
}
.chooseModel <- function(npars){
  switch(as.character(npars),
         "2" = {nPL <- .nPL2},
         "3" = {nPL <- .nPL3},
         "4" = {nPL <- .nPL4},
         "5" = {nPL <- .nPL5}
  )
  return(nPL)
}
.estimScal <- function(x, y){
  bottom <- min(y, na.rm=TRUE); top <- max(y, na.rm=TRUE)
  z <- (y - bottom)/(top - bottom)
  z[z==0] <- 0.05; z[z==1] <- 0.95
  lz <- log(z/(1-z))
  scal <- coef(lm(x ~ lz))[2]
  if(scal>1) scal <- 1/scal
  return(as.numeric(scal))
}
.initPars <- function(x, y, npars){
  if(npars<4) bottom <- 0 else bottom = min(y, na.rm=TRUE)
  if(npars<3) top <-1 else top = max(y, na.rm=TRUE)
  xmid = (max(x)+min(x))/2
  scal <- .estimScal(x, y)
  return(c(bottom, top, xmid, scal, s=1))
}
.getPars <- function(model){
  bottom <- model$estimate[1]
  top <- model$estimate[2]
  xmid<-model$estimate[3]
  scal <- model$estimate[4]
  s <- model$estimate[5]
  return(cbind.data.frame(bottom=bottom, top=top, xmid=xmid, scal=scal, s=s))
}
.fit <- function(x, y, npars, nPL, .sce, LPweight){
  best <- nlm(f=.sce, p=.initPars(x, y, npars), x=x, yobs=y, Weights=rep(1, length(x)), wcoef=LPweight, nPL)
  pars <- best$estimate
  return(nPL(pars[1], pars[2], pars[3], pars[4], pars[5], unique(x)))
}
.inflPoint <- function(pars){
  x = pars$xmid + (1/pars$scal)*log10(pars$s)
  y = pars$bottom + (pars$top - pars$bottom)*(pars$s/(pars$s+1))^pars$s
  return(cbind.data.frame(x=x, y=y))
}
.getPerf <- function(y, yfit){
  lmtest <- summary(lm(y ~ yfit))
  fstat <- lmtest$fstatistic
  p <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  goodness <- lmtest$adj.r.squared
  stdErr <- sqrt(1/(length(yfit)-2)*sum((yfit-y)^2))
  return(cbind.data.frame(goodness=goodness, stdErr=stdErr, p=p))
}
.testAll <- function(.sce, x, y, weights, LPweight){
  cat("Testing pars\n")
  err <- sapply(1, function(p){
    test2 <- try(nlm(f=.sce, p=.initPars(x, y, 2), x=x, yobs=y, Weights=weights, wcoef=LPweight, .nPL2), silent=TRUE)
    test3 <- try(nlm(f=.sce, p=.initPars(x, y, 3), x=x, yobs=y, Weights=weights, wcoef=LPweight, .nPL3), silent=TRUE)
    test4 <- try(nlm(f=.sce, p=.initPars(x, y, 4), x=x, yobs=y, Weights=weights, wcoef=LPweight, .nPL4), silent=TRUE)
    test5 <- try(nlm(f=.sce, p=.initPars(x, y, 5), x=x, yobs=y, Weights=weights, wcoef=LPweight, .nPL5), silent=TRUE)
    scores <- sapply(list(test2, test3, test4, test5), function(t){
      if(class(t)!="try-error") return(t$minimum)
      else return(Inf) 
    })
    return(scores)
  })
  return(which.min(err) + 1)
}
.AUC <- function(x, y){
  auc <- lapply(2:length(x), function(i){
    da <- x[i]-x[i-1]
    db <- y[i]-y[i-1]
    y[i]*da +1/2*db*da
  })
  return(do.call(sum, auc))
}
.Simpson <- function(x, y){
  dx <- mean(diff(x, lag = 1), na.rm = TRUE)
  n <- length(y)
  if(n%%2 != 0){
    x <- x[-n]
    y <- y[-n]
    n <- length(x)
  }
  f1 <- y[1]
  fn <- y[n]
  fy <- y[2:(n-1)]*rep(c(4, 2), (n-2)/2)
  return(dx/3*(f1 + sum(fy) + fn))
}
.confInt <- function(stdErr, yobs, yfit, newy){
  n <- length(yobs)
  ybar <- mean(yobs, na.rm = TRUE)
  t <- qt(.975, n-2)
  ci <- t*stdErr*sqrt((1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
  lo <- newy - ci
  hi <- newy + ci
  return(list(lo = lo, hi = hi))
}
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
    Dmin <- signif(min(estimates), 2)
    D <- signif(estimates[2], 2)
    Dmax <- signif(max(estimates), 2)
  }
  return(as.numeric(c(Dmin, D, Dmax)))
}
