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
.wsqRes <- function(x, yobs, yfit, LPweight) {
  residuals <- (yobs - yfit)^2
  return((1/residuals)^(LPweight))
}
.sdWeight <- function(x, yobs, yfit, LPweight){
  v <- as.numeric(by(yobs, x, var, na.rm=TRUE))
  v <- ifelse(is.na(v), 1, v)
  return(1/rep(v, times=table(x)))
}
.generalWeight <- function(x, yobs, yfit, LPweight){
  return(1/yfit^LPweight)
}
.Y2 <- function(x, yobs, yfit, LPweight){
  return(1/yfit^2)
}
.poissonWeight <- function(x, yobs, yfit, LPweight){
  return(1/abs(yfit))
}
.uniWeight <- function(yobs){
  return(rep(1, length(yobs)))
}
.sce <- function(pars, x, yobs, .weight, LPweight, nPL){
  bottom <- pars[1]
  top <- pars[2]
  xmid <- pars[3]
  scal <- pars[4]
  s <- pars[5]
  yfit <- nPL(bottom, top, xmid, scal, s, x)
  residuals <- (yobs - yfit)^2
  w <- .weight(x, yobs, yfit, LPweight)
  return(sum(w*residuals))
}
.chooseWeight <- function(method){
  switch(method,
         res = {.weight <- .wsqRes},
         sdw = {.weight <- .sdWeight},
         gw = {.weight <- .generalWeight},
         Y2 = {.weight <- .Y2},
         pw = {.weight <- .poissonWeight}
  )
  return(.weight)
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
  bottom <- as.numeric(quantile(y, .025, na.rm=TRUE))
  top <- as.numeric(quantile(y, .975, na.rm=TRUE))
  z <- (y - bottom)/(top - bottom)
  z[z<=0] <- 0.05; z[z>=1] <- 0.95
  lz <- log(z/(1-z))
  scal <- coef(lm(lz~x))[2]
  return(as.numeric(scal))
}
.initPars <- function(x, y, npars){
  if(npars<4) bottom <- 0 else bottom <- quantile(y, .05,na.rm=TRUE)
  if(npars<3) top <-1 else top <- quantile(y, .95,na.rm=TRUE)
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
.getPerf <- function(y, yfit, w){
  lmtest <- summary(lm(y ~ yfit, weights=w))
  fstat <- lmtest$fstatistic
  p <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  goodness <- lmtest$adj.r.squared
  n <- sum(w!=0)
  W <- n/((n-1)*sum(w))
  stdErr <- sqrt(W*sum(w*(yfit-y)^2))
  
  return(cbind.data.frame(goodness=goodness, stdErr=stdErr, p=p))
}
.testAll <- function(.sce, x, y, .weight, LPweight, silent){
  if(!silent) cat("Testing pars\n")
  err <- sapply(1, function(p){
    test2 <- try(nlm(f=.sce, p=.initPars(x, y, 2), x=x, yobs=y, .weight, LPweight, .nPL2), silent=TRUE)
    test3 <- try(nlm(f=.sce, p=.initPars(x, y, 3), x=x, yobs=y, .weight, LPweight, .nPL3), silent=TRUE)
    test4 <- try(nlm(f=.sce, p=.initPars(x, y, 4), x=x, yobs=y, .weight, LPweight, .nPL4), silent=TRUE)
    test5 <- try(nlm(f=.sce, p=.initPars(x, y, 5), x=x, yobs=y, .weight, LPweight, .nPL5), silent=TRUE)
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
  if(is.na(Xtarget)) x05 <- x50 <- x95 <- NA
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
  }
  return(as.numeric(c(x05, x50, x95)))
}
