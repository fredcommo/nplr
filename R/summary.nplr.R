summary.nplr <- function(object, ...){
    
    pars <- unlist(getPar(object))
    pars[-1] <- format(pars[-1], digits = 3, scientific = TRUE)
    
    gof <- format(getGoodness(object), digits = 3, scientific = TRUE)
    auc <- signif(getAUC(object), 3)
    infl <- as.numeric(getInflexion(object))
    inflpt <- cbind.data.frame(
        xInfl = format(infl[1], digits = 3),
        yInfl = format(infl[2], digits = 3)
    )
    estim <- getEstimates(object, .5)
    interv <- sprintf("[%s | %s]",
                      format(estim$x.025, digits=3, scientific = TRUE),
                      format(estim$x.975, digits=3, scientific = TRUE)
    )
    ic50 <- format(estim$x, digits = 3, scientific = TRUE)
    IC50 <- cbind(IC50 = ic50, "[95%]" = interv)
    nplrVersion <- as.character(packageVersion("nplr"))
    nplrDate <- as.character(packageDescription("nplr")["Date"])
    Rversion <- as.character(version["version.string"])
    
    out <- cbind.data.frame(t(pars),
                            GOF = gof,
                            auc,
                            inflpt,
                            "Log10(IC50)" = log10(estim$x),
                            IC50,
                            "date (Y-m-d)" = format(Sys.Date(), "%Y-%m-%d"),
                            "nplr version" = sprintf("%s (%s)",nplrVersion, nplrDate),
                            "R version" = gsub("R version ", "", Rversion)
    )

    out <- data.frame(summary = t(out))
    print(out)
}

