%\VignetteEngine{knitr::knitr}
# nplr Vignette
### **AUTHORS**: _Frederic Commo (fredcommo@gmail.com) & Brian M. Bot_




## Overview

This function computes a weighted n-parameters logistic regression, n from 2 to 5.  
Typical applications are drug-response or progression curve fitting.

The n-parameter logistic regression used by $nplr()$ is of the form:

$$y = B + \frac{(T - B)}{[1 + 10^{b.(xmid - x)}]^s}$$

Where $B$ and $T$ are the bottom and top asymptotes, respectively, $b$ and $xmid$ are the Hill slope and the x-coordinate at the inflexion point, respectively, and $s$ is an asymetric coefficient. This equation is sometimes refered to as a 5-parameter logistic regression, or the Richards equation.  
The $npars$ argument allows a user to run simplest models, while the default value $npars = "all"$ asks the function to test which model fits the best the data, with respect to a Goodness-of-Fit estimator. See the $nplr$ documentation for more details.

In a drug-response (or progression) curve fitting context, typical needs are to invert the function in order to estimate the x-value, e.g. an IC50, given a y-value (the 0.5 survival rate). To do so, the implemented $predict()$ method takes 2 arguments: the built model (an instance of the class nplr), and one (or a vector of) target(s), and returns the corresponding x-values and their 95% confidence intervals.

The $nplr()$ function has been optimized for fitting curves on y-values passed as proportions, from 0 to 1. If data are provided as original response values, e.g. optic density measurements, the $convertToProp()$ function may be helpful. If the x-values are the original drug concentrations, the $useLog=TRUE$ argument in $nplr()$ will transform it in $Log_{10}(conc.)$. Other arguments are described in the $nplr$ documentation.

A specific $plot()$ function has been implemented in order to visualize the results.

Several self-explanatory $get$ functions give an easy access to the results stored in the $nplr()$ output.



```r
sessionInfo()
```

```
## R version 3.0.1 (2013-05-16)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] RCurl_1.95-4.1 bitops_1.0-6   nplr_0.1       knitr_1.5     
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.1 formatR_0.10   stringr_0.6.2  tools_3.0.1
```



```r
# A simple example: drug-response fitting without replicates
url <- getURL("https://raw.github.com/fredcommo/IC50new/master/demo_files/demo1.tsv")
dat <- read.delim(text = url)
np1 <- nplr(dat$x, dat$y)
```

```
## Warning: 25% of your y values fall outside the range [0, 1]
```

```
## 	- any results output may not be representative.
## 	- be sure you are using y-values as proportions.
```

```
## Testing pars
## 5-Parameters model seems to have better performance.
```

```r
plot(np1)
```

<img src="figure/simpleExample.png" title="plot of chunk simpleExample" alt="plot of chunk simpleExample" style="display: block; margin: auto;" />

```r
getAUC(np1)
```

```
##   trapezoide Simpson
## 1      1.949   1.969
```

```r
getEstimates(np1)
```

```
##     y   xmin     x  xmax
## 1 0.9 0.0053 0.011 0.018
## 2 0.8 0.0160 0.024 0.034
## 3 0.7 0.0310 0.042 0.056
## 4 0.6 0.0520 0.068 0.086
## 5 0.5 0.0810 0.100 0.130
## 6 0.4 0.1200 0.150 0.190
## 7 0.3 0.1800 0.220 0.280
## 8 0.2 0.2600 0.340 0.440
## 9 0.1 0.4200 0.580 0.980
```



```r
# Here a similar example, with replicated measurements
url <- getURL("https://raw.github.com/fredcommo/IC50new/master/demo_files/demo2.tsv")
dat <- read.delim(text = url)
np2 <- nplr(dat$x, dat$y)
```

```
## Warning: 25% of your y values fall outside the range [0, 1]
```

```
## 	- any results output may not be representative.
## 	- be sure you are using y-values as proportions.
```

```
## Testing pars
## 5-Parameters model seems to have better performance.
```

```r
plot(np2)
```

<img src="figure/duplicates.png" title="plot of chunk duplicates" alt="plot of chunk duplicates" style="display: block; margin: auto;" />

```r
getAUC(np2)
```

```
##   trapezoide Simpson
## 1      2.013   2.032
```

```r
getEstimates(np2)
```

```
##     y   xmin      x   xmax
## 1 0.9 0.0042 0.0077  0.070
## 2 0.8 0.0069 0.0290  0.099
## 3 0.7 0.0130 0.0520  0.150
## 4 0.6 0.0260 0.0820  0.230
## 5 0.5 0.0440 0.1300  0.360
## 6 0.4 0.0720 0.1900  0.730
## 7 0.3 0.1100 0.3100  1.800
## 8 0.2 0.1700 0.5700  7.200
## 9 0.1 0.2600 1.3000 32.000
```



```r
# An example of progression curve fitting where the x-values Log10
# transformation is not required and the y-values are in some arbitrary
# scale.
url <- getURL("https://raw.github.com/fredcommo/IC50new/master/demo_files/demo4.tsv")
dat <- read.delim(text = url)

# convert the y-values to proportions
x <- dat$x
yp <- convertToProp(dat$y, 5, 110)
np3 <- nplr(x, yp, useLog = FALSE)
```

```
## Testing pars
## 5-Parameters model seems to have better performance.
```

```
## Warning: production de NaN
```

```r
plot(np3, showTarget = FALSE, xlab = "Time (hrs)", ylab = "Progression")
```

<img src="figure/progressCurve.png" title="plot of chunk progressCurve" alt="plot of chunk progressCurve" style="display: block; margin: auto;" />

```r
getInflexion(np3)
```

```
##       x      y
## 1 21.83 0.4705
```



```r
# On the same data: how the number of parameters can affect the fitting
url <- getURL("https://raw.github.com/fredcommo/IC50new/master/demo_files/demo4.tsv")
dat <- read.delim(text = url)

# convert the y-values to proportions
x <- dat$x
yp <- convertToProp(dat$y, 5, 110)
models <- lapply(2:5, function(p) {
    tmp <- nplr(x, yp, useLog = FALSE, npars = p)
    list(x = getXcurve(tmp), y = getYcurve(tmp), infp = getInflexion(tmp), gof = getGoodness(tmp))
})
```

```
## Warning: production de NaN
## Warning: production de NaN
## Warning: production de NaN
```

```r

plot(x, yp, col = "grey", pch = 19, cex = 1.25)
for (i in 1:length(models)) {
    tmp <- models[[i]]
    lines(tmp$x, tmp$y, lwd = 5, col = i)
    points(tmp$infp, pch = 19, col = i, cex = 2)
    legend(0, 1 - i/10, legend = sprintf("%s-par: gof=%s", i + 1, round(tmp$gof, 
        3)), lwd = 2, col = i, bty = "n")
}
```

<img src="figure/nparsEffect.png" title="plot of chunk nparsEffect" alt="plot of chunk nparsEffect" style="display: block; margin: auto;" />

