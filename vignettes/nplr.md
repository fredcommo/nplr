%\VignetteEngine{knitr::knitr}
# nplr Vignette
### **AUTHORS**: _Frederic Commo (fredcommo@gmail.com) & Brian M. Bot_




## Overview

This function computes a weighted n-parameters logistic regression, n from 2 to 5.  
Typical applications are drug-response curve or progression curve fitting.

The n-parameter logistic regression used by $nplr()$ is of the form:

$$y = B + \frac{(T - B)}{[1 + 10^{b.(c - x)}]^s}$$

Where $B$ and $T$ are the bottom and top asymptotes, respectively, $b$ and $c$ are the Hill slope and the x-coordinate at the inflexion point, respectively, and $s$ is an asymetric coefficient.
This equation is sometimes refered to as a 5-parameter logistic regression, or the Richards equation.  
The $npars$ argument allows a user to run simplest models, while the default value $npars = "all"$ asks the function to test which model fits the best the data, with respect to a Goodness-of-Fit estimator. See the $nplr$ documentation for more details.

In a drug-response (or progression) curve fitting context, typical needs are to invert the function in order to estimate the x-value, e.g. an IC50, given a y-value (the 0.5 survival). To do so, the implemented $predict()$ method takes 2 arguments: the built model (an instance of the class nplr), and one (or a vector of) target(s), and returns the corresponding x-values and their 95% confidence intervals.



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
## [1] nplr_0.1  knitr_1.5
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.1 formatR_0.10   stringr_0.6.2  tools_3.0.1
```

