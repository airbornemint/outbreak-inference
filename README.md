## Summary

This R package lets you make estimates (with confidence intervals) of outcomes using a generalized additive model (GAM) with penalized basis splines (P-splines), which is useful for non-linear systems such as sporadic and seasonal infectious disease outbreaks.

To accomplish this, you:

1. Use the `mgcv` package to obtain a generalized additive model (GAM) or generalized additive mixed model (GAMM) of your system. 
2. Write a function that calculates the desired outcome measure
3. Call the appropriate function in `pspline.inference` with the model from step 1 and the function from step 2.

With this information, `pspline.inference` samples the outcome measure of interest from the model to calculate the estimates.

This package can handle two types of estimates:

 * Scalar estimates, in which a single numerical outcome is calculated from the system (for example, the time when an infectious disease outbreak reaches its peak)
 * Time series estimates, in which a time series of outcomes is calculated from the system (for example, the cumulative incidence during an outbreak).
 
 See package documentation in R — `?pspline.inference` — for more information, as well as the package vignette — `vignette('seasonal')`.
 
## Installation
 
```
install.packages("devtools")
library(devtools)
install_gitlab("airbornemint/pspline.inference@main")
?pspline.inference
```
