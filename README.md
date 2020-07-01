## Summary

This R package lets you make estimates (with confidence intervals) using a generalized additive model (GAM) with penalized basis splines (P-splines). It was initiall developed for infectious disease outbreak research, but it applies to other domains as well.

This package can handle two types of estimates:

 * Scalar estimates, in which a single numerical outcome is calculated from the system (for example, the time when an infectious disease outbreak reaches its peak)
 * Time series estimates, in which a time series of outcomes is calculated from the system (for example, the cumulative incidence during an outbreak).
 
## Installation
 
```
install.packages('pspline.inference')
library(pspline.inference)
```

## How to use

See package documentation in R — `?pspline.inference` — for more information, as well as the package vignette — `vignette('seasonal', package="pspline.inference")`.
 
## Revision history

* v0.24: Documentation improvements
* v0.23: Initial CRAN release