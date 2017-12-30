## Summary

This R package lets you make estimates (with confidence intervals) of characteristics of infectious disease outbreaks.

To accomplish this, you:

1. Use the `mgcv` package to obtain a generalized additive model (GAM) or generalized additive mixed model (GAMM) of your outbreak. 
2. Write a function that calculates the desired characteristic of a *potential* outbreak (for example, time of outbreak peak)
3. Call the appropriate function in `outbreakpredict` with the model from step 1 and the function from step 2.

With this information, `outbreakpredict` simulates a series of potential outbreaks (sampled from the distribution provided by the model obtained in step 1), calls your function from step 2 to calculate the desired characteristic of each of those outbreaks, and returns to you the confidence interval of the calculated characteristic across all simulated outbreaks.

This package can handle two types of estimated characteristics:

 * Scalar estimates, in which a numerical descriptor is calculated from an outbreak (for example, time of outbreak peak, time of outbreak start, time of outbreak end, peak incidence, outbreak duration)
 * Time series estimates, in which a time series is calculated from the incidence time series of an outbreak (for example, cumulative incidence time series).
 
 See package documentation in R — `help(outbreakpredict)` — for more information.
 
## Installation
 
```
install.packages("devtools")
library(devtools)
install_github("airbornemint/outbreak-predict")
```