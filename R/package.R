# Copyright 2017-2019 Ben Artin
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Inference using penalized basis splines (P-splines) in a generalized additive model (GAM), with applications in
#' infectious disease outbreak modeling
#'
#' This package lets you make point and interval estimates of outcomes modeled with a non-linear P-spline GAM.
#'
#' Applications in infectious disease outbreak modeling include estimating of outbreak onset, peak, or offset,
#' as well as outbreak cumulative incidence over time.
#'
#' The package can model two types of outcomes: scalar outcomes, which are single-value
#' outcome measures (for example, timing of outbreak peak) and time series characteristics, which are functions of
#' time (for example, infection incidence over time)
#'
#' For each outcome measure, the package produces median and confidence interval estimates.
#'
#' Typical use of this package begins by using the package \code{\link{mgcv}} to obtain
#' a GAM/GAMM model of the process under investigation (such as an infectious disease outbreak), followed by
#' calling either \code{\link{pspline.estimate.scalars}} or \code{\link{pspline.estimate.timeseries}} to
#' obtain confidence intervals on the desired outcome measure
#'
#' Both \code{\link{pspline.estimate.scalars}} and \code{\link{pspline.estimate.timeseries}}
#' allow computation of arbitrary outcome measures, by passing a function that calculates
#' the desired outcome measure into \code{\link{pspline.estimate.scalars}} or \code{\link{pspline.estimate.timeseries}}.
#'
#' For convenience, this package also includes several utilities specifically aimed at modeling of
#' infectious disease outbreaks, such as \code{\link{pspline.outbreak.cases}} and \code{\link{pspline.outbreak.cumcases}}
#' (for estimation of incidence and cumulative incidence), and \code{\link{pspline.outbreak.thresholds}}, for estimation
#' of outbreak onset and offset.
#'
#' @name pspline.inference
#' @docType package
#' @author Ben Artin \email{ben@@artins.org}
#'
#' @examples
#' # Simulate an outbreak for analysis
#' cases = data.frame(
#'   time=seq(0, 51),
#'   cases=rpois(52, c(rep(1, 13), seq(1, 50, length.out=13), seq(50, 1, length.out=13), rep(1, 13)))
#' )
#' 
#' # Generate GAM model for outbreak; see mgcv for details
#' library(mgcv)
#' model = gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=cases)
#' 
#' # Generate time series at which model will be evaluated for estimates
#' # Usually you want this to be the same as the time interval that your observations are in, except
#' # divided into small increments (here, eps). Using a smaller eps gives more accurate estimates, 
#' # but takes longer to run. A value smaller than 0.25 would be better for final analysis
#' eps = .1
#' estTimes = data.frame(time=seq(min(cases$time) - 0.5, max(cases$time) + 0.5 - eps, by=eps))
#' 
#' # Estimate incidence
#' estCases = pspline.estimate.timeseries(
#'   model, estTimes,
#'   pspline.outbreak.cases,
#'   # Using a large number of samples makes the analysis more robust; 
#'   # using only 40 samples makes this example run fast (default is 2000)
#'   samples=40, 
#'   level=.95
#' )
#' 
#' # Estimate time when outbreak crosses 5\% and 95\% of cumulative case count
#' onsetThreshold = 0.025
#' offsetThreshold = 1 - onsetThreshold
#' thresholds = pspline.estimate.scalars(
#'   model, estTimes,
#'   pspline.outbreak.thresholds(onset=onsetThreshold, offset=offsetThreshold), 
#'   # Using a large number of samples makes the analysis more robust; 
#'   # using only 40 samples makes this example run fast (default is 2000)
#'   samples=40, 
#'   level=.95
#' )
#' 
#' # Plot cumulative incidence estimates and threshold estimates
#' library(ggplot2)
#' ggplot() +
#'   geom_ribbon(data=estCases, aes(x=time, ymin=cases.lower, ymax=cases.upper), fill=grey(.75)) +
#'   geom_line(data=estCases, aes(x=time, y=cases.median)) +
#'   geom_point(data=cases, aes(x=time, y=cases)) +
#'   annotate("rect",
#'     xmin=thresholds$onset.lower,
#'     xmax=thresholds$onset.upper,
#'     ymin=-Inf, ymax=Inf, alpha=.25) +
#'   annotate("rect",
#'     xmin=thresholds$offset.lower,
#'     xmax=thresholds$offset.upper,
#'     ymin=-Inf, ymax=Inf, alpha=.25) +
#'  labs(x="Time", y="Incidence")
#'
#' @importFrom stats coef na.omit predict quantile rnorm ecdf
#' @importFrom utils head tail
#' @importFrom mgcv mroot
#' @importFrom dplyr bind_rows rename rename_at arrange mutate select do group_by_at ungroup first summarize_all select_at vars contains rename_all summarize across
#' @importFrom reshape2 melt dcast
#' @importFrom plyr adply ldply
#' @importFrom stats setNames
#' @importFrom magrittr %>% %<>%
#' @importFrom assertthat assert_that
#' @importFrom plotrix std.error
#' @importFrom rlang .data
NULL
