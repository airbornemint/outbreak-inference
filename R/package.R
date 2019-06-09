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

#' Stochastic modeling of infectious disease outbreak characteristics
#'
#' This package lets you made estimates (with confidence intervals) of characteristics of infectious
#' disease outbreaks.
#'
#' It does this by modeling disease outbreaks using generalized additive (mixed) models (GAM/GAMM), then
#' using stochastic simulation based on the GAM/GAMM model to estimate outbreak characteristics.
#'
#' For example, this package can be used to model timing of outbreak onset, peak, or offset,
#' as well as outbreak cumulative incidence over time.
#'
#' The package can model two types of outbreak characteristics: scalar characteristics, which are single-value
#' measures of an outbreak (for example, timing of peak) and time series characteristics, which are functions of
#' time (for example, cumulative incidence count over time)
#'
#' For each outbreak characteristic, the package produces median and confidence interval estimates.
#'
#' Typical use of this package begins by using the package \code{\link{mgcv}} to obtain
#' a GAM/GAMM model of the outbreak, followed by calling either \code{\link{pspline.estimate.scalars}} or
#' \code{\link{pspline.estimate.timeseries}} to obtain confidence intervals on the
#' desired scalar/timeseries outcomes of the outbreak
#'
#' Both \code{\link{pspline.estimate.scalars}} and \code{\link{pspline.estimate.timeseries}}
#' allow simulation of arbitrary outbreak characteristics, by passing a function that calculates
#' the desired characteristics into \code{\link{pspline.estimate.scalars}} or \code{\link{pspline.estimate.timeseries}}.
#'
#' For convenience, this package also includes \code{\link{pspline.calc.thresholds}}, which can
#' be use in conjunction with \code{\link{pspline.estimate.scalars}} to estimate timing of outbreak onset and
#' offset, as well as \code{\link{pspline.calc.cumcases}}, which can be used in conjuction with
#' \code{\link{pspline.estimate.timeseries}} to estimate cumulative incidence vs time for the pspline.
#'
#' @name pspline.inference
#' @docType package
#' @author Ben Artin \email{ben@@artins.org}
#'
#' @examples
#' # Simulate an outbreak for analysis
#' cases = rpois(52, c(rep(1, 13), seq(1, 50, length.out=13), seq(50, 1, length.out=13), rep(1, 13)))
#' data = data.frame(cases=cases, time=seq(0, 51))
#'
#' # Generate GAM model for outbreak; see mgcv for details
#' library(mgcv)
#' model = gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=data)
#'
#' # Generate time series at which model will be evaluated for estimates
#' # For the most part, you want the time span of the outbreak data divided into
#' # small increments (here, eps)
#' eps = .05
#' modelTime = seq(min(data$time) - 0.5, max(data$time) + 0.5 - eps, eps)
#'
#' # Estimate cumulative incidence count time series
#' cumCases = pspline.estimate.timeseries(
#'   model, modelTime,
#'   pspline.calc.cumcases, level=.95
#')
#'
#' # Estimate time when outbreak crosses 5% and 95% of cumulative case count
#' thresholds = pspline.estimate.scalars(
#'   model, modelTime,
#'   pspline.calc.thresholds(onset=0.05, offset=0.95), level=.95
#' )
#'
#' # Plot cumulative incidence estimates and threshold estimates
#' library(ggplot2)
#' ggplot(cumCases) +
#'   geom_ribbon(aes(x=time, ymin=lower, ymax=upper), fill=grey(.75)) +
#'   geom_line(aes(x=time, y=median)) +
#'   annotate("rect",
#'     xmin=thresholds$onset.lower,
#'     xmax=thresholds$onset.upper,
#'     ymin=-Inf, ymax=Inf, alpha=.25) +
#'   annotate("rect",
#'     xmin=thresholds$offset.lower,
#'     xmax=thresholds$offset.upper,
#'     ymin=-Inf, ymax=Inf, alpha=.25) +
#'   annotate("segment", x=-Inf, xend=Inf, y=0.05, yend=0.05) +
#'   annotate("segment", x=-Inf, xend=Inf, y=0.95, yend=0.95) +
#'  labs(x="Time", y="Relative cumulative incidence")
#'
#' @importFrom stats coef na.omit predict quantile rnorm
#' @importFrom utils head tail
#' @importFrom mgcv mroot
#' @importFrom dplyr bind_rows rename rename_at arrange mutate select do group_by_at ungroup first
#' @importFrom reshape2 melt dcast
#' @importFrom plyr adply ldply
#' @importFrom stats setNames
#' @importFrom magrittr %>% %<>%
#' @importFrom assertthat assert_that
NULL
