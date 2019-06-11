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

######################################################################
### Scalar estimation
######################################################################

#' @keywords internal
pspline.calc.scalars = function(samples, model, predictors, outcomes) {
  samples %>% adply(1, function(params) {
    eval.model(params, model, predictors, outcomes)
  }, .id="pspline.sample")
}

#' Runs simulations on an outbreak GAM/GAMM for the purpose of estimating
#' scalar outbreak outcomes, and returns estimated scalar outcome values for each simulation.
#'
#' This is mainly used internally by \code{pspline.estimate.scalars}, but it's
#' useful if you want to calculate summary statistics of simulation results other than
#' the ones returned by \code{pspline.estimate.scalars}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param time time values at which the model will be evaluated during simulation
#' @param outcomes scalar outcome generator function; see \code{\link{pspline.estimate.scalars}} for more info
#' @param samples number of simulations to run
#' @return data frame with one row for each simulation; each row has \code{pspline.sample} column giving
#' a unique ID for the simulation, and one column for each scalar outcome returned by \code{outcomes}
#'
#' @export
#' @keywords internal
pspline.sample.scalars = function(model, predictors, calc.outcomes, samples=100) {
  # Sample model params -> calculate outcomes
  model %>%
    sample.params(samples) %>%
    pspline.calc.scalars(model, predictors, calc.outcomes)
}

#' Calculates confidence intervals and medians for scalar samples obtained from \code{\link{pspline.sample.scalars}}
#'
#' @param samples data frame of samples as returned by \code{\link{pspline.sample.scalars}}
#' @param level confidence level for calculated confidence intervals
#' @return data.frame of confidence intervals and medians
#' @export
#' @keywords internal
pspline.confints.scalars = function(samples, model, level=.95) {
  samples %>% confints.outcomes(model, level)
}

#' Calculates confidence intervals for scalars estimated from generalized additive (mixed) model of an outbreak
#'
#' This function performs Monte Carlo sampling of a GAM/GAMM outbreak model.
#' For each sampled curve, it calls \code{outcomess} to calculate scalar outcomes
#' It then calculates and returns the confidence interval of each scalar outcome
#'
#' The \code{outcomes} function must accept (\code{model}, \code{params}, \code{predictors}) and return a one-row data frame
#' in which each column lists the value of a single scalar outcome calculated from the model
#' estimates.
#'
#' A typical implementation of the \code{outcomes} function would call \code{predict} on
#' \code{model} and \code{predictors} to obtain model variable estimates at predictor values, then
#' calculate the scalar outcomes of interest and return them in a data frame.
#'
#' For example, to calculate the time of outbreak peak, you might use this function for \code{outcomes}:
#'
#' \code{
#' calc_peak = function(model, params, time) {
#'   incidence = predict(model, data.frame(time=time), type="response")
#'   data.frame(peak=time[which.max(incidence)])
#' }
#' }
#'
#' The data frame returned by \code{pspline.estimate.scalars} contains three columns for each
#' outcome calculated by \code{outcomes}: for outcome \code{x} returned by \code{outcomes},
#' \code{pspline.estimate.scalars} returns columns \code{x.lower}, \code{x.median}, and \code{x.upper}, corresponding
#' to lower confidence limit, median, and upper confidence limit of \code{x}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param predictors data.frame of predictor values at which the model will be evaluated
#' @param outcomes function returning calculated scalar outcomes, as described above
#' @param samples number of samples of outcomes to draw
#' @param level confidence level for estimates
#' @return data frame of estimates, as described above
#' @export
pspline.estimate.scalars = function(model, predictors, outcomes, samples=100, level=.95) {
  model %>%
    pspline.sample.scalars(predictors, outcomes, samples) %>%
    pspline.confints.scalars(model, level)
}
