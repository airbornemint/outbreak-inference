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
### Time series estimation
######################################################################

#' @keywords internal
pspline.calc.timeseries = function(samples, model, predictors, outcome) {
  samples %>% adply(1, function(params) {
    params %>%
      eval.model(model, predictors, outcome)
  }, .id="pspline.sample")
}

#' Runs simulations on an outbreak GAM/GAMM for the purpose of estimating time series outbreak outcomes, and
#' returns estimated time series outcomes for each simulation.
#'
#' This is mainly used internally by \code{pspline.estimate.timeseries}, but
#' it's useful if you want to calculate summary statistics of simulation results other
#' than the ones returned by \code{pspline.estimate.timeseries}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param time time values at which the model will be evaluated during simulation
#' @param outcome time series outcome generator function; see \code{\link{pspline.estimate.timeseries}} for more info
#' @param samples number of simulations to run
#' @return matrix with one row for each simulation; each row contains the
#' time series calculated by \code{outcome} for the corresponding simulation run
#'
#' @export
#' @keywords internal
pspline.sample.timeseries = function(model, predictors, outcome, samples=100) {
  # Sample model parameters -> sample outcome
  model %>%
    sample.params(samples) %>%
    pspline.calc.timeseries(model, predictors, outcome)
}

#' Calculates confidence intervals for results of simulation performed by \code{\link{pspline.sample.timeseries}}
#'
#' @param samples data frame of samples as returned by \code{\link{pspline.sample.timeseries}}
#' @param level confidence level for calculated confidence intervals
#' @return data.frame of confidence intervals
#' @export
#' @keywords internal
pspline.confints.timeseries = function(samples, model, level=0.95) {
  samples %>% confints.outcomes(model, level)
}

#' Calculates confidence intervals for time series sampled from generalized additive (mixed) model of an outbreak
#'
#' This function performs a series of Monte Carlo simulations of a GAM/GAMM outbreak model.
#' For each simulated outbreak, it calls \code{outcome} to calculate a time series for the
#' simulated outbreak (for example, the number of cumulative cases vs time).
#' It then calculates and returns the confidence interval of the simulated time series at
#' each time point across all simulations
#'
#' The \code{outcome} function must accept (\code{model}, \code{params}, \code{time}) and return a vector
#' containing the outcome time series obtained by evaluating the model at the time points given in \code{time} and
#' using the model parameters given in \code{params}.
#'
#' A typical implementation of the \code{outcome} function would call \code{predict} on
#' \code{model} and \code{time} to obtain the linear predictor matrix, and then post-multiply
#' that matrix by \code{params}. Having thus obtained model prediction at every time point,
#' it would calculate the desired time series outcome and return it in a vector.
#'
#' For example, to calculate the time series of the first derivative of incidence,
#' you might use this function for \code{outcome}:
#'
#' \code{
#' calc_deriv = function(model, params, time) {
#'   eps = 0.001
#'   predictors = predict(model, data.frame(time=time), type="lpmatrix")
#'   fit = model$family$linkinv(predictors %*% params)
#'   predictors_eps = predict(model, data.frame(time=time + eps), type="lpmatrix")
#'   fit_eps = model$family$linkinv(predictors_eps %*% params)
#'   (fit_eps - fit) / eps
#' }
#' }
#'
#' The data frame returned by \code{pspline.estimate.timeseries} contains three columns and
#' one row for each time point in \code{time}. The columns are \code{lower}, \code{median}, and
#' \code{upper}, containing the median and the confidence interval for the computed
#' outcome time series at each time point.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param predictors data frame of predictor values at which the model will be evaluated
#' @param outcome function returning calculated outcome time series, as described above
#' @param samples number of simulations to run
#' @param level confidence level for returned estimates
#' @return data frame of estimates, as described above
#' @export
pspline.estimate.timeseries = function(model, predictors, outcome, samples=1000, level=.95) {
  model %>%
    pspline.sample.timeseries(predictors, outcome, samples) %>%
    pspline.confints.timeseries(model, level)
}

