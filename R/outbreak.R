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
### Outbreak inference
######################################################################

#' Calculate incidence for an outbreak
#'
#' This is useful as \code{outcome} for \code{\link{pspline.estimate.scalars}}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param params model parameter matrix
#' @param predictors data frame of predictor values at which the model will be evaluated
#' @return data frame of predictor values with corresponding incidence estimates in \code{$cases}
#' @export
pspline.outbreak.cases = function(model, params, predictors) {
  # Get model predictions for given param values
  fit = predict(model, predictors, type="lpmatrix") %*% params

  # Map spline fit back to data
  data.frame(cases=model$family$linkinv(fit)) %>% cbind(predictors)
}

#' Calculate cumulative incidence for an outbreak
#'
#' This is useful as \code{outcome} for \code{\link{pspline.estimate.timeseries}}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param params model parameter matrix
#' @param predictors data frame of predictor values at which the model will be evaluated
#' @return data frame of predictor values with corresponding cumulative incidence estimates in \code{$cumcases}
#' @export
pspline.outbreak.cumcases = function(model, params, predictors) {
  assert_that(length(all.vars(model$pred.formula)) == 1, msg="Cumulative incidence currently requires time as the only predictor")
  pred.time = all.vars(model$pred.formula)[1]
  model %>%
    pspline.outbreak.cases(params, predictors) %>%
    rename(pspline.time=pred.time) %>%
    arrange(pspline.time) %>%
    mutate(cumcases=pspline.outbreak.calc.cumcases(pspline.time, cases)) %>%
    select(-cases) %>%
    rename_at("pspline.time", function(.) pred.time)
}

#' Calculate relative incidence for an outbreak
#'
#' This is useful as \code{outcome} for \code{\link{pspline.estimate.timeseries}}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param params model parameter matrix
#' @param predictors data frame of predictor values at which the model will be evaluated
#' @return data frame of predictor values with corresponding relative cumulative incidence estimates in \code{$cumcases.relative}
#' @export
pspline.outbreak.cumcases.relative = function(model, params, predictors) {
  model %>%
    pspline.outbreak.cumcases(params, predictors) %>%
    mutate(cumcases.relative=cumcases / max(cumcases)) %>%
    select(-cumcases)
}

#' Calculate outbreak thresholds for an outbreak
#'
#' The result of calling this is useful as \code{outcomes} for \code{\link{pspline.estimate.scalars}}.
#'
#' @param onset onset threshold (as fraction of total outbreak case count)
#' @param offset offset threshold (as fraction of total outbreak case count)
#' @return function suitable as outcome estimator parameter of \code{\link{pspline.estimate.scalars}}
#' @export
pspline.outbreak.thresholds = function(onset=NA, offset=NA) {
  function(model, params, predictors) {
    # Calculate cumulative case counts from the model and parameters
    cum.rel = pspline.outbreak.cumcases.relative(model, params, predictors)

    data.frame(
      onset=threshold.ts(cum.rel$time, cum.rel$cumcases.relative, onset),
      offset=threshold.ts(cum.rel$time, cum.rel$cumcases.relative, offset)
    )
  }
}

#' @keywords internal
threshold.ts = function(time, cumfrac, threshold) {
  if (is.na(threshold)) {
    return(NA)
  }

    # Linear interpolation across the time interval that where cumfrac crosses threshold
  idx = sum(cumfrac < threshold)
  timeLow = time[idx]
  timeStep = time[idx + 1] - time[idx]
  cumLow = cumfrac[idx]
  cumStep = cumfrac[idx + 1] - cumfrac[idx]
  return(timeLow + timeStep / cumStep * (threshold - cumLow))
}

#' Calculate cumulative incidence time series from incidence time series
#'
#' Correctly handles accumulating over time intervals different from 1
#'
#' @param time vector of times
#' @param cases vector of corresponding incidences
#' @return vector of corresponding cumulative incidences
#' @export
pspline.outbreak.calc.cumcases = function(time, cases) {
  cumsum(c(1, diff(time)) * cases)
}
