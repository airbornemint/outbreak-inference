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

#' Calculate case count for an outbreak
#'
#' This is useful as \code{outcome} for pspline.estimate.scalars.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param params model parameter matrix
#' @param time vector of time values at which the model will be evaluated
#' @return time series of cumulative case counts
#' @export
pspline.outbreak.cases = function(model, params, predictors) {
  # Get model predictions for given param values
  fit = predict(model, predictors, type="lpmatrix") %*% params

  # Map spline fit back to data
  data.frame(cases=model$family$linkinv(fit)) %>% cbind(predictors)
}

#' Calculate cumulative case count for an outbreak
#'
#' This is useful as \code{outcome} for pspline.estimate.timeseries.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param params model parameter matrix
#' @param time vector of time values at which the model will be evaluated
#' @return time series of relative cumulative cases
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

#' @export
pspline.outbreak.cumcases.relative = function(model, params, predictors) {
  model %>%
    pspline.outbreak.cumcases(params, predictors) %>%
    mutate(cumcases.relative=cumcases / max(cumcases)) %>%
    select(-cumcases)
}

# Calculate outbreak thresholds for an outbreak
#' This is useful as \code{outcomes} for pspline.estimate.scalars.
#'
#' @param onset onset threshold as fraction of total outbreak case count
#' @param offset offset threshold as fraction of total outbreak case count
#' @return data frame with columns \code{onset} and \code{offset} representing time
#' when the outbreak crossed onsed and offset thresholds
#' @export
pspline.outbreak.thresholds = function(onset=NA, offset=NA) {
  function(model, params, predictors) {
    # Calculate cumulative case counts from the model and parameters
    cumfrac = pspline.outbreak.cumcases.relative(model, params, predictors)

    data.frame(
      onset=threshold.ts(predictors$time, cumfrac, onset),
      offset=threshold.ts(predictors$time, cumfrac, offset)
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

#' @export
pspline.outbreak.calc.cumcases = function(time, cases) {
  cumsum(c(1, diff(time)) * cases)
}
