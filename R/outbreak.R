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

#' Calculate cumulative incidence for an outbreak
#'
#' This is useful as \code{outcome} for \code{\link{pspline.estimate.timeseries}}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param data data frame of predictor values at which the model will be evaluated
#' @return data frame of predictor values with corresponding cumulative incidence estimates in \code{$cumcases}
#' @export
pspline.outbreak.cases = function(model, data) {
  data
}

#' Calculate cumulative incidence for an outbreak
#'
#' This is useful as \code{outcome} for \code{\link{pspline.estimate.timeseries}}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param data data frame of predictor values at which the model will be evaluated
#' @return data frame of predictor values with corresponding cumulative incidence estimates in \code{$cumcases}
#' @export
pspline.outbreak.cumcases = function(model, data) {
  time.name = pred.var(model)
  cases.name = model.var(model)
  cumcases.name = sprintf("%s.cum", cases.name)

  data %>%
    rename_at(time.name, function(x) "pspline.time") %>%
    rename_at(cases.name, function(x) "pspline.cases") %>%
    arrange(pspline.time) %>%
    mutate(pspline.cumcases=pspline.outbreak.calc.cumcases(pspline.time, pspline.cases)) %>%
    select(-pspline.cases) %>%
    rename_at("pspline.cumcases", function(x) cumcases.name) %>%
    rename_at("pspline.time", function(x) time.name)
}

#' Calculate relative incidence for an outbreak
#'
#' This is useful as \code{outcome} for \code{\link{pspline.estimate.timeseries}}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param data data frame of predictor values at which the model will be evaluated
#' @return data frame of predictor values with corresponding relative cumulative incidence estimates in \code{$cumcases.relative}
#' @export
pspline.outbreak.cumcases.relative = function(model, data) {
  time.name = pred.var(model)
  cases.name = model.var(model)
  cumcases.name = sprintf("%s.cumrel", cases.name)

  data %>%
    rename_at(time.name, function(x) "pspline.time") %>%
    rename_at(cases.name, function(x) "pspline.cases") %>%
    arrange(pspline.time) %>%
    mutate(pspline.cumcases=pspline.outbreak.calc.cumcases(pspline.time, pspline.cases)) %>%
    mutate(pspline.cumcases=pspline.cumcases / max(pspline.cumcases)) %>%
    select(-pspline.cases) %>%
    rename_at("pspline.cumcases", function(x) cumcases.name) %>%
    rename_at("pspline.time", function(x) time.name)
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
  function(model, data) {
    # Calculate cumulative case counts from the model and parameters
    data.cumrel = pspline.outbreak.cumcases.relative(model, data)
    time.name = pred.var(model)
    cases.name = model.var(model)
    cumcases.name = sprintf("%s.cumrel", cases.name)

    data.frame(
      onset=threshold.ts(data.cumrel[[time.name]], data.cumrel[[cumcases.name]], onset),
      offset=threshold.ts(data.cumrel[[time.name]], data.cumrel[[cumcases.name]], offset)
    )
  }
}

#' @keywords internal
threshold.ts = function(time, cumfrac, threshold) {
  if (is.na(threshold) || is.nan(min(cumfrac)) || threshold < min(cumfrac)) {
    return(NA)
  }

  # Linear interpolation across the time interval that where cumfrac crosses threshold
  idx = sum(cumfrac <= threshold)
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
