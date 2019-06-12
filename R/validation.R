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

#' Run a simulation study to validate a scalar estimator
#'
#' @param fun.truth function that generates a true state of the system. Takes no arguments, returns data frame of true values for model variables
#' @param n.truths number of different truths to generate for simulation study
#' @param fun.observations function that generates a set of observations from truth. Takes one argument (truth data frame) and returns data frame of observations
#' @param n.observations number of sets of observations to generate for each truth in the simulation study
#' @param fun.model function that returns a model to be used for estimation. Takes one argument (observations data frame) and returns the model
#' @param fun.outcomes function that calculates the outcomes of interest. Same as outcomes function in \code{\link{pspline.estimate.scalars}}.
#' @param n.samples number of samples to use for estimation. See \code{\link{pspline.estimate.scalars}}.
#' @param level confidence level to use for estimation. See \code{\link{pspline.estimate.scalars}}.
#' @return data frame specifying the fraction of true values that were contained in their estimated confidence interval
#' @export
pspline.validate.scalars <- function(fun.truth, n.truths, fun.observations, n.observations, fun.model, fun.outcomes, n.samples, level) {
  1:n.truths %>%
    ldply(function(idx) {
      fun.truth() %>%
        truth.results.scalars(fun.observations, n.observations, fun.model, fun.outcomes, n.samples, level) %>%
        mutate(idx.truth=idx)
    }, .progress="text") %>%
    select(-idx.truth, -idx.observations)%>%
    summarize_all(function(col) mean(col, na.rm=TRUE))
}

#' Assess whether scalar estimates fall within their intended confidence intervals
#' @keywords internal
check.scalars <- function(estimates, expected) {
  results = estimates
  for (outcome in names(expected)) {
    lclName = sprintf("%s.lower", outcome)
    uclName = sprintf("%s.upper", outcome)
    medianName = sprintf("%s.median", outcome)
    checkName = sprintf("%s.good", outcome)
    results %<>%
      rename_at(lclName, function(x) "pspline.lower") %>%
      rename_at(uclName, function(x) "pspline.upper") %>%
      mutate(pspline.good = (pspline.lower < expected[outcome]) & pspline.upper > expected[outcome]) %>%
      rename_at("pspline.good", function(x) checkName) %>%
      select(-medianName) %>%
      select(-pspline.lower, -pspline.upper)
  }
  results
}

#' Run simulation study on one set of observations
#' @keywords internal
observed.results.scalars <- function(observed, truth, expected, fun.model, fun.outcome, n.samples, level) {
  model = observed %>% fun.model()
  truth %<>% select_at(pred.vars(model))
  model %>%
    pspline.estimate.scalars(truth, fun.outcome, n.samples, level) %>%
    check.scalars(expected)
}

#' Run simulation study on one truth
#' @keywords internal
truth.results.scalars <- function(truth, fun.observations, n.observations, fun.model, fun.outcome, n.samples, level) {
  # This model and set of observations are just used for looking up variables in fun.outcome
  expected = truth %>% fun.observations() %>% fun.model() %>% fun.outcome(truth)
  1:n.observations %>%
    ldply(function(idx) {
      truth %>%
        fun.observations() %>%
        observed.results.scalars(truth, expected, fun.model, fun.outcome, n.samples, level) %>%
        mutate(idx.observations=idx)
    }, .parallel=TRUE)
}

