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

#' Sample spline parameters
#' @param model spline GAM
#' @param samples number of samples
#' @keywords internal
sample.params = function(model, samples) {
  mu = stats::coef(model)
  sig = model$Vp
  sigL = mgcv::mroot(sig)
  colL = ncol(sigL)

  t(mu + sigL %*% matrix(rnorm(colL * samples), colL, samples))
}

#' @keywords internal
quantile.multi = function(data, probs, prob.names) {
  names(data) %>%
    ldply(function(name) {
      data.frame(t(quantile(data[[name]], probs=probs, names=FALSE, na.rm=TRUE))) %>%
        setNames(prob.names) %>%
        mutate(name=name)
    }) %>%
    melt(id.vars="name") %>%
    mutate(variable=sprintf(as.character(variable), name)) %>%
    dcast(. ~ variable) %>%
    select(-.)
}

#' Calculate quantiles of model outcome variables
#' @param data data.frame containing predictor variables and outcome variables
#' @param model associated model (its formula is used to determine which variables in \code{data} are predictors and which are outcomes)
#' @param probs probabilities at which quantiles will be calculated (same as for `quantile`)
#' @param prob.names format strings used to construct names of the new data columns. For example, if \code{probs} is 0.5 and
#' \code{prob.names} is \code{\%s.median}, then the median of each outcome variable \code{var} will be returned in \code{var.median}.
#' @return data frame of predictors and associated outcome quantiles.
#' @keywords internal
quantile.outcomes = function(data, model, probs, prob.names) {
  predictors = intersect(names(data), all.vars(model$pred.formula))

  confints = function(samples) {
    outcomes = samples %>% select(-pspline.sample)
    if (!is.null(predictors)) {
      common = outcomes %>% select(predictors) %>% first()
      outcomes %<>% select(-predictors)
    }
    outcomes %<>% quantile.multi(probs, prob.names)

    if (!is.null(predictors)) {
      outcomes %<>% cbind(common)
    }

    return(outcomes)
  }

  if (is.null(predictors)) {
    data %>% do(confints(.))
  } else {
    data %>%
      group_by_at(predictors) %>%
      do(confints(.)) %>%
      ungroup()
  }
}

#' Calculate confidence intervals of model outcome variables
#' @param data data frame containing predictor and outcome variables
#' @param model associated model
#' @param level confidence level
#' @return data frame of predictors and associated lower CL, median, and upper CL
#' @keywords internal
confints.outcomes = function(data, model, level) {
  quantile.outcomes(
    data, model,
    c((1 - level) / 2, .5, (1 + level) / 2),
    c("%s.lower", "%s.median", "%s.upper")
  )
}
