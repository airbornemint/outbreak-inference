# Copyright 2017 Ben Artin
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
#' using stochastic simulation based on the GAM/GAMM model to predict outbreak characteristics.
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
#' a GAM/GAMM model of the outbreak, followed by calling either \code{\link{outbreak.predict.scalars}} or
#' \code{\link{outbreak.predict.timeseries}} to obtain confidence intervals on the
#' desired scalar/timeseries parameters of the outbreak
#'
#' Both \code{\link{outbreak.predict.scalars}} and \code{\link{outbreak.predict.timeseries}}
#' allow simulation of arbitrary outbreak characteristics, by passing a function that calculates
#' the desired characteristics into \code{\link{outbreak.predict.scalars}} or \code{\link{outbreak.predict.timeseries}}.
#'
#' For convenience, this package also includes \code{\link{outbreak.calc.thresholds}}, which can
#' be use in conjunction with \code{\link{outbreak.predict.scalars}} to predict timing of outbreak onset and 
#' offset, as well as \code{\link{outbreak.calc.cum}}, which can be used in conjuction with
#' \code{\link{outbreak.predict.timeseries}} to #' predict cumulative incidence vs time for the outbreak.
#'
#' @name outbreakpredict-package
#' @docType package
#' @author Ben Artin \email{ben@@artins.org}
#'
#' @examples
#' data = ... # Import data as data frame of newcases vs time
#'
#' # Generate GAM model for outbreak; see mgcv for details
#' model = gam(newcases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=data)
#'
#' # Generate time series at which model will be evaluated for predictions
#' eps = .05
#' modelTime = data.frame(time=seq(min(data$time) - 1 + eps, max(data$time), eps))
#'
#' # Predict cumulative case count time series
#' cumCases = outbreak.predict.timeseries(model, modelTime, function(model, params, newdata) {
#'   outbreak.calc.cum(model, params, newdata$time, timedelta=1)
#' })
#'
#' # Predict time when outbreak crosses 5% and 95% of cumulative case count
#' thresholds = outbreak.predict.scalars(model, modelTime, function(model, params, newdata) {
#'   outbreak.calc.thresholds(model, params, newdata$time, onset=0.05, offset=0.95)
#' })
#' 
#' @importFrom stats coef na.omit predict quantile rnorm
#' @importFrom utils head tail
#' @importFrom mgcv mroot
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
NULL

#' Sample from a multivariate normal distribution
#'
#' @param mu matrix of means
#' @param sig matrix of covariances
#' @param nsim number of samples to draw
#' @return matrix of sampled values
#' @keywords internal
randomMVN = function(mu, sig, nsim) {
  L = mroot(sig)
  m = ncol(L)
  t(mu + L %*% matrix(rnorm(m * nsim), m, nsim))
}

#' Runs simulations on an outbreak GAM/GAMM for the purpose of predicting
#' scalar outbreak parameters, and returns predicted scalar parameter values for each simulation.
#'
#' This is mainly used internally by \code{outbreak.predict.scalars}, but it's
#' useful if you want to calculate summary statistics of simulation results other than
#' the ones returned by \code{outbreak.predict.scalars}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param newdata time values at which the model will be evaluated during simulation
#' @param quant scalar parameter generator function; see \code{\link{outbreak.predict.scalars}} for more info
#' @param nsim number of simulations to run
#' @return data frame with one row for each simulation; each row has \code{outbreak.sim} column giving
#' a unique ID for the simulation, and one column for each scalar parameter returned by \code{quant}
#'
#' @export
#' @keywords internal
outbreak.predict.scalars.sim = function(model, newdata, quant, nsim=100) {
  # Multivariate normal random generator
  # Generate random model parameters
  randomParams = randomMVN(coef(model), model$Vp, nsim)

  predictions = bind_rows(apply(randomParams, 1, function(params) { quant(model, params, newdata) }))
  if (any(is.na(predictions))) {
    warning("Some predictions are NA")
    return()
  }

  predictions %>% cbind(data.frame(outbreak.sim=seq(1, nsim)))
}

#' Calculates confidence intervals for results of simulation performed by \code{\link{outbreak.predict.scalars.sim}}
#'
#' @param predictions data frame of predictions as returned by \code{\link{outbreak.predict.scalars.sim}}
#' @param level confidence level for calculated confidence intervals
#' @return data.frame of confidence intervals, as returned by \code{\link{outbreak.predict.scalars}}
#' @keywords internal
outbreak.predict.scalars.confints = function(predictions, level=.95) {
  confints = matrix(NA, nrow=1, ncol=0) %>% as.data.frame()

  for (name in names(predictions)) {
    if (name == "outbreak.sim") {
      next
    }
    nameConfInts = as.data.frame(matrix(NA, 1, 3))
    names(nameConfInts) = paste(name, c("lower", "median", "upper"), sep=".")
    nameConfInts[1,] = quantile(predictions[[name]], c(1-level, 0.5, level), names=FALSE)
    confints = confints %>% cbind(nameConfInts)
  }

  confints
}

#' Calculates confidence intervals for scalar predicted from generalized additive (mixed) model of an outbreak
#'
#' This function performs a series of Monte Carlo simulations of a GAM/GAMM outbreak model.
#' For each simulated outbreak, it calls \code{quant} to calculate a scalar
#' parameter for the simulated outbreak (for example, the time of outbreak peak).
#' It then calculates and returns the confidence interval of each simulated scalar
#' parameter across all simulations.
#'
#' The \code{quant} function must accept (\code{model}, \code{params}, \code{newdata}) and return a one-row data frame
#' in which each column lists the value of a single scalar parameter calculated by the model
#' evaluated at the time points given in \code{newdata} and using the model parameters
#' given in \code{params}.
#'
#' A typical implementation of the \code{quant} function would call \code{predict} on
#' \code{model} and \code{newdata} to obtain the linear predictor matrix, and then post-multiply
#' that matrix by \code{params}. Having thus obtained model prediction at every time point,
#' it would calculate the desired scalar parameters and return them in a data frame.
#'
#' For example, to calculate the time of outbreak peak, you might use this function for \code{quant}:
#'
#' \code{
#' calc_peak(model, params, newdata) {
#'   predictors = predict(model, data.frame(time=newdata), type="lpmatrix")
#'   fit = model$family$linkinv(predictors %*% params)
#'   data.frame(peak=newdata[which.max(fit)])
#' }
#' }
#'
#' The data frame returned by \code{outbreak.predict.scalars} contains three columns for each
#' parameter calculated by \code{quant}: for parameter \code{x} returned by \code{quant},
#' \code{outbreak.predict.scalars} returns columns \code{x.lower}, \code{x.median}, and \code{x.upper}, corresponding
#' to lower confidence limit, median, and upper confidence limit.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param newdata vector of time values at which the model will be evaluated
#' @param quant function returning calculated scalar parameters, as described above
#' @param nsim number of simulations to run
#' @param level confidence level for predictions
#' @return data frame of predictions, as described above
#' @export
outbreak.predict.scalars = function(model, newdata, quant, nsim=100, level=.95) {
  model %>%
    outbreak.predict.scalars.sim(newdata, quant, nsim) %>%
    outbreak.predict.scalars.confints(level)
}

#' Runs simulations on an outbreak GAM/GAMM for the purpose of predicting time series outbreak parameters, and 
#' returns predicted time series parameter values for each simulation.
#'
#' This is mainly used internally by \code{outbreak.predict.timeseries}, but
#' it's useful if you want to calculate summary statistics of simulation results other
#' than the ones returned by \code{outbreak.predict.timeseries}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param newdata time values at which the model will be evaluated during simulation
#' @param quant time series parameter generator function; see \code{\link{outbreak.predict.timeseries}} for more info
#' @param nsim number of simulations to run
#' @return matrix with one row for each simulation; each row contains the
#' time series calculated by \code{quant} for the corresponding simulation run
#'
#' @export
#' @keywords internal
outbreak.predict.timeseries.sim = function(model, newdata, quant, nsim=100) {
  # Generate random model parameters
  randomParams = randomMVN(coef(model), model$Vp, nsim)

  predictions = apply(randomParams, 1, function(params) { quant(model, params, newdata) })
  if (any(is.na(predictions))) {
    warning("Model did not converge, try reducing the number of spline knots in GAM")
    return()
  }

  predictions
}

#' Calculates confidence intervals for results of simulation performed by \code{\link{outbreak.predict.timeseries.sim}}
#'
#' @param predictions data frame of predictions as returned by \code{\link{outbreak.predict.timeseries.sim}}
#' @param level confidence level for calculated confidence intervals
#' @return data.frame of confidence intervals, as returned by \code{\link{outbreak.predict.timeseries}}
#' @keywords internal
outbreak.predict.timeseries.confints = function(predictions, level=0.95) {
  resultNames = c("lower", "median", "upper")
  if (is.null(predictions)) {
    confints = as.data.frame(matrix(NA, 0, 3))
    names(confints) = resultNames
  } else {
    confints = apply(predictions, 1, function(prediction) {
      ci = as.data.frame(matrix(NA, 0, 3))
      names(ci) = resultNames

      ci[1,] = quantile(prediction, c(1-level, 0.5, level), names=FALSE)
      ci
    }) %>% bind_rows()
  }

  confints
}

#' Calculates confidence intervals for time series predicted from generalized additive (mixed) model of an outbreak
#'
#' This function performs a series of Monte Carlo simulations of a GAM/GAMM outbreak model.
#' For each simulated outbreak, it calls \code{quant} to calculate a time series  for the
#' simulated outbreak (for example, the number of cumulative cases vs time).
#' It then calculates and returns the confidence interval of the simulated time series at
#' each time point across all simulations
#'
#' The \code{quant} function must accept (\code{model}, \code{params}, \code{newdata}) and return a vector
#' with the time series obtaines by evaluating the model at the time points given in \code{newdata} and
#' using the model parameters given in \code{params}.
#'
#' A typical implementation of the \code{quant} function would call \code{predict} on
#' \code{model} and \code{newdata} to obtain the linear predictor matrix, and then post-multiply
#' that matrix by \code{params}. Having thus obtained model prediction at every time point,
#' it would calculate the desired time series parameters and return it in a vector.
#'
#' For example, to calculate the time series of the first derivative of incidence,
#' you might use this function for \code{quant}:
#'
#' \code{
#' calc_deriv(model, params, newdata) {
#'   eps = 0.001
#'   predictors = predict(model, data.frame(time=newdata), type="lpmatrix")
#'   fit = model$family$linkinv(predictors %*% params)
#'   predictors_eps = predict(model, data.frame(time=newdata + eps), type="lpmatrix")
#'   fit_eps = model$family$linkinv(predictors_eps %*% params)
#'   (fit_eps - fit) / eps
#' }
#' }
#'
#' The data frame returned by \code{outbreak.predict.timeseries} contains three columns and
#' one row for each time point in \code{newdata}. The columns are \code{lower}, \code{median}, and
#' \code{upper}, containing the median and the confidence interval for the computed
#' time series at each time point.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param newdata vector of time values at which the model will be evaluated
#' @param quant function returning calculated time series, as described above
#' @param nsim number of simulations to run
#' @param level confidence level for returned predictions
#' @return data frame of predictions, as described above
#' @export
outbreak.predict.timeseries = function(model, newdata, quant, nsim=1000, level=.95) {
  model %>%
    outbreak.predict.timeseries.sim(newdata, quant, nsim) %>%
    outbreak.predict.timeseries.confints(level)
}

#' Predict cumulative case count for an outbreak
#'
#' This is useful as \code{quant} for outbreak.predict.timeseries.
#'
#' @param model GAM/GAMM model for the outbreak
#' @param params model parameters for the outbreak model
#' @param time vector of time values at which the model will be evaluated
#' @param timedelta time step of the data described by \code{model}; note that this is the
#' time step of the original time series from which \code{model} was obtained, not the
#' (potentially different) time step at which model predictions are being evaluated
#' @return time series of cumulative case counts
#' @export
outbreak.calc.cum = function(model, params, time, timedelta=1) {
  # Get model predictions for given (randomized) param values
  predictors = model %>% predict(data.frame(time=time + timedelta / 2), type="lpmatrix")
  fit = predictors %*% params

  # Map spline fit back to data
  fit = fit %>% model$family$linkinv()

  # Calculate cumulative case fraction
  cumsum(fit) / sum(fit)
}

# Predict outbreak thresholds for an outbreak
#' This is useful as \code{quant} for outbreak.predict.scalars.
#'
#' @param model GAM/GAMM model for the outbreak
#' @param params model parameters for the outbreak model
#' @param time vector of time values at which the model will be evaluated
#' @param onset onset threshold as fraction of total outbreak case count
#' @param offset offset threshold as fraction of total outbreak case count
#' @param timedelta time step of the data described by \code{model}; note that this is the
#' time step of the original time series from which \code{model} was obtained, not the
#' (potentially different) time step at which model predictions are being evaluated
#' @return data frame with columns \code{onset} and \code{offset} representing time
#' when the outbreak crossed onsed and offset thresholds
#' @export
outbreak.calc.thresholds = function(model, params, time, onset=0.05, offset=0.95, timedelta=1) {
  # Calculate cumulative case counts from the model and parameters
  cumfit = outbreak.calc.cum(model, params, time, timedelta)

  # Compare to thresholds
  onsetTime = function(threshold) {
    found = any(cumfit < threshold)
    if (is.na(found)) {
      browser()
      head(time, 1)
    } else if (found) {
      tail(na.omit(time[cumfit < threshold]), 1)
    } else {
      head(time, 1)
    }
  }

  offsetTime = function(threshold) {
    found = any(cumfit > threshold)
    if (is.na(found)) {
      browser()
      tail(time, 1)
    } else if (found) {
      head(na.omit(time[cumfit > threshold]), 1)
    } else {
      tail(time, 1)
    }
  }

  data.frame(onset=onsetTime(onset), offset=offsetTime(offset))
}

