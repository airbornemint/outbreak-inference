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
#' a GAM/GAMM model of the outbreak, followed by calling either \code{\link{outbreak.estimate.scalars}} or
#' \code{\link{outbreak.estimate.timeseries}} to obtain confidence intervals on the
#' desired scalar/timeseries outcomes of the outbreak
#'
#' Both \code{\link{outbreak.estimate.scalars}} and \code{\link{outbreak.estimate.timeseries}}
#' allow simulation of arbitrary outbreak characteristics, by passing a function that calculates
#' the desired characteristics into \code{\link{outbreak.estimate.scalars}} or \code{\link{outbreak.estimate.timeseries}}.
#'
#' For convenience, this package also includes \code{\link{outbreak.calc.thresholds}}, which can
#' be use in conjunction with \code{\link{outbreak.estimate.scalars}} to estimate timing of outbreak onset and
#' offset, as well as \code{\link{outbreak.calc.cumcases}}, which can be used in conjuction with
#' \code{\link{outbreak.estimate.timeseries}} to estimate cumulative incidence vs time for the outbreak.
#'
#' @name outbreakinference-package
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
#' cumCases = outbreak.estimate.timeseries(
#'   model, modelTime,
#'   outbreak.calc.cumcases, level=.95
#')
#'
#' # Estimate time when outbreak crosses 5% and 95% of cumulative case count
#' thresholds = outbreak.estimate.scalars(
#'   model, modelTime,
#'   outbreak.calc.thresholds(onset=0.05, offset=0.95), level=.95
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
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
NULL

#' Sample from a multivariate normal distribution
#'
#' @param mu matrix of means
#' @param sig matrix of covariances
#' @param samples number of samples to draw
#' @return matrix of sampled values
#' @keywords internal
randomMVN = function(mu, sig, samples) {
  L = mroot(sig)
  m = ncol(L)
  t(mu + L %*% matrix(rnorm(m * samples), m, samples))
}

#' @keywords internal
outbreak.estimate.params = function(model, samples) {
  randomMVN(coef(model), model$Vp, samples)
}

#' @keywords internal
outbreak.estimate.scalars.sampleapply = function(params, model, time, outcomes) {
  params %>% apply(1, function(p) { outcomes(model, p, time) }) %>% bind_rows()
}

#' Runs simulations on an outbreak GAM/GAMM for the purpose of estimating
#' scalar outbreak outcomes, and returns estimated scalar outcome values for each simulation.
#'
#' This is mainly used internally by \code{outbreak.estimate.scalars}, but it's
#' useful if you want to calculate summary statistics of simulation results other than
#' the ones returned by \code{outbreak.estimate.scalars}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param time time values at which the model will be evaluated during simulation
#' @param outcomes scalar outcome generator function; see \code{\link{outbreak.estimate.scalars}} for more info
#' @param samples number of simulations to run
#' @return data frame with one row for each simulation; each row has \code{outbreak.sample} column giving
#' a unique ID for the simulation, and one column for each scalar outcome returned by \code{outcomes}
#'
#' @export
#' @keywords internal
outbreak.estimate.scalars.sample = function(model, time, outcomes, samples=100) {
  # Multivariate normal random generator
  # Generate random model parameters
  estimates = model %>% outbreak.estimate.params(samples) %>%
    outbreak.estimate.scalars.sampleapply(model, time, outcomes)

  if (any(is.na(estimates))) {
    warning("Some estimates are NA")
    return()
  }

  estimates %>% cbind(data.frame(outbreak.sample=seq(1, samples)))
}

#' Calculates confidence intervals for results of simulation performed by \code{\link{outbreak.estimate.scalars.sample}}
#'
#' @param samples data frame of samples as returned by \code{\link{outbreak.estimate.scalars.sample}}
#' @param level confidence level for calculated confidence intervals
#' @return data.frame of confidence intervals, as returned by \code{\link{outbreak.estimate.scalars}}
#' @export
#' @keywords internal
outbreak.estimate.scalars.confints = function(samples, level=.95) {
  confints = matrix(NA, nrow=1, ncol=0) %>% as.data.frame()

  for (name in names(samples)) {
    if (name == "outbreak.sample") {
      next
    }
    nameConfInts = as.data.frame(matrix(NA, 1, 3))
    names(nameConfInts) = paste(name, c("lower", "median", "upper"), sep=".")
    nameConfInts[1,] = quantile(samples[[name]], c((1 - level) / 2, 0.5, (1 + level) / 2), names=FALSE)
    confints = confints %>% cbind(nameConfInts)
  }

  confints
}

#' Calculates confidence intervals for scalar estimated from generalized additive (mixed) model of an outbreak
#'
#' This function performs a series of Monte Carlo simulations of a GAM/GAMM outbreak model.
#' For each simulated outbreak, it calls \code{outcomess} to calculate scalar
#' outcomes for the simulated outbreak (for example, the time of outbreak peak).
#' It then calculates and returns the confidence interval of each simulated scalar
#' outcome across all simulations.
#'
#' The \code{outcomes} function must accept (\code{model}, \code{params}, \code{time}) and return a one-row data frame
#' in which each column lists the value of a single scalar outcome calculated from the model
#' evaluated at the time points given in \code{time} and using the model parameters
#' given in \code{params}.
#'
#' A typical implementation of the \code{outcomes} function would call \code{predict} on
#' \code{model} and \code{time} to obtain the linear predictor matrix, and then post-multiply
#' that matrix by \code{params}. Having thus obtained model prediction at every time point,
#' it would calculate the desired scalar outcomes and return them in a data frame.
#'
#' For example, to calculate the time of outbreak peak, you might use this function for \code{outcomes}:
#'
#' \code{
#' calc_peak = function(model, params, time) {
#'   predictors = predict(model, data.frame(time=time), type="lpmatrix")
#'   fit = model$family$linkinv(predictors %*% params)
#'   data.frame(peak=time[which.max(fit)])
#' }
#' }
#'
#' The data frame returned by \code{outbreak.estimate.scalars} contains three columns for each
#' outcome calculated by \code{outcomes}: for outcome \code{x} returned by \code{outcomes},
#' \code{outbreak.estimate.scalars} returns columns \code{x.lower}, \code{x.median}, and \code{x.upper}, corresponding
#' to lower confidence limit, median, and upper confidence limit of \code{x}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param time vector of time values at which the model will be evaluated
#' @param outcomes function returning calculated scalar outcomes, as described above
#' @param samples number of simulations to run
#' @param level confidence level for estimates
#' @return data frame of estimates, as described above
#' @export
outbreak.estimate.scalars = function(model, time, outcomes, samples=100, level=.95) {
  model %>%
    outbreak.estimate.scalars.sample(time, outcomes, samples) %>%
    outbreak.estimate.scalars.confints(level)
}

#' @keywords internal
outbreak.estimate.timeseries.sampleapply = function(params, model, time, outcome) {
  apply(params, 1, function(p) { outcome(model, p, time) })
}

#' Runs simulations on an outbreak GAM/GAMM for the purpose of estimating time series outbreak outcomes, and
#' returns estimated time series outcomes for each simulation.
#'
#' This is mainly used internally by \code{outbreak.estimate.timeseries}, but
#' it's useful if you want to calculate summary statistics of simulation results other
#' than the ones returned by \code{outbreak.estimate.timeseries}.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}
#' @param time time values at which the model will be evaluated during simulation
#' @param outcome time series outcome generator function; see \code{\link{outbreak.estimate.timeseries}} for more info
#' @param samples number of simulations to run
#' @return matrix with one row for each simulation; each row contains the
#' time series calculated by \code{outcome} for the corresponding simulation run
#'
#' @export
#' @keywords internal
outbreak.estimate.timeseries.sample = function(model, time, outcome, samples=100) {
  # Generate random model parameters
  estimates = model %>% outbreak.estimate.params(samples) %>%
    outbreak.estimate.timeseries.sampleapply(model, time, outcome)

  if (any(is.na(estimates))) {
    warning("Model did not converge, try reducing the number of spline knots in GAM")
    return()
  }

  estimates
}

#' Calculates confidence intervals for results of simulation performed by \code{\link{outbreak.estimate.timeseries.sample}}
#'
#' @param samples data frame of samples as returned by \code{\link{outbreak.estimate.timeseries.sample}}
#' @param level confidence level for calculated confidence intervals
#' @return data.frame of confidence intervals, as returned by \code{\link{outbreak.estimate.timeseries}}
#' @export
#' @keywords internal
outbreak.estimate.timeseries.confints = function(samples, level=0.95) {
  resultNames = c("lower", "median", "upper")
  if (is.null(samples)) {
    confints = as.data.frame(matrix(NA, 0, 3))
    names(confints) = resultNames
  } else {
    confints = apply(samples, 1, function(sample) {
      ci = as.data.frame(matrix(NA, 0, 3))
      names(ci) = resultNames

      ci[1,] = quantile(sample, c(1-level, 0.5, level), names=FALSE)
      ci
    }) %>% bind_rows()
  }

  confints
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
#' The data frame returned by \code{outbreak.estimate.timeseries} contains three columns and
#' one row for each time point in \code{time}. The columns are \code{lower}, \code{median}, and
#' \code{upper}, containing the median and the confidence interval for the computed
#' outcome time series at each time point.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param time vector of time values at which the model will be evaluated
#' @param outcome function returning calculated outcome time series, as described above
#' @param samples number of simulations to run
#' @param level confidence level for returned estimates
#' @return data frame of estimates, as described above
#' @export
outbreak.estimate.timeseries = function(model, time, outcome, samples=1000, level=.95) {
  model %>%
    outbreak.estimate.timeseries.sample(time, outcome, samples) %>%
    outbreak.estimate.timeseries.confints(level) %>%
		cbind(time=time)
}

#' Calculate case count for an outbreak
#'
#' This is useful as \code{outcome} for outbreak.estimate.scalars.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param params model parameter matrix
#' @param time vector of time values at which the model will be evaluated
#' @return time series of cumulative case counts
#' @export
outbreak.calc.cases = function(model, params, time) {
  # Get model predictions for given (randomized) param values
  predictors = model %>% predict(data.frame(time=time), type="lpmatrix")
  fit = predictors %*% params

  # Map spline fit back to data
  fit %>% model$family$linkinv()
}

#' Calculate cumulative case count for an outbreak
#'
#' This is useful as \code{outcome} for outbreak.estimate.timeseries.
#'
#' @param model model returned by \code{\link[mgcv]{gam}} or \code{\link[mgcv]{gamm}}, with a single parameter (time)
#' @param params model parameter matrix
#' @param time vector of time values at which the model will be evaluated
#' @return time series of relative cumulative cases
#' @export
outbreak.calc.cumcases = function(model, params, time) {
  # Get model predictions for given (randomized) param values
  predictors = model %>% predict(data.frame(time=time), type="lpmatrix")
  fit = predictors %*% params

  # Map spline fit back to data
  fit = fit %>% model$family$linkinv()

  # Calculate cumulative case fraction
  result = cumsum(fit) / sum(fit)

  if (any(is.nan(result))) {
    stop("cumcases diverged")
  }

  result
}

# Calculate outbreak thresholds for an outbreak
#' This is useful as \code{outcomes} for outbreak.estimate.scalars.
#'
#' @param onset onset threshold as fraction of total outbreak case count
#' @param offset offset threshold as fraction of total outbreak case count
#' @return data frame with columns \code{onset} and \code{offset} representing time
#' when the outbreak crossed onsed and offset thresholds
#' @export
outbreak.calc.thresholds = function(onset=NA, offset=NA) {
  function(model, params, time) {
    # Calculate cumulative case counts from the model and parameters
    cumfrac = outbreak.calc.cumcases(model, params, time)

    # Compare to threshold
    thresholdTime = function(threshold) {
      # Linear interpolation across the time interval that where cumfit crosses threshold
      idx = sum(cumfit < threshold)
      timeLow = time[idx]
      timeHigh = time[idx + 1]
      cumLow = cumfit[idx]
      cumHigh = cumfit[idx + 1]
      return(timeLow + (timeHigh - timeLow) / (cumHigh - cumLow) * (threshold - cumLow))
    }

		if (is.na(offset)) {
	    data.frame(
	      onset=outbreak.calc.threshold.ts(time, cumfrac, onset)
      )
		} else if (is.na(onset)) {
	    data.frame(
	      offset=outbreak.calc.threshold.ts(time, cumfrac, offset)
      )
		} else {
	    data.frame(
	      onset=outbreak.calc.threshold.ts(time, cumfrac, onset),
	      offset=outbreak.calc.threshold.ts(time, cumfrac, offset)
      )
		}
  }
}

outbreak.calc.threshold.ts = function(time, cumfrac, threshold) {
  # Linear interpolation across the time interval that where cumfrac crosses threshold
  idx = sum(cumfrac < threshold)
  timeLow = time[idx]
  timeStep = time[idx + 1] - time[idx]
  cumLow = cumfrac[idx]
  cumStep = cumfrac[idx + 1] - cumfrac[idx]
  return(timeLow + timeStep / cumStep * (threshold - cumLow))
}

