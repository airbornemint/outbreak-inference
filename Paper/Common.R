import::from(dplyr, rename, mutate)
import::from(data.table, setnames)
import::from(reshape, melt)

predCases = function(params, model, modelTime) {
  pred = params %>%
    outbreakinference:::outbreak.estimate.timeseries.sampleapply(model, modelTime, outbreak.calc.cases) %>%
    data.frame()

  col = ncol(pred)

  pred = pred %>%
    setnames(as.character(seq(1:col))) %>%
    cbind(time=modelTime) %>%
    melt(c("time")) %>%
    dplyr::rename(sim=variable, cases=value) %>%
    dplyr::mutate(sim=as.numeric(sim))
}

predCum = function(params, model, modelTime) {
  pred = params %>%
    outbreakinference:::outbreak.estimate.timeseries.sampleapply(model, modelTime, outbreak.calc.cumcases) %>%
    data.frame()

  col = ncol(pred)

  pred = pred %>%
    setnames(as.character(seq(1:col))) %>%
    cbind(time=modelTime) %>%
    melt(c("time")) %>%
    dplyr::rename(sim=variable, cases.cum.frac=value) %>%
    dplyr::mutate(sim=as.numeric(sim))
}

predThresholds = function(params, model, modelTime, seasonThreshold) {
  thresholds = params %>%
    outbreakinference:::outbreak.estimate.scalars.sampleapply(model, modelTime, outbreak.calc.thresholds(seasonThreshold, 1-seasonThreshold)) %>%
    select(onset, offset)
  thresholds$onset.cases = params %>%
    outbreakinference:::outbreak.estimate.timeseries.sampleapply(model, thresholds$onset, outbreak.calc.cases) %>%
    diag()
  thresholds$offset.cases = params %>%
    outbreakinference:::outbreak.estimate.timeseries.sampleapply(model, thresholds$offset, outbreak.calc.cases) %>%
    diag()
  thresholds
}

ppxStart = 20
ppxEnd = 20 + 24

aapStrat = function(t) {
  as.numeric(t >= ppxStart & t <= ppxEnd)
}

calcFraction = function(model, params, time) {
  predictors = predict(model, data.frame(time=time), type="lpmatrix")
  cases = model$family$linkinv(predictors %*% params)

  cases = data.frame(time=time, cases=cases)
  totalCases = sum(cases$cases)
  preventableCases = sum(cases$cases[cases$time >= ppxStart & cases$time <= ppxEnd])

  data.frame(preventable=c(preventableCases / totalCases))
}

predFraction = function(params, model, modelTime) {
  params %>%
    outbreakinference:::outbreak.estimate.scalars.sampleapply(model, modelTime, calcFraction)
}

