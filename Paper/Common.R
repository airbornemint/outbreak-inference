predCases = function(params, model, modelTime) {
  pred = params %>%
    outbreakinference:::outbreak.estimate.timeseries.sampleapply(model, modelTime, outbreak.calc.cases) %>%
    data.frame()

  col = ncol(pred)

  pred = pred %>%
    setnames(as.character(seq(1:col))) %>%
    cbind(time=modelTime) %>%
    melt(c("time")) %>%
    rename(sim=variable, cases=value) %>%
    mutate(sim=as.numeric(sim))
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
    rename(sim=variable, cases.cum.frac=value) %>%
    mutate(sim=as.numeric(sim))
}

predOnset = function(params, model, modelTime) {
  params %>%
    outbreakinference:::outbreak.estimate.scalars.sampleapply(model, modelTime, outbreak.calc.thresholds(seasonThreshold, 1-seasonThreshold)) %>%
    select(onset)
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

