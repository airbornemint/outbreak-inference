# ---- paper ----
import::from(outbreakinference, outbreak.calc.cum)

import::from(reshape, melt)
import::from(dplyr, "%>%", rename, mutate, bind_rows, select, filter, group_by, ungroup, do, arrange, first, last, inner_join)
import::from(mgcv, gam, mroot)
import::from(data.table, setnames)

sampleObs = read.csv("../vignettes/rsv.csv")
sampleObs$rsv.cum = cumsum(sampleObs$rsv)
sampleObs$rsv.cum.frac = sampleObs$rsv.cum / max(sampleObs$rsv.cum)

simulations = 2000
eps = .05
seasonThreshold = 0.025
ppxDuration = 24 # weeks
startYear = 1996
endYear = 2013

rsvTime = sampleObs$time
modelTime = seq(min(rsvTime) - 1 + eps, max(rsvTime), eps)

sampleModel = gam(rsv ~ s(time, k=20, bs="cp", m=3), family=poisson, data=sampleObs)

predCases = function(params, model) {
  pred = params %>%
    outbreakinference:::outbreak.predict.timeseries.simapply(model, modelTime, outbreak.calc.cases()) %>%
    data.frame()

  col = ncol(pred)

  pred = pred %>%
    setnames(as.character(seq(1:col))) %>%
    cbind(time=modelTime) %>%
    melt(c("time")) %>%
    rename(sim=variable, rsv=value) %>%
    mutate(sim=as.numeric(sim))
}

predCum = function(params, model) {
  pred = params %>%
    outbreakinference:::outbreak.predict.timeseries.simapply(model, modelTime, outbreak.calc.cum()) %>%
    data.frame()

  col = ncol(pred)

  pred = pred %>%
    setnames(as.character(seq(1:col))) %>%
    cbind(time=modelTime) %>%
    melt(c("time")) %>%
    rename(sim=variable, rsv.cum.frac=value) %>%
    mutate(sim=as.numeric(sim))
}

predOnset = function(params, model) {
  params %>%
    outbreakinference:::outbreak.predict.scalars.simapply(model, modelTime, outbreak.calc.thresholds(seasonThreshold, 1-seasonThreshold)) %>%
    select(onset)
}

evalStrategy = function(start, end, time) {
  as.numeric((time >= start) & (time < end))
}
aapStart = 20
aapStrat = function(time) {
  evalStrategy(aapStart, aapStart + ppxDuration, time)
}

# .01 to let minor gridlines show through
epiWeekBreaks = c(3.01, 7.25, 11.75, 16.01, 20.25, 24.75, 29.01, 33.01, 37.5, 42.01, 46.25, 50.75)
epiWeekLabels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun")

breaks.df = data.frame(i=seq(1, 12), mid=epiWeekBreaks)
monthBoundaries = breaks.df %>%
  inner_join(
    breaks.df %>% mutate(i=i %% 12 + 1) %>% rename(prevMid=mid),
    by="i"
  ) %>%
  inner_join(
    breaks.df %>% mutate(i=(i - 2) %% 12 + 1) %>% rename(nextMid=mid),
    by="i"
  ) %>%
  mutate(
    min=(prevMid + (mid - prevMid) %% 52 / 2) %% 52,
    max=mid + (nextMid - mid) %% 52 / 2
  ) %>%
  select(i, min, mid, max)

calcFraction = function(model, params, time) {
  predictors = predict(model, data.frame(time=time), type="lpmatrix")
  cases = model$family$linkinv(predictors %*% params)

  cases = data.frame(time=time, cases=cases)
  totalCases = sum(cases$cases)
  preventableCases = sum(cases$cases[cases$time >= ppxStart & cases$time <= ppxEnd])

  data.frame(preventable=c(preventableCases / totalCases))
}

predFraction = function(params, model) {
  params %>%
    outbreakinference:::outbreak.predict.scalars.simapply(model, modelTime, calcFraction)
}

sampleParamsSingle = sampleModel %>% outbreakinference:::outbreak.predict.params(1)
predCumSingle = sampleParamsSingle %>% predCum(sampleModel)
predFractionSingle = sampleParamsSingle %>% predFraction(sampleModel)

sampleNsimMini = 5
sampleParamsMini = sampleModel %>% outbreakinference:::outbreak.predict.params(sampleNsimMini)
predCasesMini = sampleParamsMini %>% predCases(sampleModel)
predCumMini = sampleParamsMini %>% predCum(sampleModel)
predOnsetMini = sampleParamsMini %>% predOnset(sampleModel)

zoomedStartWeek = min(monthBoundaries$min) + 5
zoomedEndWeek = zoomedStartWeek + 21

sampleNsimFull = simulations
sampleParamsFull = sampleModel %>% outbreakinference:::outbreak.predict.params(sampleNsimFull)
predCumFull = sampleParamsFull %>% predCum(sampleModel)
predOnsetFull = sampleParamsFull %>% predOnset(sampleModel)
predFractionFull = sampleParamsFull %>% predFraction(sampleModel)

sampleDisplayNsim = 100
