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

randomMVN = function(mu, sig, nsim) {
  L = mroot(sig)
  m = ncol(L)
  t(mu + L %*% matrix(rnorm(m * nsim), m, nsim))
}

sampleCalcPred = function(sampleParams, sampleNsim) {
  sampleParams %>%
    apply(1, function(params) { outbreak.calc.cum(1)(sampleModel, params, modelTime) } ) %>%
    data.frame() %>%
    setnames(as.character(seq(1:sampleNsim))) %>%
    cbind(time=modelTime) %>%
    melt(c("time")) %>%
    rename(sim=variable, rsv.cum.frac=value) %>%
    mutate(sim=as.numeric(sim))  
}

sampleCalcOnset = function(sampleParams) {
  sampleParams %>% 
    apply(1, function(params) { outbreakinference::outbreak.calc.thresholds(seasonThreshold, 1-seasonThreshold)(sampleModel, params, modelTime) } ) %>%
    bind_rows() %>%
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

sampleCalcFraction = function(samplePred) {
  samplePred %>%
    mutate(ppx=aapStrat(time)) %>%
    filter(ppx > 0) %>%
    group_by(sim) %>%
    do((function(df) {
      df = df %>% arrange(time)
      data.frame(
        ppx.start=min(df$time),
        ppx.end=max(df$time),
        unprotected.start=first(df$rsv.cum.frac),
        unprotected.end=last(df$rsv.cum.frac)
      ) %>% mutate(
        unprotected=unprotected.end-unprotected.start
      )
    })(.)) %>%
    ungroup() %>%
    as.data.frame()  
}

sampleParamsSingle = randomMVN(coef(sampleModel), sampleModel$Vp, 1)
samplePredSingle = sampleParamsSingle %>% sampleCalcPred(1)
sampleFractionSingle = samplePredSingle %>% sampleCalcFraction()

sampleNsimMini = 5
sampleParamsMini = randomMVN(coef(sampleModel), sampleModel$Vp, sampleNsimMini)
samplePredMini = sampleParamsMini %>% sampleCalcPred(sampleNsimMini)
sampleOnsetMini = sampleParamsMini %>% sampleCalcOnset()

zoomedStartWeek = min(monthBoundaries$min) + 5
zoomedEndWeek = zoomedStartWeek + 21

sampleNsimFull = simulations
sampleParamsFull = randomMVN(coef(sampleModel), sampleModel$Vp, sampleNsimFull)
samplePredFull = sampleParamsFull %>% sampleCalcPred(sampleNsimFull)
sampleOnsetFull = sampleParamsFull %>% sampleCalcOnset()
sampleFractionFull = samplePredFull %>% sampleCalcFraction()

sampleDisplayNsim = 100