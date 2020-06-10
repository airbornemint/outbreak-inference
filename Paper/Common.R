import::from(dplyr, rename, mutate)
import::from(plyr, adply)
import::from(data.table, setnames)
import::from(reshape, melt)
import::from(scales, percent)
import::from(stringr, str_replace_all)

predCases = function(params, model, predictors) {
  params %>%
    pspline.inference:::pspline.calc.timeseries(model, predictors, pspline.outbreak.cases)
}

predCum = function(params, model, predictors) {
  params %>%
    pspline.inference:::pspline.calc.timeseries(model, predictors, pspline.outbreak.cumcases.relative)
}

predThresholds = function(params, model, predictors, seasonThreshold) {
  thresholds = params %>%
    pspline.inference:::pspline.calc.scalars(model, predictors, pspline.outbreak.thresholds(seasonThreshold, 1 - seasonThreshold)) %>%
    adply(1, function(row) {
      cases = params[row$pspline.sample,] %>%
        t() %>%
        pspline.inference:::pspline.calc.timeseries(model, data.frame(time=c(row$onset, row$offset)), pspline.outbreak.cases)
      data.frame(onset.cases=cases$cases[1], offset.cases=cases$cases[2]) %>%
        cbind(row %>% select(onset, offset))
    })
}

aapStrat = function(t) {
  as.numeric(t >= 20 & t <= 44)
}

calcFraction = function(model, data) {
  totalCases = sum(data$cases)
  preventableCases = sum(data$cases[aapStrat(data$time)])
  data.frame(preventable=c(preventableCases / totalCases))
}

predFraction = function(params, model, predictors) {
  params %>%
    pspline.inference:::pspline.calc.scalars(model, predictors, calcFraction)
}

latexPercent = function(...) {
  percent(...) %>% str_replace_all("%", "\\\\%")
}
