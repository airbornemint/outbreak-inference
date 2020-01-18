# ---- sim ----

import::from(pspline.inference, pspline.outbreak.thresholds, pspline.validate.scalars)

import::from(dplyr, "%>%", mutate, inner_join)
import::from(mgcv, gam)
import::from(doParallel, registerDoParallel)
import::from(parallel, detectCores)

registerDoParallel(detectCores())
message(paste("foreach::dopar", foreach::getDoParRegistered(), foreach::getDoParName(), foreach::getDoParVersion(), foreach::getDoParWorkers(), detectCores()))
set.seed(NULL)

source("./Common.R")

onsetSimMin = 5
onsetSimMax = 15
offsetSimMin = 30
offsetSimMax = 45
peakSimMin = 50
peakSimMax = 500

# Calculate estimates of outcome measures
outcomes = function(onsetThreshold, offsetThreshold) {
  thresholds = pspline.outbreak.thresholds(onset=onsetThreshold, offset=offsetThreshold)
  function(model, data) {
    thresholds(model, data)
  }
}

generateTruth = function(time.min, time.max, time.delta) {
  function() {
    time = seq(time.min - 0.5, time.max + 0.5 - time.delta, time.delta)

    # Pick a random onset/offset time and peak amplitude
    onsetT = runif(1, onsetSimMin, onsetSimMax)
    offsetT = runif(1, offsetSimMin, offsetSimMax)
    peak = runif(1, peakSimMin, peakSimMax)

    # Then use a single period of a sine wave for the outbreak peak
    cases0 = log(peak) / 2 * (1 - cos(2 * pi * (time - onsetT) / (offsetT - onsetT))) * (onsetT <= time & time < offsetT)
    cases0 = round(exp(cases0) - 1)
    data0 = data.frame(time=time, cases=cases0)

    # And approximate it with a spline
    model = gam(cases ~ s(time, k=4, bs="cp", m=3), family=poisson, data=data0)
    data.frame(
      time=time,
      cases=model %>% predict(type="response")
    )
  }
}

generateObservations = function(truth) {
  seq(ceiling(min(truth$time)), floor(max(truth$time))) %>%
    data.frame(time=.) %>%
    inner_join(truth, by="time") %>%
    mutate(cases=rpois(length(cases), cases))
}

makeModel = function(data) {
  gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=data)
}

simNTruth = 60
simNObs = 60

set.seed(NULL)
if(getOption("pspline.paper.validation.run", FALSE)) {
  validationResults = pspline.validate.scalars(
    generateTruth(1, 52, 0.05), simNTruth,
    generateObservations, simNObs,
    makeModel, outcomes(onsetThreshold=0.025, offsetThreshold=0.975), 2000, 0.95
  )

  saveRDS(validationResults, "ValidationResults.rds")
} else {
  # Even when skipping full validation, run one cycle of it just to make sure it still works
  validationResults = pspline.validate.scalars(
    generateTruth(1, 52, 0.05), 20,
    generateObservations, 20,
    makeModel, outcomes(onsetThreshold=0.025, offsetThreshold=0.975), 20, 0.95
  )

#  validationResults = readRDS("ValidationResults.rds")
}

tikz(sprintf("%s/onsetQuantilesFull.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

ggplot(validationResults$results) +
  theme_light(base_size=plotTextBaseSize) +
  geom_histogram(aes(x=onset.quantile), binwidth=0.025, fill="gray") +
  coord_cartesian(xlim=c(0, 1)) +
  scale_x_continuous(labels = latexPercent) +
  labs(x=NULL, y="Count")


dev.off()

tikz(sprintf("%s/offsetQuantilesFull.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

ggplot(validationResults$results) +
  theme_light(base_size=plotTextBaseSize) +
  geom_histogram(aes(x=offset.quantile), binwidth=0.025, fill="gray") +
  coord_cartesian(xlim=c(0, 1)) +
  scale_x_continuous(labels = latexPercent) +
  labs(x=NULL, y="Count")


dev.off()


