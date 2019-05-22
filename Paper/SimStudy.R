# ---- sim ----

import::from(outbreakinference, outbreak.calc.cumcases, outbreak.calc.cases, outbreak.calc.thresholds, outbreak.estimate.scalars)

import::from(dplyr, "%>%", rename, mutate, bind_rows, select, filter, group_by, ungroup, do, arrange, first, last, inner_join)
import::from(plyr, rlply, llply, ldply)
import::from(mgcv, gam)
import::from(doParallel, registerDoParallel)
import::from(parallel, detectCores)

registerDoParallel(detectCores())
set.seed(NULL)

source("./Paper/Common.R")

# Generate observed values of outcome measures
simExpectedOutcome = function(truth, onsetThreshold, offsetThreshold) {
  peakRange = range(truth$time[truth$cumfrac > onsetThreshold & truth$cumfrac < offsetThreshold])
  data.frame(
    onset=outbreakinference:::outbreak.calc.threshold.ts(truth$time, truth$cumfrac, onsetThreshold),
    offset=outbreakinference:::outbreak.calc.threshold.ts(truth$time, truth$cumfrac, offsetThreshold)
  )
}

# Calculate estimated values of outcome measures
runObsOutcome = function(model, modelTime, onsetThreshold, offsetThreshold, n, level) {
  tryCatch({
    model %>% outbreak.estimate.scalars(
      modelTime,
      outbreak.calc.thresholds(onset=onsetThreshold, offset=offsetThreshold),
      samples=n, level=level)
  }, error=function(e) {
    data.frame(onset.lower=NA, onset.upper=NA, offset.lower=NA, offset.upper=NA)
  })
}

runObsOutcomeDetail = function(model, modelTime, onsetThreshold, offsetThreshold, n, level, seed) {
  assign('.Random.seed', seed, envir = .GlobalEnv)
  params = model %>% outbreakinference:::outbreak.estimate.params(n)
  cases = params %>% predCases(model, modelTime)
  cumcases = params %>% predCum(model, modelTime)
  estimates = params %>% outbreakinference:::outbreak.estimate.scalars.sampleapply(model, modelTime, outbreak.calc.thresholds(onset=onsetThreshold, offset=offsetThreshold)) %>%
    cbind(data.frame(outbreak.sample=seq(1, n)))
  outcome = estimates  %>% outbreak.estimate.scalars.confints(level)
  list(cases=cases, cumcases=cumcases, estimates=estimates, outcome=outcome)
}

simResults = function(simRunCount, simTime, obsTime) {
  # Generate something that is roughly outbreak-shaped. We'll smooth it with splines and then use it as outbreak truth.
  outbreakTemplate = function(simTime) {
    # For our purposes, pick a random onset/offset time and peak amplitude
    onsetT = runif(1, 5, 15)
    offsetT = runif(1, 30, 45)
    peakT = (onsetT + offsetT) / 2
    peak = runif(1, 50, 500)

    # Then use a single period of a sine wave for the outbreak peak
    cases = log(peak) / 2 * (1 - cos(2 * pi * (simTime - onsetT) / (offsetT - onsetT))) * (onsetT <= simTime & simTime < offsetT)
    cases = round(exp(cases) - 1)

    data.frame(time=simTime, cases=cases)
  }

  makeModel = function(data) {
    gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=data)
  }

  trueCases = function(simTime) {
    seed = .Random.seed
    # Generate a template outbreak
    template = outbreakTemplate(simTime)

    # Approximate it with splines to get "real" outbreak truth
    model = makeModel(template)
    cases = model %>% predict(type="response")

    cumcases = cumsum(cases)
    list(seed=seed, truth=data.frame(time=simTime, cases=cases, cumfrac=cumcases / tail(cumcases, 1)))
  }

  # Classify outcomes based on whether expected outcome is within the CI of the estimated outcome
  classifyOutcomes = function(expected, observed) {
    expected %>% cbind(observed) %>%
      mutate(
        onsetGood = observed$onset.lower <= expected$onset & expected$onset <= observed$onset.upper,
        offsetGood = observed$offset.lower <= expected$offset & expected$offset <= observed$offset.upper
      )
  }

  # Generate observed number of cases vs time
  obsCases = function(truth, obsTime) {
    data.frame(time=truth$time, cases=truth$cases) %>%
      inner_join(data.frame(time=obsTime), by=c("time")) %>%
      mutate(cases=rpois(length(cases), cases))
  }

  onsetThreshold = 0.025
  offsetThreshold = 1 - onsetThreshold
  n = 2000
  cl = 0.95

  simRuns = function(truth, simTime, obsTime) {
    expected = truth %>% simExpectedOutcome(onsetThreshold, offsetThreshold)
    llply(
      1:simRunCount, function(idx) {
        observedCases = obsCases(truth, obsTime)
        model = gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=observedCases)
        seed = .Random.seed
        observed = model %>% runObsOutcome(simTime, onsetThreshold, offsetThreshold, n, cl)
        list(runIdx=idx, seed=seed,
             observedCases=observedCases, model=model,
             n=n, cl=cl,
             outcomes=classifyOutcomes(expected, observed)
        )
      }, .parallel=TRUE
    )
  }

  sim = simTime %>% trueCases
  results = sim$truth %>% simRuns(simTime, obsTime)
  c(sim, list(onsetThreshold=onsetThreshold, offsetThreshold=offsetThreshold, runs=results))
}

simStudy = function(simCount, simRunCount, eps) {
	obsTime = seq(1, 52)
	simTime = seq(min(obsTime), max(obsTime) + 1 - eps, by=eps)

	# Main result calculation
	results = llply(1:simCount, function(idx) { c(list(simIdx=idx), simResults(simRunCount, simTime, obsTime)) }, .progress="text")

	# Some useful summaries of results
	byRun = results %>%
	  llply(function(studyResult) {
	    studyResult$runs %>%
	      llply(function(runResult) {
	        runResult$outcomes %>% mutate(runIdx=runResult$runIdx, simIdx=studyResult$simIdx)
	      }) %>%
	      bind_rows()
  	}) %>%
	  bind_rows()

	summary = data.frame(
	  onset=mean(byRun$onsetGood[!is.na(byRun$onsetGood)]),
	  offset=mean(byRun$offsetGood[!is.na(byRun$offsetGood)])
	)

	list(results=results, byRun=byRun, summary=summary)
}

# Use this to repeat a specific run for debugging purpose. Returns more detailed information about the run.
repeatRun = function(results, simIdx, runIdx) {
  sim = results$results[[simIdx]]
  run = sim$runs[[runIdx]]

  expected = sim$truth %>% simExpectedOutcome(sim$onsetThreshold, sim$offsetThreshold)
  observedDetail = runObsOutcomeDetail(run$model, sim$truth$time, sim$onsetThreshold, sim$offsetThreshold, run$n, run$cl, run$seed)
  c(list(expected=expected), observedDetail)
}

plotRun = function(results, simIdx, runIdx) {
  results = results$results
  study = results[[simIdx]]
  run = study$runs[[runIdx]]

  plot(study$truth$time, study$truth$cases, type="l")
  observed = run$observedCases
  points(observed$time, observed$cases, new=TRUE)
}

studyResults = simStudy(1, 5, 0.05)
#studyResults = simStudy(60, 60, 0.05)

onsetFailResults = studyResults$byRun %>% filter(!onsetGood) %>% mutate(onsetError=ifelse(onset > onset.upper, onset-onset.upper, onset-onset.lower))
# hist(onsetFailResults$onset - onsetFailResults$onset.median, breaks=40)

offsetFailResults = studyResults$byRun %>% filter(!offsetGood) %>% mutate(offsetError=ifelse(offset > offset.upper, offset-offset.upper, offset-offset.lower))
# hist(offsetFailResults$offset - offsetFailResults$offset.median, breaks=40)

tryadj = function(delta) {
  adj = studyResults$byRun %>% mutate(onset.lower=onset.lower+delta, onset.upper=onset.upper+delta) %>% mutate(onsetGood=(onset>=onset.lower & onset<=onset.upper))
  data.frame(delta=delta, onset=mean(adj$onsetGood[!is.na(adj$onsetGood)]))
}

# ldply(seq(-.1, .1, length.out=201), tryadj)

#runRepeat = repeatRun(studyResults1, 1, 2)

studyResults$summary

saveRDS(studyResults, "studyResults.rds")
