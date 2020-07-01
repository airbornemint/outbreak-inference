# ---- paper ----
import::from(pspline.inference, pspline.outbreak.cumcases.relative, pspline.outbreak.cases, pspline.outbreak.thresholds, pspline.estimate.scalars, pspline.estimate.timeseries, pspline.confints.scalars)

import::from(dplyr, "%>%", rename, mutate, select, filter, group_by, ungroup, do, arrange, first, last, inner_join)
import::from(mgcv, gam)

source("./Common.R")

obs = read.csv(system.file("extdata", "seasonal.csv", package="pspline.inference"))
obs$cases.cum = cumsum(obs$cases)
obs$cases.cumrel = obs$cases.cum / max(obs$cases.cum)

simulations = 2000
eps = .05
seasonThreshold = 0.025
ppxDuration = 24 # weeks
startYear = 1996
endYear = 2013
level = 0.95

obsTime = obs$time

model = gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=obs)
predictors = data.frame(time=seq(min(obsTime) - 0.5, max(obsTime) + 0.5 - eps, eps))

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

predCasesBest = model %>% pspline.estimate.timeseries(predictors, pspline.outbreak.cases)

sampleParamsSingle = model %>% pspline.inference:::sample.params(1)
predCumSingle = sampleParamsSingle %>% predCum(model, predictors)
predFractionSingle = sampleParamsSingle %>% predFraction(model, predictors)

sampleNsimMini = 5
sampleParamsMini = model %>% pspline.inference:::sample.params(sampleNsimMini)
predCasesMini = sampleParamsMini %>% predCases(model, predictors)
predCumMini = sampleParamsMini %>% predCum(model, predictors)
predThresholdsMini = sampleParamsMini %>% predThresholds(model, predictors, seasonThreshold)

zoomedStartWeek = min(monthBoundaries$min) + 5
zoomedEndWeek = zoomedStartWeek + 21

sampleNsimFull = simulations
sampleParamsFull = model %>% pspline.inference:::sample.params(sampleNsimFull)
predCumFull = sampleParamsFull %>% predCum(model, predictors)
predThresholdsFull = sampleParamsFull %>% predThresholds(model, predictors, seasonThreshold)
estThresholdsFull = predThresholdsFull %>% select(-onset.cases, -offset.cases) %>% pspline.confints.scalars(model, level)
predFractionFull = sampleParamsFull %>% predFraction(model, predictors)
estFractionFull = predFractionFull %>% pspline.confints.scalars(model, level)

sampleDisplayNsim = 100

# ---- figures ----

import::from(ggplot2, ggplot, theme_light, geom_point, aes, scale_x_continuous, scale_y_continuous, geom_line, geom_segment, theme, coord_cartesian, element_blank, element_text, geom_rect, geom_text, geom_violin, geom_density, labs, sec_axis, geom_histogram, scale_color_identity, scale_size_identity, scale_shape_identity, margin, unit, element_rect, scale_fill_identity)
import::from(dplyr, "%>%", filter, mutate, do, arrange, first, last, ungroup)
import::from(ggstance, geom_violinh)
import::from(gridExtra, grid.arrange)
import::from(grDevices, dev.off)
import::from(tikzDevice, tikz)
import::from(magrittr, "%>%")

inlinePlotWidth = 3.1
inlinePlotHeight = 2.5
pagePlotWidth = inlinePlotWidth * 1.8
pagePlotHeight = inlinePlotHeight * 1.8 * 2 / 5
plotTextBaseSize = 8

figuresDir = paste0(getOption("pspline.paper.output", "."), "/figures")
if (!dir.exists(figuresDir)) {
  dir.create(figuresDir, recursive=TRUE)
}
options(
  tikzMetricsDictionary='./tikzDictionary.dat',
  tikzDocumentDeclaration='\\documentclass[tikz]{standalone}',
  tikzLatexPackages=c(
    "\\usepackage{tikz}\n"
  )
)

commonOptions = list(
  scale_x_continuous(expand=c(0, 0)),
  scale_y_continuous(),
  labs(x="Time (weeks)", y="RSV incidence"),
  coord_cartesian(xlim=range(c(monthBoundaries$min, monthBoundaries$max))),
  theme_light(base_size=plotTextBaseSize) +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position=c(0.025, 0.975),
    legend.justification = c("left", "top"),
    legend.margin = margin(0, 2, 0, 2),
    legend.box.margin = margin(0, 2, 0, 2),
    legend.box.spacing = unit(0, "in"),
    legend.spacing = unit(0, "in"),
    legend.key.height = unit(0.25, "line")
  )
)

tikz(sprintf("%s/sampleCases.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

ggplot(obs) +
  geom_point(aes(x=time, y=cases, shape=16), size=.5) +
  scale_shape_identity(guide="legend", breaks=c(16), labels=("Observed"), name=NULL) +
  commonOptions

dev.off()

tikz(sprintf("%s/sampleBestFit.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

ggplot(obs) +
  geom_line(data=predCasesBest, aes(x=time, y=cases.median, color="black"), size=0.2) +
  geom_point(aes(x=time, y=cases, size=0.5), color="black", shape=1) +
  scale_color_identity(guide="legend", breaks=c("black"), labels=c("Model"), name=NULL) +
  scale_size_identity(guide="legend", breaks=c(0.5), labels=c("Observed"), name=NULL) +
  commonOptions

dev.off()

tikz(sprintf("%s/samplePredMini.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

ggplot(predCasesMini) +
  geom_line(aes(x=time, y=cases, group=pspline.sample, color="black"), size=0.2) +
  geom_point(data=obs, aes(x=time, y=cases, size=0.5), color="black", shape=1) +
  scale_color_identity(guide="legend", breaks=c("black"), labels=c("Model"), name=NULL) +
  scale_size_identity(guide="legend", breaks=c(0.5), labels=c("Observed"), name=NULL) +
  commonOptions

dev.off()

tikz(sprintf("%s/samplePredOnsetMini.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

ggplot(predCasesMini) +
  geom_line(aes(x=time, y=cases, group=pspline.sample, size=0.2), color="black") +
  geom_point(data=predThresholdsMini, aes(x=onset, y=0, shape=17), size=0.5) +
  geom_segment(data=predThresholdsMini, aes(x=onset, xend=onset, y=onset.cases, yend=0), size=0.1) +
  geom_point(data=predThresholdsMini, aes(x=offset, y=0, shape=17), size=0.5) +
  geom_segment(data=predThresholdsMini, aes(x=offset, xend=offset, y=offset.cases, yend=0), size=0.1) +
  scale_shape_identity(guide="legend", breaks=c(17), labels=c("Onset / offset"), name=NULL) +
  scale_size_identity(guide="legend", breaks=c(0.2), labels=c("Model"), name=NULL) +
  commonOptions

dev.off()

tikz(sprintf("%s/samplePredOnsetFull.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10, standAlone = TRUE)

data = predCumFull
splineData = data %>% filter(pspline.sample < sampleDisplayNsim)

ggplot(predThresholdsFull) +
  # geom_line(data=splineData, aes(x=time, y=cases, group=pspline.sample), color="gray") +
  # geom_segment(aes(x=-Inf, y=seasonThreshold, xend=+Inf, yend=seasonThreshold), linetype="11", data=data.frame(), size=.375) +
  geom_line(data=predCasesBest, aes(x=time, y=cases.median, size=0.2), color="black") +
  geom_point(data=obs, aes(x=time, y=cases, shape=1), size=0.25, color="black") +
  geom_density(
    data=predThresholdsFull,
    aes(x=onset, y=5 * ..density.., fill="black"),
    # width=1,
    color=NA, trim=TRUE
    # trim=FALSE, draw_quantiles=c(0.025, 0.5, 0.975)
  ) +
  geom_density(
    data=predThresholdsFull,
    aes(x=offset, y=5 * ..density.., fill="black"),
    # width=1,
    color=NA, trim=TRUE
    # trim=FALSE, draw_quantiles=c(0.025, 0.5, 0.975)
  ) +
  # coord_cartesian(xlim=c(zoomedStartWeek, zoomedEndWeek), ylim=c(0, 4*seasonThreshold)) +
  scale_shape_identity(guide="legend", breaks=c(1), labels=c("Observed"), name=NULL) +
  scale_fill_identity(guide="legend", breaks=c("black"), labels=c("Onset / offset"), name=NULL) +
  scale_size_identity(guide="legend", breaks=c(0.2), labels=c("Model"), name=NULL) +
  commonOptions

dev.off()

tikz(sprintf("%s/samplePredSingle.tex", figuresDir), width=0.8*pagePlotWidth, height=2*pagePlotHeight, pointsize=10, standAlone = TRUE)

data = predCumSingle

thresholds = data %>%
	mutate(ppx=aapStrat(time)) %>%
	filter(ppx > 0) %>%
	do((function(df) {
	  df = df %>% arrange(time)
		data.frame(
			ppx.start=min(df$time),
			ppx.end=max(df$time),
			unprotected.start=first(df$cases.cumrel),
			unprotected.end=last(df$cases.cumrel)
		)
	})(.)) %>%
	ungroup() %>%
	as.data.frame()

ggplot(data) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(aes(x=time, y=cases.cumrel, group=pspline.sample), color="grey") +
  geom_segment(data=thresholds, aes(x=ppx.start, y=unprotected.start, xend=ppx.start, yend=+Inf), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.end, y=unprotected.start, xend=ppx.end, yend=+Inf), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.start, y=unprotected.start, xend=+Inf, yend=unprotected.start), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.start, y=unprotected.end, xend=+Inf, yend=unprotected.end), linetype="11", size=.375) +
  geom_point(data=obs, aes(x=time, y=cases.cumrel, size=0.5), size=0.5) +
	geom_rect(data=thresholds, aes(xmin=ppx.start, xmax=ppx.end, ymin=1, ymax=1.1), size=.375) +
	geom_text(data=thresholds, aes(x=mean(c(ppx.start, ppx.end)), y=1.05, label="Prophylaxis window"), color="white", size=3.5) +
	geom_rect(data=thresholds, aes(xmin=max(monthBoundaries$max) - 3, xmax=max(monthBoundaries$max), ymin=unprotected.start, ymax=unprotected.end), size=.375) +
	geom_text(data=thresholds, aes(x=max(monthBoundaries$max) - 1.5, y=mean(c(unprotected.start, unprotected.end)), label="$\\Delta$ = preventable fraction", angle=90, lineheight=0.5), color="white", size=3.5) +
  scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
  coord_cartesian(xlim=range(c(monthBoundaries$min, monthBoundaries$max))) +
  labs(x="Time", y="RSV relative cumulative incidence") +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
    legend.position="bottom"
  )

dev.off()

tikz(sprintf("%s/samplePreventableFractionFull.tex", figuresDir), width=0.8*pagePlotWidth, height=2*pagePlotHeight, pointsize=10, standAlone = TRUE)

data = predCumFull
splineData = data %>% filter(pspline.sample < sampleDisplayNsim)

splinePlot = ggplot(data) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(data=splineData, aes(x=time, y=cases.cumrel, group=pspline.sample), color="gray") +
  geom_segment(data=thresholds, aes(x=ppx.start, y=-Inf, xend=ppx.start, yend=+Inf), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.end, y=-Inf, xend=ppx.end, yend=+Inf), linetype="11", size=.375) +
	geom_rect(data=thresholds, aes(xmin=ppx.start, xmax=ppx.end, ymin=1, ymax=1.1), size=.375) +
	geom_text(data=thresholds, aes(x=mean(c(ppx.start, ppx.end)), y=1.05, label="Prophylaxis window"), color="white", size=3.5) +
  geom_point(data=obs, aes(x=time, y=cases.cumrel), size=0.5) +
  scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
  scale_y_continuous() +
  labs(x="Time", y="RSV relative cumulative incidence") +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
    legend.position="bottom"
  )

violinPlot <- ggplot(predFractionFull) +
  theme_light(base_size=plotTextBaseSize) +
  geom_violin(
    aes(x=0, y=preventable),
    width=.01,
    fill="gray50", color="black",
    trim=FALSE, draw_quantiles=c(0.025, 0.5, 0.975)
  ) +
  scale_x_continuous() +
  scale_y_continuous() +
  coord_cartesian(ylim=c(0, 1)) +
  labs(y="RSV preventable fraction") +
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    legend.position="bottom"
  )
grid.arrange(splinePlot, violinPlot, layout_matrix=rbind(c(rep(1, 4), 2)))

dev.off()

