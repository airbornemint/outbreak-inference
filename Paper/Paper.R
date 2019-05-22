# ---- paper ----
import::from(outbreakinference, outbreak.calc.cumcases, outbreak.calc.cases, outbreak.calc.thresholds, outbreak.estimate.scalars)

import::from(reshape, melt)
import::from(dplyr, "%>%", rename, mutate, bind_rows, select, filter, group_by, ungroup, do, arrange, first, last, inner_join)
import::from(plyr, rlply, llply)
import::from(mgcv, gam, mroot)
import::from(data.table, setnames)
import::from(tidyr, unnest)

source("./Common.R")

sampleObs = read.csv("../vignettes/seasonal.csv")
sampleObs$cases.cum = cumsum(sampleObs$cases)
sampleObs$cases.cum.frac = sampleObs$cases.cum / max(sampleObs$cases.cum)

simulations = 2000
eps = .05
seasonThreshold = 0.025
ppxDuration = 24 # weeks
startYear = 1996
endYear = 2013

obsTime = sampleObs$time
modelTime = seq(min(obsTime) - 1 + eps, max(obsTime), eps)

sampleModel = gam(cases ~ s(time, k=20, bs="cp", m=3), family=poisson, data=sampleObs)

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

sampleParamsSingle = sampleModel %>% outbreakinference:::outbreak.estimate.params(1)
predCumSingle = sampleParamsSingle %>% predCum(sampleModel, modelTime)
predFractionSingle = sampleParamsSingle %>% predFraction(sampleModel, modelTime)

sampleNsimMini = 5
sampleParamsMini = sampleModel %>% outbreakinference:::outbreak.estimate.params(sampleNsimMini)
predCasesMini = sampleParamsMini %>% predCases(sampleModel, modelTime)
predCumMini = sampleParamsMini %>% predCum(sampleModel, modelTime)
predOnsetMini = sampleParamsMini %>% predOnset(sampleModel, modelTime)

zoomedStartWeek = min(monthBoundaries$min) + 5
zoomedEndWeek = zoomedStartWeek + 21

sampleNsimFull = simulations
sampleParamsFull = sampleModel %>% outbreakinference:::outbreak.estimate.params(sampleNsimFull)
predCumFull = sampleParamsFull %>% predCum(sampleModel, modelTime)
predOnsetFull = sampleParamsFull %>% predOnset(sampleModel, modelTime)
predFractionFull = sampleParamsFull %>% predFraction(sampleModel, modelTime)

sampleDisplayNsim = 100

# ---- figures ----

import::from(ggplot2, ggplot, theme_light, geom_point, aes, scale_x_continuous, scale_y_continuous, geom_line, geom_segment, theme, coord_cartesian, element_blank, element_text, geom_rect, geom_text, geom_violin, labs)
import::from(ggstance, geom_violinh)
import::from(gridExtra, grid.arrange)
import::from(grDevices, dev.off)
import::from(tikzDevice, tikz)

inlinePlotWidth = 3.1
inlinePlotHeight = 2.5
pagePlotWidth = inlinePlotWidth * 2.1
pagePlotHeight = inlinePlotHeight * 2.1 * 2 / 5
plotTextBaseSize = 8

figuresDir = "./figures"
if (!dir.exists(figuresDir)) {
  dir.create(figuresDir, recursive=TRUE)
}
options(tikzMetricsDictionary='./tikzDictionary.dat')

tikz(sprintf("%s/sampleCases.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10)

ggplot(sampleObs) +
    theme_light(base_size=plotTextBaseSize) +
    geom_point(aes(x=time, y=cases), size=.5) +
    scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
    scale_y_continuous() +
    coord_cartesian(xlim=range(c(monthBoundaries$min, monthBoundaries$max))) +
    labs(x="Time", y="RSV incidence") +
    theme(
        panel.grid.major.x=element_blank(),
        legend.title=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
        legend.position="bottom"
    )

dev.off()

tikz(sprintf("%s/samplePredMini.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10)

ggplot(predCasesMini) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(aes(x=time, y=cases, group=sim), color="grey") +
  geom_point(data=sampleObs, aes(x=time, y=cases), size=0.5) +
  scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
  scale_y_continuous() +
  coord_cartesian(xlim=range(c(monthBoundaries$min, monthBoundaries$max))) +
  labs(x="Time", y="RSV incidence") +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
    legend.position="bottom"
  )

dev.off()

tikz(sprintf("%s/samplePredCumMini.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10)

ggplot(predCumMini) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(aes(x=time, y=cases.cum.frac, group=sim), color="grey") +
  geom_segment(aes(x=-Inf, y=seasonThreshold, xend=+Inf, yend=seasonThreshold), linetype="11", data=data.frame(), size=.375) +
  geom_point(data=sampleObs, aes(x=time, y=cases.cum.frac), size=0.5) +
  scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
  coord_cartesian(xlim=range(c(monthBoundaries$min, monthBoundaries$max))) +
  labs(x="Time", y="Relative cumulative incidence") +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
    legend.position="bottom"
  )

dev.off()

tikz(sprintf("%s/samplePredOnsetMini.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10)

ggplot(predCumMini) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(aes(x=time, y=cases.cum.frac, group=sim), color="grey") +
  geom_segment(aes(x=-Inf, y=seasonThreshold, xend=+Inf, yend=seasonThreshold), linetype="11", data=data.frame(), size=.375) +
  geom_point(data=sampleObs, aes(x=time, y=cases.cum.frac), size=0.5) +
  geom_point(data=predOnsetMini, aes(x=onset, y=seasonThreshold), size=0.75, shape=17) +
  geom_segment(data=predOnsetMini, aes(x=onset, xend=onset, y=seasonThreshold, yend=0), linetype="11") +
  scale_y_continuous() +
  scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
  coord_cartesian(xlim=c(zoomedStartWeek, zoomedEndWeek), ylim=c(0, 2*seasonThreshold)) +
  labs(x="Time", y="Relative cumulative incidence") +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
    legend.position="bottom"
  )

dev.off()

tikz(sprintf("%s/samplePredOnsetFull.tex", figuresDir), width=pagePlotWidth * 0.4, height=pagePlotHeight, pointsize=10)

data = predCumFull
splineData = data %>% filter(sim < sampleDisplayNsim)

ggplot(predOnsetFull) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(data=splineData, aes(x=time, y=cases.cum.frac, group=sim), color="gray") +
  geom_segment(aes(x=-Inf, y=seasonThreshold, xend=+Inf, yend=seasonThreshold), linetype="11", data=data.frame(), size=.375) +
  geom_point(data=sampleObs, aes(x=time, y=cases.cum.frac), size=0.5) +
  geom_violinh(
    data=predOnsetFull,
    aes(x=onset, y=seasonThreshold),
    width=.01,
    fill="gray50", color="black",
    trim=FALSE, draw_quantiles=c(0.025, 0.5, 0.975)
  ) +
  scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) +
  scale_y_continuous() +
  coord_cartesian(xlim=c(zoomedStartWeek, zoomedEndWeek), ylim=c(0, 4*seasonThreshold)) +
  labs(x="Time", y="Relative cumulative incidence") +
  theme(
    panel.grid.major.x=element_blank(),
    legend.title=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
    legend.position="bottom"
  )

dev.off()

tikz(sprintf("%s/samplePredSingle.tex", figuresDir), width=0.8*pagePlotWidth, height=2*pagePlotHeight, pointsize=10)

data = predCumSingle

thresholds = data %>%
	mutate(ppx=aapStrat(time)) %>%
	filter(ppx > 0) %>%
	do((function(df) {
	  df = df %>% arrange(time)
		data.frame(
			ppx.start=min(df$time),
			ppx.end=max(df$time),
			unprotected.start=first(df$cases.cum.frac),
			unprotected.end=last(df$cases.cum.frac)
		)
	})(.)) %>%
	ungroup() %>%
	as.data.frame()

ggplot(data) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(aes(x=time, y=cases.cum.frac, group=sim), color="grey") +
  geom_segment(data=thresholds, aes(x=ppx.start, y=unprotected.start, xend=ppx.start, yend=+Inf), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.end, y=unprotected.start, xend=ppx.end, yend=+Inf), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.start, y=unprotected.start, xend=+Inf, yend=unprotected.start), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.start, y=unprotected.end, xend=+Inf, yend=unprotected.end), linetype="11", size=.375) +
  geom_point(data=sampleObs, aes(x=time, y=cases.cum.frac, size=0.5), size=0.5) +
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

tikz(sprintf("%s/samplePreventableFractionFull.tex", figuresDir), width=0.8*pagePlotWidth, height=2*pagePlotHeight, pointsize=10)

data = predCumFull
splineData = data %>% filter(sim < sampleDisplayNsim)

splinePlot = ggplot(data) +
  theme_light(base_size=plotTextBaseSize) +
  geom_line(data=splineData, aes(x=time, y=cases.cum.frac, group=sim), color="gray") +
  geom_segment(data=thresholds, aes(x=ppx.start, y=-Inf, xend=ppx.start, yend=+Inf), linetype="11", size=.375) +
  geom_segment(data=thresholds, aes(x=ppx.end, y=-Inf, xend=ppx.end, yend=+Inf), linetype="11", size=.375) +
	geom_rect(data=thresholds, aes(xmin=ppx.start, xmax=ppx.end, ymin=1, ymax=1.1), size=.375) +
	geom_text(data=thresholds, aes(x=mean(c(ppx.start, ppx.end)), y=1.05, label="Prophylaxis window"), color="white", size=3.5) +
  geom_point(data=sampleObs, aes(x=time, y=cases.cum.frac), size=0.5) +
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

