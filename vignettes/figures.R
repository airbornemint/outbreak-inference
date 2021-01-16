library(ggplot2)

# These are utilities used for plots in the vignette

# Correspondence between epi weeks and month labels. 
# Fudge by .01 to let minor gridlines show through
epiWeekBreaks <- c(3.01, 7.25, 11.75, 16.01, 20.25, 24.75, 29.01, 33.01, 37.5, 42.01, 46.25, 50.75)
epiWeekLabels <- c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun")

breaks.df <- data.frame(i=seq(1, 12), mid=epiWeekBreaks)
monthBoundaries <- breaks.df %>%
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

ggplot <- function(data, legends) {
  plotFills = c(
    casesObs=NA, 
    casesEst=grey(0.75), 
    casesEstMedian=grey(0.75), 
    casesEstSamples=grey(0.75), 
    seriousEstMedian=grey(0.5), 
    supplyEnd=grey(0.5)
  )
  plotColors = c(
    casesObs="black", 
    casesEst=grey(0.75), 
    casesEstMedian=grey(0.75), 
    casesEstSamples=grey(0.75), 
    seriousEstMedian=grey(0.5), 
    supplyEnd=grey(0.5)
  )
  plotSizes = c(
    casesObs=0.5, 
    casesEst=1, 
    casesEstMedian=1, 
    casesEstSamples=1, 
    seriousEstMedian=1, 
    supplyEnd=1
  )
  legendLabels = c(
    casesObs="Observed cases", 
    casesEst="Predicted cases (95% CI)", 
    casesEstMedian="Predicted cases (median)", 
    casesEstSamples="Predicted cases (sampled)", 
    seriousEstMedian="Predicted serious cases (median)", 
    supplyEnd="End of supply (sampled)"
  )
  
  ggplot2::ggplot(data) +
    theme_light() +
    scale_y_continuous() +
    scale_x_continuous(breaks=monthBoundaries$mid, labels=epiWeekLabels, expand=c(0, 0)) + 
    labs(x=NULL, y=NULL) + 
    coord_cartesian(xlim=range(c(monthBoundaries$min, monthBoundaries$max))) +
    scale_fill_manual(
      name=NULL, 
      values=plotFills, 
      limits=legends, 
      labels=legendLabels,
      drop=TRUE,
      guide=guide_legend(ncol=1)
    ) +
    scale_colour_manual(
      name=NULL, 
      values=plotColors, 
      limits=legends, 
      labels=legendLabels,
      drop=TRUE,
      guide=guide_legend(ncol=1)
    ) +
    scale_size_manual(
      name=NULL, 
      values=plotSizes, 
      limits=legends, 
      labels=legendLabels,
      drop=TRUE,
      guide=guide_legend(ncol=1)
    ) +
    theme(
      panel.grid.major.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_text(hjust=0.5, vjust=0.5),
      legend.key = element_rect(colour="white", fill="white"),
      legend.position = "bottom"
    )
}

knitr::opts_chunk$set(fig.width=4, fig.height=3, fig.align="center", dev=ifelse(capabilities('cairo'), 'svg', 'png'))
