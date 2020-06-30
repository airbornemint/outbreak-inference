library(doParallel)
library(parallel)

# Set up parallel processing, used for simulation study
registerDoParallel(min(detectCores(), vignette.eval$limitCores, na.rm=TRUE))
