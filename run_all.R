# Run the whole analysis in order. In RStudio: open the .Rproj, then
# source("run_all.R"). Roughly 1 minute for the one-at-a-time sweeps,
# 15 minutes for the global analysis, 5 minutes for the stochastic check.
source("R/setup.R")
source("R/run_oat.R")
source("R/run_global.R")
source("R/run_stochastic.R")
cat("Done. Figures in output/figs, tables in output/tables.\n")
