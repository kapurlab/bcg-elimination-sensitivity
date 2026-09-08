# Run the whole analysis in order. In RStudio: open the .Rproj, then
# source("run_all.R"). Roughly 1 minute for the one-at-a-time sweeps,
# 15 minutes for the global analysis, 5 minutes for the stochastic check,
# about 20 minutes for the revaccination scenarios, a few minutes for onset.
source("R/setup.R")
source("R/run_oat.R")
source("R/run_global.R")
source("R/run_stochastic.R")
source("R/run_revaccination.R")
source("R/run_onset.R")
source("R/run_memo_figure.R")
source("R/run_tiles.R")
cat("Done. Figures in output/figs, tables in output/tables.\n")
