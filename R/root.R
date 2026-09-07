# Every script sources this first. It checks that the working directory is
# the project root (where the .Rproj file lives), which is what RStudio sets
# when the project is open, and moves there if a script was run from R/.
if (!file.exists("Model_Sensitivity_Tradeoffs.Rproj")) {
  if (file.exists("../Model_Sensitivity_Tradeoffs.Rproj")) setwd("..") else
    stop("Set the working directory to the project root (open the .Rproj in RStudio).")
}
dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
