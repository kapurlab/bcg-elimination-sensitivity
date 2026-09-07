# Install the packages the analysis needs. Run once.
pkgs <- c("deSolve", "SimInf", "lhs", "sensitivity", "ggplot2", "dplyr", "tidyr")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("All packages available:", paste(pkgs, collapse = ", "), "\n")
