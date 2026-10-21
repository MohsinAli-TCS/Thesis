# =============================================================================
# run_all.R — reproduces every table and figure in the manuscript.
# Set the working directory to the repository root, then: source("run_all.R")
# Runtime: about 45-60 minutes on a laptop (dominated by 500-rep bootstraps).
# =============================================================================
pkgs <- c("fixest", "dplyr", "splines", "ggplot2", "jsonlite", "scales")
miss <- setdiff(pkgs, rownames(installed.packages())); if (length(miss)) install.packages(miss)
.scripts <- c("01_build_data.R", "02_main_results.R", "03_instrument_designs.R", "04_robustness.R", "05_placebo.R")
for (.k in seq_along(.scripts)) {
  .t0 <- Sys.time(); cat("\n==== running", .scripts[.k], "====\n")
  source(file.path("code", .scripts[.k]), local = new.env(), echo = FALSE)   # each script runs in its own environment
  cat("finished", .scripts[.k], "in", round(as.numeric(Sys.time() - .t0, units = "mins"), 1), "min\n")
}
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
