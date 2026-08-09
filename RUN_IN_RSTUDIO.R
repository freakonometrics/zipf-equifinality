# Full analysis for the Beyond Zipf's Law reproducibility companion.
# No terminal is required: source this file from the repository root in RStudio.

source("R/generators.R")
source("R/diagnostics.R")
source("R/calibration.R")
source("R/figures.R")
source("R/run_all.R")

# Optional short check:
# run_smoke_test()

# Full calibration, held-out benchmark and robustness analyses:
run_paper_analysis()

# Then render index.qmd with analysis_mode: "results" to rebuild the HTML page
# from the raw simulation outputs without rerunning the Monte Carlo study.
