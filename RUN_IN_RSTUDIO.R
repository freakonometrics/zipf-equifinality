# v0.5.1 — run the paper analysis directly from the RStudio console.
# No terminal or Anaconda is required.

source("R/generators.R")
source("R/diagnostics.R")
source("R/calibration.R")
source("R/figures.R")
source("R/run_all.R")

# Full paper analysis:
run_paper_analysis()

# Then open index.qmd and click Render. Because results now exist, the QMD will
# regenerate tables and figures without rerunning the simulations unless
# force_rerun: true is set in the YAML header.
