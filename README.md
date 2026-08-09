# Beyond Zipf's Law — reproducible benchmark v0.5.1

This repository is the reproducibility companion for the working JRSI paper on
Zipf-law equifinality and mechanistic inference.

## What changed in v0.5.1

The benchmark architecture is unchanged from v0.5. The only calibration change
is prespecified: the latent-mixture parameter is selected on a 50-seed
calibration panel rather than 10 seeds, while held-out evaluation remains on
seeds 2001 onward. This addresses the larger Monte Carlo variability of the
latent generator without using held-out results for tuning.

The notebook also produces paper-ready tables and figures for:

- OLS versus conditional-MLE exponent estimates;
- a non-SSR stationary Zipf Markov control;
- bootstrap AUC / signed rank-effect uncertainty;
- joint matching on `(alpha, V_n)`;
- excess-NMI sensitivity for `K = 25, 50, 100`;
- observation-window sensitivity and latent conditioning.

## RStudio: one-click route

Open `index.qmd` and click **Render**.

On a fresh repository, the default YAML settings run the full analysis because
required result files are absent. Later renders reuse `results/` and only
rebuild the notebook, tables and figures. To force a complete rerun, set
`force_rerun: true` for one render.

Alternatively, from the RStudio console:

```r
source("RUN_IN_RSTUDIO.R")
```

and then Render `index.qmd`.

## Required R packages

- ggplot2
- patchwork
- knitr

No terminal workflow is required.

## Main outputs

Paper-facing tables are written to `tables/` as both CSV and LaTeX fragments.
Paper-facing figures are written to `figures/` as PDFs. Raw calibration and
held-out outputs are written to `results/`.
