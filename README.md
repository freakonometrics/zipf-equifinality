# Beyond Zipf's Law — reproducibility companion

This repository contains the computational companion to **Beyond Zipf's Law: Equifinality and Mechanistic Inference from Scaling Laws**.

The landing page is `index.qmd`. It is designed both as a readable HTML companion and as an executable record of the simulations.

## Render the online companion

The default mode uses the paper-facing precomputed outputs bundled with the repository:

```bash
quarto render
```

The website is written to `docs/`, which can be published directly with GitHub Pages.

## Rebuild from existing simulation outputs

If `results/` already contains a full run, change the YAML parameter in `index.qmd` to:

```yaml
analysis_mode: "results"
```

and render again. Figures and tables are then rebuilt from the raw outputs.

## Rerun the full analysis

Set:

```yaml
analysis_mode: "recompute"
```

and render `index.qmd`. For a short smoke test, also set `quick: true`.

Alternatively, from RStudio:

```r
source("RUN_IN_RSTUDIO.R")
```

then render `index.qmd`.

## Main design

The primary benchmark compares:

- i.i.d. finite Zipf sampling;
- a sticky Markov process with the same stationary finite-Zipf marginal;
- canonical sample-space reduction (SSR), also with the same stationary finite-Zipf marginal;
- a latent-scale mixture that obtains a similar marginal through pooling.

Random segmentation and Simon reinforcement are used separately for observation-window sensitivity.

## Required R packages

- ggplot2
- patchwork
- knitr

The simulation source is under `R/`. Precomputed website assets are stored under `figures/precomputed/` and `tables/precomputed/` so the public companion can be rendered without rerunning the Monte Carlo study.
