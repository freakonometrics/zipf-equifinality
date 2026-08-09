# v0.5.1 notes

## Prespecified numerical change

The latent-mixture calibration was the only part of v0.5 showing a systematic
difference between a 10-seed calibration average and the 100-seed held-out
average. v0.5.1 therefore stabilizes latent calibration without touching the
held-out panel:

1. broad latent grid screened on a small calibration subset;
2. top regions locally refined;
3. final parameter selected on 50 calibration seeds (1001:1050 by default);
4. held-out evaluation remains seeds 2001:2100.

The same strategy is used for exponent-only and joint `(alpha, V_n)` latent
calibration.

## No changes to the scientific diagnostics

The Markov control, OLS/MLE comparison, bootstrap AUC, joint-matching benchmark,
and NMI sensitivity remain as in v0.5. The goal is calibration stability, not
post-hoc retuning of the benchmark.

## Reproducibility

`index.qmd` is now rerun-safe: on a fresh repository it generates simulations,
tables and figures; on later renders it detects existing paper outputs and
skips simulation unless `force_rerun: true` is set.
