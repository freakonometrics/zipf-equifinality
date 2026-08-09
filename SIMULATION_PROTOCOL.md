# Simulation protocol — v0.5.1

## Fixed design

- Corpus size: `n = 100000` tokens.
- Primary rank window: ranks 10--500, with empirical count at least 5.
- Held-out paper replications: 100, beginning at seed 2001.
- Within-sequence shuffle surrogates for excess NMI: 20 per replication.
- Primary collapsed alphabet: top 50 frequent types + `OTHER`.
- NMI sensitivity: top 25, 50 and 100 types + `OTHER`.
- Bootstrap rank-effect replicates: 2000.
- Markov persistence: `rho = 0.35`.

## Calibration / evaluation separation

Core calibration models use seeds beginning at 1001. The latent mixture uses a
50-seed final calibration panel, 1001:1050 by default. The broad latent grid is
screened on a smaller subset, after which the retained regions and local
refinements are evaluated on the full 50-seed calibration panel.

No held-out seed beginning at 2001 is used to select calibration parameters.

## Strict benchmark

1. finite Zipf reference;
2. stationary finite-Zipf Markov control;
3. latent mixture;
4. sample-space reduction.

## Compatible / mimicking benchmark

1. random segmentation;
2. Simon reinforcement.

## Primary diagnostics

- OLS Zipf exponent;
- conditional-MLE Zipf exponent on the same empirical rank window;
- final vocabulary size;
- global and early/late Heaps slopes;
- excess lag-1 normalized mutual information;
- transition direction `P(next rank is more frequent | rank changes)`.

## Joint robustness benchmark

The second calibration targets `(OLS alpha, final V) = (1, 5000)` and evaluates
all higher-order diagnostics out of calibration.
